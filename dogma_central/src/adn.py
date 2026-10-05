"""
adn.py
------
Representa una molécula de ADN de doble hebra y simula su REPLICACIÓN.

Conceptos que se modelan (según el guion de la práctica):
    - Apertura de la doble hélice (helicasa; topoisomerasa y proteínas SSB).
    - Cebadores de ARN (primasa) y ADN polimerasa III.
    - Hebra líder: síntesis continua.
    - Hebra rezagada: síntesis discontinua mediante fragmentos de Okazaki,
      unidos después por la ADN ligasa.

Convenio de orientación usado en todo el proyecto:
    - `hebra_codificante` se guarda siempre en sentido 5' -> 3' (izquierda a
      derecha), tal y como se introduce o se lee de un fichero FASTA/GenBank.
    - `hebra_molde` es su complementaria, alineada base a base con la
      anterior, por lo que queda en sentido 3' <- 5' (antiparalela).
    - La horquilla de replicación avanza de izquierda a derecha.
"""

from __future__ import annotations

from dataclasses import dataclass, field

BASES_ADN = set("ACGT")

# Reglas de complementariedad de Watson-Crick para el ADN
COMPLEMENTO_ADN = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Longitud (nt) de los cebadores de ARN. Es una simplificación didáctica.
TAM_CEBADOR = 6


def es_adn_valido(secuencia: str) -> bool:
    """Comprueba que la secuencia contenga únicamente A, C, G, T."""
    return len(secuencia) > 0 and set(secuencia.upper()) <= BASES_ADN


def complementaria_adn(secuencia: str) -> str:
    """Devuelve la hebra complementaria (SIN invertir), base a base."""
    return "".join(COMPLEMENTO_ADN[b] for b in secuencia.upper())


def revcomp_adn(secuencia: str) -> str:
    """Complementaria inversa: la hebra antiparalela leída en sentido 5' -> 3'."""
    return complementaria_adn(secuencia)[::-1]


def _cebador_arn(adn_5_3: str) -> str:
    """Cebador de ARN: mismos nucleótidos que el ADN al que sustituye, con U en vez de T."""
    return adn_5_3.replace("T", "U")


@dataclass
class FragmentoOkazaki:
    """Fragmento de la hebra rezagada, en sentido 5'->3': cebador de ARN + ADN."""
    indice: int
    cebador_arn: str
    fragmento_adn: str
    inicio_molde: int = 0   # bloque de la hebra molde (0-indexado) que copia
    fin_molde: int = 0

    @property
    def longitud(self) -> int:
        return len(self.cebador_arn) + len(self.fragmento_adn)


@dataclass
class ResultadoReplicacion:
    """Resultado completo de simular la replicación de una molécula de ADN."""
    hebra_parental_codificante: str
    hebra_parental_molde: str
    cebador_hebra_lider: str
    hebra_lider_nueva: str
    fragmentos_okazaki: list[FragmentoOkazaki] = field(default_factory=list)
    hebra_rezagada_nueva: str = ""

    @property
    def molecula_hija_1(self) -> tuple[str, str]:
        """(hebra molde parental, hebra líder nueva)."""
        return self.hebra_parental_molde, self.hebra_lider_nueva

    @property
    def molecula_hija_2(self) -> tuple[str, str]:
        """(hebra codificante parental, hebra rezagada nueva)."""
        return self.hebra_parental_codificante, self.hebra_rezagada_nueva


class ADN:
    """Molécula de ADN de doble hebra."""

    def __init__(self, secuencia_codificante: str, nombre: str = "ADN"):
        secuencia_codificante = secuencia_codificante.upper().replace(" ", "").replace("\n", "")
        if not es_adn_valido(secuencia_codificante):
            raise ValueError(
                "La secuencia de ADN solo puede contener las bases A, C, G, T."
            )
        self.nombre = nombre
        self.hebra_codificante = secuencia_codificante  # 5' -> 3'
        self.hebra_molde = complementaria_adn(secuencia_codificante)  # 3' <- 5'

    def __len__(self) -> int:
        return len(self.hebra_codificante)

    def __repr__(self) -> str:
        return f"ADN('{self.nombre}', {len(self)} pb)"

    # ------------------------------------------------------------------ #
    #  REPLICACIÓN
    # ------------------------------------------------------------------ #
    def replicar(self, tamano_fragmento_okazaki: int = 10) -> ResultadoReplicacion:
        """
        Simula la replicación semiconservativa con la horquilla avanzando de
        izquierda a derecha.

        1) La helicasa abre la doble hélice (hebra_codificante / hebra_molde);
           la topoisomerasa alivia la tensión y las SSB estabilizan las hebras.
        2) HEBRA LÍDER: molde = hebra_molde (se lee 3'->5' en el sentido de la
           horquilla). La ADN polimerasa III sintetiza de forma continua, con
           un único cebador de ARN puesto por la primasa.
        3) HEBRA REZAGADA: molde = hebra_codificante (5'->3' en el sentido de
           la horquilla). Cada fragmento de Okazaki se sintetiza 5'->3' en
           sentido contrario a la horquilla, con su propio cebador de ARN en
           el extremo 5'.
        4) La ADN pol I sustituye los cebadores por ADN y la ADN ligasa une
           los fragmentos en una hebra continua.
        """
        if tamano_fragmento_okazaki <= TAM_CEBADOR:
            raise ValueError(
                f"El tamaño de fragmento de Okazaki debe ser mayor que {TAM_CEBADOR} nt."
            )

        n = len(self)

        # --- Hebra líder: síntesis continua, un único cebador ---
        # Complementaria del molde: coincide con la codificante, en 5'->3'.
        lider = complementaria_adn(self.hebra_molde)
        cebador_lider = _cebador_arn(lider[: min(TAM_CEBADOR, n)])

        # --- Hebra rezagada: fragmentos de Okazaki ---
        # Se parte el molde en bloques; si el último queda con TAM_CEBADOR nt o
        # menos (solo cabría el cebador, sin ADN), se une al bloque anterior.
        limites: list[list[int]] = []
        for inicio in range(0, n, tamano_fragmento_okazaki):
            limites.append([inicio, min(inicio + tamano_fragmento_okazaki, n)])
        if len(limites) > 1 and limites[-1][1] - limites[-1][0] <= TAM_CEBADOR:
            ultimo = limites.pop()
            limites[-1][1] = ultimo[1]

        fragmentos: list[FragmentoOkazaki] = []
        for k, (inicio, fin) in enumerate(limites, start=1):
            bloque = self.hebra_codificante[inicio:fin]
            nuevo = revcomp_adn(bloque)  # fragmento nuevo en sentido 5'->3'
            pc = min(TAM_CEBADOR, len(nuevo))
            fragmentos.append(
                FragmentoOkazaki(
                    indice=k,
                    cebador_arn=_cebador_arn(nuevo[:pc]),
                    fragmento_adn=nuevo[pc:],
                    inicio_molde=inicio,
                    fin_molde=fin,
                )
            )

        # Tras retirar cebadores y unir con la ligasa, la hebra rezagada queda
        # completa. Se muestra alineada con su molde (3'<-5' de izda. a dcha.).
        hebra_rezagada_completa = complementaria_adn(self.hebra_codificante)

        return ResultadoReplicacion(
            hebra_parental_codificante=self.hebra_codificante,
            hebra_parental_molde=self.hebra_molde,
            cebador_hebra_lider=cebador_lider,
            hebra_lider_nueva=lider,
            fragmentos_okazaki=fragmentos,
            hebra_rezagada_nueva=hebra_rezagada_completa,
        )