"""
adn.py
------
Representa una molécula de ADN de doble hebra y simula su REPLICACIÓN.

Conceptos que se modelan (según el guion de la práctica):
    - Apertura de la doble hélice (helicasa).
    - Cebadores de ARN (primasa) y ADN polimerasa.
    - Hebra líder: síntesis continua.
    - Hebra rezagada: síntesis discontinua mediante fragmentos de Okazaki,
      unidos después por la ADN ligasa.

Convenio de orientación usado en todo el proyecto:
    - `hebra_codificante` se guarda siempre en sentido 5' -> 3' (izquierda a
      derecha), tal y como se introduce o se lee de un fichero FASTA/GenBank.
    - `hebra_molde` es su complementaria, alineada base a base con la
      anterior, por lo que queda en sentido 3' <- 5' (antiparalela).
"""

from __future__ import annotations

from dataclasses import dataclass, field

BASES_ADN = set("ACGT")

# Reglas de complementariedad de Watson-Crick para el ADN
COMPLEMENTO_ADN = {"A": "T", "T": "A", "C": "G", "G": "C"}


def es_adn_valido(secuencia: str) -> bool:
    """Comprueba que la secuencia contenga únicamente A, C, G, T."""
    return len(secuencia) > 0 and set(secuencia.upper()) <= BASES_ADN


def complementaria_adn(secuencia: str) -> str:
    """Devuelve la hebra complementaria (SIN invertir), base a base."""
    return "".join(COMPLEMENTO_ADN[b] for b in secuencia.upper())


@dataclass
class FragmentoOkazaki:
    """Un fragmento de la hebra rezagada, con su propio cebador de ARN."""
    indice: int
    cebador_arn: str
    fragmento_adn: str

    @property
    def longitud(self) -> int:
        return len(self.fragmento_adn)


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
    def molecula_hija_1(self) -> str:
        """Hebra parental codificante + hebra líder nueva (complementaria)."""
        return self.hebra_lider_nueva

    @property
    def molecula_hija_2(self) -> str:
        """Hebra parental molde + hebra rezagada nueva (complementaria)."""
        return self.hebra_rezagada_nueva


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
        Simula la replicación semiconservativa de la molécula de ADN.

        1) La helicasa "abre" la doble hélice (ya la tenemos separada en
           hebra_codificante / hebra_molde).
        2) La primasa coloca un cebador de ARN al inicio de cada hebra nueva.
        3) La ADN polimerasa sintetiza:
             - la HEBRA LÍDER de forma continua, usando `hebra_codificante`
               como molde (avanza en el mismo sentido que la horquilla).
             - la HEBRA REZAGADA de forma discontinua, usando `hebra_molde`
               como molde, generando fragmentos de Okazaki que luego una
               ADN ligasa uniría en una única hebra continua.
        """
        if tamano_fragmento_okazaki < 1:
            raise ValueError("El tamaño de fragmento de Okazaki debe ser >= 1")

        # --- Hebra líder: síntesis continua, un único cebador ---
        cebador_lider = "AUAAGC"[: min(6, len(self))]  # cebador corto simbólico
        hebra_lider_nueva = complementaria_adn(self.hebra_codificante)

        # --- Hebra rezagada: fragmentos de Okazaki ---
        fragmentos: list[FragmentoOkazaki] = []
        molde_rezagada = self.hebra_molde
        n = len(molde_rezagada)
        indice_fragmento = 1
        for inicio in range(0, n, tamano_fragmento_okazaki):
            trozo_molde = molde_rezagada[inicio: inicio + tamano_fragmento_okazaki]
            cebador = "AUAAGC"[: min(6, len(trozo_molde))]
            trozo_nuevo = complementaria_adn(trozo_molde)
            fragmentos.append(FragmentoOkazaki(indice_fragmento, cebador, trozo_nuevo))
            indice_fragmento += 1

        # La ADN ligasa retira los cebadores de ARN (simplificado aquí: no se
        # simula su eliminación nucleotídica) y une los fragmentos entre sí.
        hebra_rezagada_completa = "".join(f.fragmento_adn for f in fragmentos)

        return ResultadoReplicacion(
            hebra_parental_codificante=self.hebra_codificante,
            hebra_parental_molde=self.hebra_molde,
            cebador_hebra_lider=cebador_lider,
            hebra_lider_nueva=hebra_lider_nueva,
            fragmentos_okazaki=fragmentos,
            hebra_rezagada_nueva=hebra_rezagada_completa,
        )
