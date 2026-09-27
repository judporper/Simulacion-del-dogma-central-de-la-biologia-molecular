"""
proteina.py
-----------
Simula la TRADUCCIÓN: ARNm -> secuencia de aminoácidos (proteína).

    - El ribosoma lee el ARNm en tripletes (codones) a partir del primer
      codón de inicio AUG.
    - Cada codón se traduce a un aminoácido según el código genético
      estándar (tabla 1 de la IUPAC/NCBI).
    - La lectura continúa hasta encontrar un codón de parada
      (UAA, UAG, UGA) o hasta que no queden suficientes bases para
      formar otro codón completo.
"""

from __future__ import annotations

from dataclasses import dataclass, field

CODON_STOP = "*"

# Código genético estándar: codón de ARNm -> aminoácido (letra única IUPAC)
TABLA_CODONES: dict[str, str] = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": CODON_STOP, "UAG": CODON_STOP,
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": CODON_STOP, "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Nombres completos de los aminoácidos, para un informe más didáctico
NOMBRE_AMINOACIDO = {
    "A": "Alanina", "R": "Arginina", "N": "Asparagina", "D": "Ác. aspártico",
    "C": "Cisteína", "Q": "Glutamina", "E": "Ác. glutámico", "G": "Glicina",
    "H": "Histidina", "I": "Isoleucina", "L": "Leucina", "K": "Lisina",
    "M": "Metionina", "F": "Fenilalanina", "P": "Prolina", "S": "Serina",
    "T": "Treonina", "W": "Triptófano", "Y": "Tirosina", "V": "Valina",
}

CODON_INICIO = "AUG"
CODONES_PARADA = {"UAA", "UAG", "UGA"}


@dataclass
class CodonTraducido:
    posicion: int          # posición del codón dentro del ARNm (0-indexado)
    codon: str
    aminoacido: str        # letra IUPAC, o "*" si es parada


@dataclass
class ResultadoTraduccion:
    arnm_usado: str
    posicion_inicio: int
    codones: list[CodonTraducido] = field(default_factory=list)
    proteina: str = ""             # secuencia de aminoácidos, letra única
    parada_encontrada: bool = False

    @property
    def proteina_nombres(self) -> list[str]:
        return [NOMBRE_AMINOACIDO.get(aa, aa) for aa in self.proteina]


def traducir(arnm_secuencia: str) -> ResultadoTraduccion:
    """
    Traduce una secuencia de ARNm a proteína, comenzando en el primer
    codón AUG encontrado y terminando en el primer codón de parada.
    """
    arnm_secuencia = arnm_secuencia.upper()
    inicio = arnm_secuencia.find(CODON_INICIO)
    if inicio == -1:
        return ResultadoTraduccion(
            arnm_usado=arnm_secuencia, posicion_inicio=-1, parada_encontrada=False
        )

    codones: list[CodonTraducido] = []
    proteina_letras: list[str] = []
    parada = False

    pos = inicio
    while pos + 3 <= len(arnm_secuencia):
        codon = arnm_secuencia[pos: pos + 3]
        aa = TABLA_CODONES.get(codon)
        if aa is None:
            break  # base ambigua o codón inválido: detenemos la lectura
        codones.append(CodonTraducido(pos, codon, aa))
        if aa == CODON_STOP:
            parada = True
            break
        proteina_letras.append(aa)
        pos += 3

    return ResultadoTraduccion(
        arnm_usado=arnm_secuencia,
        posicion_inicio=inicio,
        codones=codones,
        proteina="".join(proteina_letras),
        parada_encontrada=parada,
    )
