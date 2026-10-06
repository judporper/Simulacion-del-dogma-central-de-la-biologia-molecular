from __future__ import annotations

from dataclasses import dataclass

from Bio.Data import IUPACData
from Bio.Seq import Seq

TAM_CEBADOR = 6

NOMBRE_AMINOACIDO = {
    "A": "Alanina", "R": "Arginina", "N": "Asparagina", "D": "Ác. aspártico",
    "C": "Cisteína", "Q": "Glutamina", "E": "Ác. glutámico", "G": "Glicina",
    "H": "Histidina", "I": "Isoleucina", "L": "Leucina", "K": "Lisina",
    "M": "Metionina", "F": "Fenilalanina", "P": "Prolina", "S": "Serina",
    "T": "Treonina", "W": "Triptófano", "Y": "Tirosina", "V": "Valina",
}


@dataclass
class Okazaki:
    indice: int
    cebador_arn: str
    fragmento_adn: str
    inicio_molde: int
    fin_molde: int

    @property
    def longitud(self) -> int:
        return len(self.cebador_arn) + len(self.fragmento_adn)


@dataclass
class Codon:
    posicion: int
    codon: str
    aminoacido: str

    @property
    def anticodon(self) -> str:
        return str(Seq(self.codon).complement_rna())


@dataclass
class Simulacion:
    nombre: str
    descripcion: str
    codificante: str
    molde: str
    cebador_lider: str
    lider: str
    okazaki: list[Okazaki]
    rezagada: str
    arnm: str
    inicio: int
    codones: list[Codon]
    proteina: str
    parada: bool


def replicar(codificante: Seq, tam_okazaki: int = 10) -> dict:
    if tam_okazaki <= TAM_CEBADOR:
        raise ValueError(f"El tamaño de fragmento de Okazaki debe ser mayor que {TAM_CEBADOR} nt.")
    molde = codificante.complement()
    n = len(codificante)

    lider = molde.complement()

    trozos = [[i, min(i + tam_okazaki, n)] for i in range(0, n, tam_okazaki)]
    if len(trozos) > 1 and trozos[-1][1] - trozos[-1][0] <= TAM_CEBADOR:
        ultimo = trozos.pop()
        trozos[-1][1] = ultimo[1]

    okazaki = []
    for k, (ini, fin) in enumerate(trozos, 1):
        nuevo = codificante[ini:fin].reverse_complement()
        okazaki.append(Okazaki(k, str(nuevo[:TAM_CEBADOR].transcribe()),
                               str(nuevo[TAM_CEBADOR:]), ini, fin))

    return dict(cebador_lider=str(lider[:TAM_CEBADOR].transcribe()), lider=str(lider),
                okazaki=okazaki, rezagada=str(codificante.complement()))


def transcribir(molde: Seq) -> Seq:
    return molde.complement_rna()


def traducir(arnm: Seq) -> dict:
    inicio = arnm.find("AUG")
    if inicio == -1:
        return dict(inicio=-1, codones=[], proteina="", parada=False)

    orf = arnm[inicio:]
    orf = orf[: len(orf) // 3 * 3]
    proteina, parada, _ = str(orf.translate()).partition("*")
    codones = [Codon(inicio + 3 * i, str(orf[3 * i: 3 * i + 3]), aa)
               for i, aa in enumerate(proteina + parada)]
    return dict(inicio=inicio, codones=codones, proteina=proteina, parada=bool(parada))


def simular(secuencia: str, nombre: str = "Secuencia sin nombre",
            descripcion: str = "", tam_okazaki: int = 10) -> Simulacion:
    secuencia = "".join(secuencia.split()).upper()
    if not secuencia or set(secuencia) - set(IUPACData.unambiguous_dna_letters):
        raise ValueError("La secuencia de ADN solo puede contener las bases A, C, G, T.")

    codificante = Seq(secuencia)
    molde = codificante.complement()
    arnm = transcribir(molde)
    return Simulacion(nombre, descripcion, secuencia, str(molde),
                      **replicar(codificante, tam_okazaki),
                      arnm=str(arnm), **traducir(arnm))
