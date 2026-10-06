from __future__ import annotations

import os
import random
import re
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Data import CodonTable

TABLA = CodonTable.unambiguous_dna_by_id[1]
CODONES_SENSE = sorted(TABLA.forward_table)
CODONES_PARADA = sorted(TABLA.stop_codons)

FORMATOS = {".fasta": "fasta", ".fa": "fasta", ".fna": "fasta",
            ".gb": "genbank", ".gbk": "genbank", ".genbank": "genbank"}


@dataclass
class SecuenciaLeida:
    identificador: str
    descripcion: str
    secuencia: str


def generar_adn_aleatorio(longitud: int = 60, semilla: int | None = None) -> SecuenciaLeida:
    if longitud < 6:
        raise ValueError("La longitud mínima es 6 pb (codón de inicio ATG + codón de parada).")
    rng = random.Random(semilla)
    cuerpo = "".join(rng.choice(CODONES_SENSE) for _ in range(-(-longitud // 3) - 2))
    secuencia = "ATG" + cuerpo + rng.choice(CODONES_PARADA)
    return SecuenciaLeida(
        "ADN_ALEATORIO",
        f"Gen aleatorio generado ({len(secuencia)} pb): ATG + codones sin parada + codón de parada",
        secuencia,
    )


def leer_secuencia(ruta: str) -> SecuenciaLeida:
    formato = FORMATOS.get(os.path.splitext(ruta)[1].lower())
    if formato is None:
        raise ValueError("Extensión no reconocida. Usa un fichero .fasta/.fa o .gb/.genbank/.gbk")

    with open(ruta, "r", encoding="utf-8") as f:
        registro = next(SeqIO.parse(f, formato), None)
    if registro is None:
        raise ValueError(f"El fichero '{ruta}' no parece un {formato.upper()} válido (sin registros).")

    secuencia = re.sub(r"[^ACGT]", "", str(registro.seq).upper())
    if not secuencia:
        raise ValueError(f"No se ha podido extraer una secuencia de ADN válida de '{ruta}'.")

    descripcion = registro.description.removeprefix(registro.id).strip()
    return SecuenciaLeida(registro.name or "SIN_ID", descripcion, secuencia)
