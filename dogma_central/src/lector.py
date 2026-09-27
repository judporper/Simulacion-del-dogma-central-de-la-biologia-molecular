"""
lector.py
---------
Funciones para obtener una secuencia de ADN con la que alimentar al
simulador:

    - leer_fasta(ruta)     -> lee el primer registro de un fichero .fasta
    - leer_genbank(ruta)   -> lee la sección ORIGIN de un fichero .gb / .genbank
    - generar_adn_aleatorio(longitud) -> secuencia aleatoria (para pruebas rápidas)

No se depende de Biopython: los parsers son deliberadamente sencillos
(formato de texto plano) para que el proyecto funcione con Python puro,
pero son suficientes para los ficheros de ejemplo de la práctica
(FASTA y GenBank descargados de NCBI).
"""

from __future__ import annotations

import random
import re
from dataclasses import dataclass

BASES = "ACGT"


@dataclass
class SecuenciaLeida:
    identificador: str
    descripcion: str
    secuencia: str
    fuente: str  # "fasta", "genbank" o "aleatoria"


def generar_adn_aleatorio(longitud: int = 60, semilla: int | None = None) -> SecuenciaLeida:
    """Genera una secuencia de ADN aleatoria de la longitud indicada."""
    if longitud < 3:
        raise ValueError("La longitud mínima razonable es 3 (un codón).")
    rng = random.Random(semilla)
    secuencia = "".join(rng.choice(BASES) for _ in range(longitud))
    return SecuenciaLeida(
        identificador="ADN_ALEATORIO",
        descripcion=f"Secuencia aleatoria generada ({longitud} pb)",
        secuencia=secuencia,
        fuente="aleatoria",
    )


def leer_fasta(ruta: str) -> SecuenciaLeida:
    """
    Lee el PRIMER registro de un fichero FASTA (.fasta / .fa).

    Formato esperado:
        >identificador descripción opcional
        ACGT...
        ACGT...
    """
    identificador = ""
    descripcion = ""
    lineas_secuencia: list[str] = []
    encabezado_leido = False

    with open(ruta, "r", encoding="utf-8") as f:
        for linea in f:
            linea = linea.rstrip("\n")
            if not linea:
                continue
            if linea.startswith(">"):
                if encabezado_leido:
                    break  # ya teníamos un registro completo: nos quedamos con el primero
                encabezado_leido = True
                cabecera = linea[1:].strip()
                partes = cabecera.split(maxsplit=1)
                identificador = partes[0] if partes else "SIN_ID"
                descripcion = partes[1] if len(partes) > 1 else ""
            else:
                if encabezado_leido:
                    lineas_secuencia.append(linea.strip())

    if not encabezado_leido:
        raise ValueError(f"El fichero '{ruta}' no parece un FASTA válido (falta '>').")

    secuencia = "".join(lineas_secuencia).upper()
    secuencia = re.sub(r"[^ACGT]", "", secuencia)  # descarta N u otros códigos IUPAC ambiguos
    if not secuencia:
        raise ValueError(f"No se ha podido extraer una secuencia de ADN válida de '{ruta}'.")

    return SecuenciaLeida(
        identificador=identificador,
        descripcion=descripcion,
        secuencia=secuencia,
        fuente="fasta",
    )


def leer_genbank(ruta: str) -> SecuenciaLeida:
    """
    Lee un fichero GenBank (.gb / .genbank) de forma sencilla:
    extrae LOCUS/DEFINITION como cabecera y la sección ORIGIN como secuencia.
    """
    identificador = ""
    descripcion_partes: list[str] = []
    en_definicion = False
    en_origen = False
    lineas_secuencia: list[str] = []

    with open(ruta, "r", encoding="utf-8") as f:
        for linea in f:
            cruda = linea.rstrip("\n")

            if cruda.startswith("LOCUS"):
                campos = cruda.split()
                if len(campos) > 1:
                    identificador = campos[1]
                continue

            if cruda.startswith("DEFINITION"):
                en_definicion = True
                descripcion_partes.append(cruda[len("DEFINITION"):].strip())
                continue
            if en_definicion:
                if cruda.startswith(" "):
                    descripcion_partes.append(cruda.strip())
                    continue
                en_definicion = False

            if cruda.startswith("ORIGIN"):
                en_origen = True
                continue

            if cruda.strip() == "//":
                en_origen = False
                continue

            if en_origen:
                # Cada línea tiene el formato: "   1 atgacc atgatt acgg..."
                sin_numeros = re.sub(r"^\s*\d+\s*", "", cruda)
                lineas_secuencia.append(sin_numeros)

    if not lineas_secuencia:
        raise ValueError(f"No se ha encontrado sección ORIGIN en '{ruta}'.")

    secuencia = "".join(lineas_secuencia).upper()
    secuencia = re.sub(r"[^ACGT]", "", secuencia)
    if not secuencia:
        raise ValueError(f"No se ha podido extraer una secuencia de ADN válida de '{ruta}'.")

    return SecuenciaLeida(
        identificador=identificador or "SIN_ID",
        descripcion=" ".join(descripcion_partes),
        secuencia=secuencia,
        fuente="genbank",
    )


def leer_secuencia_automatico(ruta: str) -> SecuenciaLeida:
    """Detecta la extensión del fichero y llama al parser adecuado."""
    ruta_lower = ruta.lower()
    if ruta_lower.endswith((".fasta", ".fa", ".fna")):
        return leer_fasta(ruta)
    if ruta_lower.endswith((".gb", ".genbank", ".gbk")):
        return leer_genbank(ruta)
    raise ValueError(
        "Extensión no reconocida. Usa un fichero .fasta/.fa o .gb/.genbank/.gbk"
    )
