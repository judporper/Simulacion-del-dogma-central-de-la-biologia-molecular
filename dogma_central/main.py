"""
main.py
-------
Punto de entrada del simulador del dogma central de la biología molecular.

Ejecutar con:
    python main.py

Ofrece un menú interactivo por consola para:
    1) Escribir una secuencia de ADN manualmente.
    2) Generar una secuencia de ADN aleatoria.
    3) Cargar una secuencia real desde un fichero FASTA (.fasta/.fa).
    4) Cargar una secuencia real desde un fichero GenBank (.gb/.gbk).
    5) Salir.

Tras elegir la secuencia, se ejecuta la simulación completa
(replicación, transcripción, traducción). El resultado no se imprime
por pantalla: se ofrece guardarlo en un informe de texto y, opcionalmente,
en un informe HTML interactivo.
"""

from __future__ import annotations

import os
import sys
import re

from src.lector import (
    SecuenciaLeida,
    generar_adn_aleatorio,
    leer_secuencia_automatico,
)
from src.simulador import SimuladorDogmaCentral
from src import visualizador as vis
from src.informe_html import generar_informe_html

CARPETA_SALIDAS = os.path.join(os.path.dirname(__file__), "salidas")
CARPETA_EJEMPLOS = os.path.join(os.path.dirname(__file__), "data", "ejemplos")


def _nombre_fichero_seguro(identificador: str) -> str:
    """
    Convierte un identificador de secuencia (que puede traer ':', '/',
    espacios u otros símbolos, p.ej. 'NC_000913.3:c366305-363231') en un
    nombre de fichero válido en Windows/macOS/Linux.

    En Windows, un ':' en un nombre de fichero no da error: crea un
    "alternate data stream" y el fichero visible se queda VACÍO, que es
    justo el síntoma de "el informe se genera vacío".
    """
    seguro = re.sub(r"[^A-Za-z0-9._-]+", "_", identificador.strip())
    seguro = seguro.strip("_.")
    return seguro or "secuencia"


def _ruta_disponible(carpeta: str, nombre_base: str, extension: str) -> str:
    """
    Devuelve una ruta libre en 'carpeta' con sufijo numérico, empezando
    siempre en '-1' (nunca sin número), y subiendo a '-2', '-3', ... si
    ya existen ficheros anteriores con ese mismo nombre base (p.ej. varias
    secuencias ALEATORIAS o MANUALES seguidas), para no sobrescribirlos.
    """
    contador = 1
    while True:
        ruta = os.path.join(carpeta, f"{nombre_base}-{contador}{extension}")
        if not os.path.exists(ruta):
            return ruta
        contador += 1


def _pedir_secuencia_manual() -> SecuenciaLeida:
    print("\nEscribe la secuencia de ADN (solo bases A, C, G, T), por ejemplo:")
    print("  ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    seq = input("Secuencia: ").strip().upper()
    return SecuenciaLeida(
        identificador="ADN_MANUAL",
        descripcion="Secuencia introducida manualmente por el usuario",
        secuencia=seq,
        fuente="manual",
    )


def _pedir_secuencia_aleatoria() -> SecuenciaLeida:
    try:
        longitud = int(input("Longitud deseada (por defecto 60): ") or "60")
    except ValueError:
        longitud = 60
    return generar_adn_aleatorio(longitud=longitud)


def _pedir_fichero(carpeta_sugerida: str) -> SecuenciaLeida:
    print(f"\nEjemplos disponibles en '{carpeta_sugerida}':")
    if os.path.isdir(carpeta_sugerida):
        for f in sorted(os.listdir(carpeta_sugerida)):
            print(f"   - {f}")
    ruta = input("Ruta del fichero (.fasta/.fa o .gb/.gbk): ").strip()
    if not os.path.isabs(ruta) and not os.path.exists(ruta):
        ruta_alternativa = os.path.join(carpeta_sugerida, ruta)
        if os.path.exists(ruta_alternativa):
            ruta = ruta_alternativa
    return leer_secuencia_automatico(ruta)


def _menu() -> str:
    print("\n" + "-" * 60)
    print(" SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR")
    print("-" * 60)
    print(" 1) Escribir una secuencia de ADN manualmente")
    print(" 2) Generar una secuencia de ADN aleatoria")
    print(" 3) Cargar una secuencia real desde un fichero FASTA")
    print(" 4) Cargar una secuencia real desde un fichero GenBank")
    print(" 5) Salir")
    return input("Elige una opción [1-5]: ").strip()


def main() -> None:
    os.makedirs(CARPETA_SALIDAS, exist_ok=True)

    while True:
        opcion = _menu()

        if opcion == "5":
            print("¡Hasta luego!")
            sys.exit(0)

        try:
            if opcion == "1":
                datos = _pedir_secuencia_manual()
            elif opcion == "2":
                datos = _pedir_secuencia_aleatoria()
            elif opcion == "3":
                datos = _pedir_fichero(CARPETA_EJEMPLOS)
            elif opcion == "4":
                datos = _pedir_fichero(CARPETA_EJEMPLOS)
            else:
                print("Opción no válida, inténtalo de nuevo.")
                continue
        except (ValueError, FileNotFoundError, OSError) as err:
            print(f"\nError al obtener la secuencia: {err}")
            continue

        try:
            simulador = SimuladorDogmaCentral(datos.secuencia, nombre=datos.identificador)
        except ValueError as err:
            print(f"\nSecuencia no válida: {err}")
            continue

        resultado = simulador.ejecutar()

        print(
            f"\nSimulación completada - ADN: {len(simulador.adn)} pb, "
            f"ARNm: {len(resultado.transcripcion.arnm)} nt, "
            f"Proteína: {len(resultado.traduccion.proteina)} aa"
        )

        respuesta = input("\n¿Guardar informe en un fichero de texto? [S/n]: ").strip().lower()
        if respuesta != "n":
            nombre_base = f"informe_{_nombre_fichero_seguro(datos.identificador)}"
            ruta_txt = _ruta_disponible(CARPETA_SALIDAS, nombre_base, ".txt")
            lineas = vis.construir_informe_texto(
                nombre=datos.identificador,
                descripcion=datos.descripcion,
                longitud_adn=len(simulador.adn),
                replicacion=resultado.replicacion,
                transcripcion=resultado.transcripcion,
                traduccion=resultado.traduccion,
            )
            try:
                vis.guardar_informe_texto(ruta_txt, lineas)
                print(f"Informe guardado en: {ruta_txt}")
            except OSError as err:
                print(f"No se pudo guardar el informe de texto: {err}")

        respuesta_html = input(
            "¿Generar un informe HTML interactivo? [S/n]: "
        ).strip().lower()
        if respuesta_html != "n":
            nombre_base_html = f"informe_{_nombre_fichero_seguro(datos.identificador)}"
            ruta_html = _ruta_disponible(CARPETA_SALIDAS, nombre_base_html, ".html")
            try:
                generar_informe_html(
                    ruta_html,
                    nombre=datos.identificador,
                    descripcion=datos.descripcion,
                    replicacion=resultado.replicacion,
                    transcripcion=resultado.transcripcion,
                    traduccion=resultado.traduccion,
                )
                print(f"Informe interactivo guardado en: {ruta_html}")
                print("   Ábrelo con cualquier navegador para verlo con pestañas e interactividad.")
            except OSError as err:
                print(f"No se pudo guardar el informe HTML: {err}")

        input("\nPulsa Enter para volver al menú...")


if __name__ == "__main__":
    main()