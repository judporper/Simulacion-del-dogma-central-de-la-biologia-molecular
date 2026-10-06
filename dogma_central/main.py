from __future__ import annotations

import os
import re
import sys

from src.dogma import simular
from src.informe_html import generar_informe_html
from src.lector import SecuenciaLeida, generar_adn_aleatorio, leer_secuencia
from src.visualizador import construir_informe_texto, guardar_informe_texto

AQUI = os.path.dirname(__file__)
CARPETA_SALIDAS = os.path.join(AQUI, "salidas")
CARPETA_EJEMPLOS = os.path.join(AQUI, "data", "ejemplos")


def _ruta_libre(identificador: str, extension: str) -> str:
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", identificador.strip()).strip("_.") or "secuencia"
    k = 1
    while os.path.exists(ruta := os.path.join(CARPETA_SALIDAS, f"informe_{base}-{k}{extension}")):
        k += 1
    return ruta


def _pedir_secuencia(opcion: str) -> SecuenciaLeida:
    if opcion == "1":
        print("\nEscribe la secuencia de ADN (solo bases A, C, G, T), por ejemplo:")
        print("  ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        return SecuenciaLeida("ADN_MANUAL", "Secuencia introducida manualmente por el usuario",
                              input("Secuencia: ").strip().upper())
    if opcion == "2":
        try:
            longitud = int(input("Longitud deseada en pb (por defecto 60, mínimo 6): ") or "60")
        except ValueError:
            longitud = 60
        datos = generar_adn_aleatorio(longitud)
        if len(datos.secuencia) != longitud:
            print(f"(La longitud se ha ajustado a {len(datos.secuencia)} pb, múltiplo de 3.)")
        return datos

    print(f"\nEjemplos disponibles en '{CARPETA_EJEMPLOS}':")
    if os.path.isdir(CARPETA_EJEMPLOS):
        for f in sorted(os.listdir(CARPETA_EJEMPLOS)):
            print(f"   - {f}")
    ruta = input("Ruta del fichero (.fasta/.fa o .gb/.gbk): ").strip()
    if not os.path.isabs(ruta) and not os.path.exists(ruta):
        ruta_ejemplo = os.path.join(CARPETA_EJEMPLOS, ruta)
        ruta = ruta_ejemplo if os.path.exists(ruta_ejemplo) else ruta
    return leer_secuencia(ruta)


def _guardar(pregunta: str, tipo: str, identificador: str, extension: str, escribir) -> bool:
    """Pregunta si guardar un informe y lo escribe; devuelve True si se guardó."""
    if input(pregunta).strip().lower() == "n":
        return False
    ruta = _ruta_libre(identificador, extension)
    try:
        escribir(ruta)
    except OSError as err:
        print(f"No se pudo guardar el informe {tipo}: {err}")
        return False
    print(f"Informe {'interactivo ' if tipo == 'HTML' else ''}guardado en: {ruta}")
    return True


def main() -> None:
    os.makedirs(CARPETA_SALIDAS, exist_ok=True)

    while True:
        print("\n" + "-" * 60)
        print(" SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR")
        print("-" * 60)
        print(" 1) Escribir una secuencia de ADN manualmente")
        print(" 2) Generar un gen de ADN aleatorio (ATG ... parada)")
        print(" 3) Cargar una secuencia real desde un fichero FASTA")
        print(" 4) Cargar una secuencia real desde un fichero GenBank")
        print(" 5) Salir")
        opcion = input("Elige una opción [1-5]: ").strip()

        if opcion == "5":
            print("¡Hasta luego!")
            sys.exit(0)
        if opcion not in ("1", "2", "3", "4"):
            print("Opción no válida, inténtalo de nuevo.")
            continue

        try:
            datos = _pedir_secuencia(opcion)
        except (ValueError, OSError) as err:
            print(f"\nError al obtener la secuencia: {err}")
            continue
        try:
            r = simular(datos.secuencia, datos.identificador, datos.descripcion)
        except ValueError as err:
            print(f"\nSecuencia no válida: {err}")
            continue

        n = len(r.codificante)
        lineas = construir_informe_texto(r)
        if n <= 300:
            print("\n" + "\n".join(lineas))
        else:
            print(f"\nSecuencia larga ({n} pb): el detalle completo se guarda en el informe.")
        print(f"\nSimulación completada - ADN: {n} pb, ARNm: {len(r.arnm)} nt, "
              f"Proteína: {len(r.proteina)} aa")

        _guardar("\n¿Guardar informe en un fichero de texto? [S/n]: ", "de texto",
                 r.nombre, ".txt", lambda ruta: guardar_informe_texto(ruta, lineas))
        if _guardar("¿Generar un informe HTML interactivo? [S/n]: ", "HTML",
                    r.nombre, ".html", lambda ruta: generar_informe_html(ruta, r)):
            print("   Ábrelo con cualquier navegador para verlo con pestañas e interactividad.")

        input("\nPulsa Enter para volver al menú...")


if __name__ == "__main__":
    main()
