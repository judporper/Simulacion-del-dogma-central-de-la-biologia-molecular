"""
visualizador.py
----------------
Construye el contenido del informe de texto plano (.txt) con el
resultado completo de las tres etapas del dogma central:

    ADN -> ADN  (replicación)
    ADN -> ARN  (transcripción)
    ARN -> proteína (traducción)

El resultado ya NO se imprime por consola (no aportaba valor verlo
volcado en la terminal): se guarda directamente en un fichero de texto
y, si se quiere, en el informe HTML interactivo de `informe_html.py`.
"""

from __future__ import annotations

import textwrap

from .adn import ResultadoReplicacion
from .arn import ResultadoTranscripcion
from .proteina import ResultadoTraduccion, NOMBRE_AMINOACIDO

ANCHO_LINEA = 70


def _envolver(secuencia: str, cada: int = 60) -> str:
    return "\n".join(textwrap.wrap(secuencia, cada)) if secuencia else "(vacía)"

def _doble_hebra(sup: str, inf: str, cada: int = 60) -> list[str]:
    """Dos hebras emparejadas, troceadas en bloques para que sigan alineadas."""
    out: list[str] = []
    for i in range(0, len(sup), cada):
        out.append(f"  5'-{sup[i:i + cada]}")
        out.append(f"  3'-{inf[i:i + cada]}")
        out.append("")
    return out

def construir_informe_texto(
    nombre: str,
    descripcion: str,
    longitud_adn: int,
    replicacion: ResultadoReplicacion,
    transcripcion: ResultadoTranscripcion,
    traduccion: ResultadoTraduccion,
) -> list[str]:
    """Construye, línea a línea, el contenido íntegro del informe de texto."""
    barra = "=" * ANCHO_LINEA
    l: list[str] = []

    l += [barra, " SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR", barra, ""]
    l.append(f"Secuencia: {nombre}")
    if descripcion:
        l.append(f"Descripción: {descripcion}")
    l.append(f"Longitud: {longitud_adn} pb")
    l.append("")

    # ETAPA 1 · Replicación
    l += [barra, " ETAPA 1 - REPLICACION DEL ADN  (ADN -> ADN)", barra]
    l.append("Enzimas/moléculas implicadas: helicasa, primasa, ADN polimerasa, ADN ligasa")
    l.append("Modelo semiconservativo; la horquilla de replicación avanza de izquierda a derecha.")
    l.append("")
    l.append("1. Apertura de la doble hélice (helicasa rompe los puentes de hidrógeno)")
    l.append("  Arriba: hebra codificante (5'->3')  |  Abajo: hebra molde (3'<-5')")
    l += _doble_hebra(replicacion.hebra_parental_codificante, replicacion.hebra_parental_molde)
    l.append("2. Hebra líder: síntesis continua")
    l.append("  Molde: hebra molde (se lee 3'->5', en el mismo sentido que avanza la horquilla)")
    l.append(f"  Cebador de ARN (primasa): 5'-{replicacion.cebador_hebra_lider}-3'")
    l.append("  La ADN polimerasa extiende el cebador sin interrupción, en sentido 5'->3'")
    l.append(f"  Nueva hebra líder (5'->3'): {_envolver(replicacion.hebra_lider_nueva)}")
    l.append("")
    l.append(f"3. Hebra rezagada: {len(replicacion.fragmentos_okazaki)} fragmentos de Okazaki")
    l.append("  Molde: hebra codificante. Cada fragmento se sintetiza en sentido 5'->3', contrario")
    l.append("  al avance de la horquilla, y empieza con su propio cebador de ARN (primasa).")
    for frag in replicacion.fragmentos_okazaki:
        l.append(
            f"   - Fragmento {frag.indice:>3} | molde pos. {frag.inicio_molde + 1}-{frag.fin_molde} "
            f"| cebador ARN 5'-{frag.cebador_arn}-3' + ADN {frag.fragmento_adn} ({frag.longitud} nt)"
        )
    l.append("  La ADN polimerasa I sustituye los cebadores por ADN y la ADN ligasa sella los huecos")
    l.append("  -> hebra rezagada completa (alineada con su molde, 3'<-5'):")
    l.append(f"  {_envolver(replicacion.hebra_rezagada_nueva)}")
    l.append("")
    l.append("4. Resultado: dos moléculas hijas semiconservativas (cada una con una hebra vieja y una nueva)")
    l.append("  Hija 1: hebra líder nueva (arriba, 5'->3') + hebra molde parental (abajo, 3'<-5')")
    l += _doble_hebra(replicacion.hebra_lider_nueva, replicacion.hebra_parental_molde)
    l.append("  Hija 2: hebra codificante parental (arriba, 5'->3') + hebra rezagada nueva (abajo, 3'<-5')")
    l += _doble_hebra(replicacion.hebra_parental_codificante, replicacion.hebra_rezagada_nueva)

    # ETAPA 2 · Transcripción
    l += [barra, " ETAPA 2 - TRANSCRIPCION  (ADN -> ARN)", barra]
    l.append("Enzimas/moléculas implicadas: ARN polimerasa, ribonucleótidos (A, U, G, C)")
    l.append("")
    l.append("Iniciación: la ARN polimerasa se une al promotor y abre una burbuja de transcripción.")
    l.append("  (Simplificación: se transcribe toda la secuencia, sin modelar el promotor.)")
    l.append("Elongación: lee la hebra molde en sentido 3'->5' y añade ribonucleótidos")
    l.append("  complementarios, sintetizando el ARNm en sentido 5'->3'.")
    l.append("")
    l.append("Molde utilizado (hebra de ADN, 3'->5'):")
    l.append(f"  {_envolver(transcripcion.hebra_molde_usada)}")
    l.append("")
    l.append("ARN mensajero sintetizado (5'->3'):")
    l.append(f"  {_envolver(transcripcion.arnm)}")
    l.append("")
    l.append("(Regla de complementariedad aplicada: A<->U, T<->A, C<->G, G<->C)")
    coincide = transcripcion.arnm == transcripcion.hebra_codificante.replace("T", "U")
    l.append(
        "Comprobación: el ARNm coincide con la hebra codificante cambiando T por U: "
        + ("sí" if coincide else "NO")
    )
    l.append("Terminación: la polimerasa se libera y el ARNm queda listo para la traducción.")
    l.append("")

    # ETAPA 3 · Traducción
    l += [barra, " ETAPA 3 - TRADUCCION  (ARN -> PROTEINA)", barra]
    l.append("Enzimas/moléculas implicadas: ribosoma, ARNt, aminoacil-ARNt sintetasa")
    l.append("")
    if traduccion.posicion_inicio == -1:
        l.append("No se ha encontrado ningún codón de inicio (AUG) en el ARNm.")
    else:
        l.append(f"Iniciación: el ribosoma localiza el codón AUG en la posición {traduccion.posicion_inicio}")
        l.append("Elongación: cada codón es reconocido por un ARNt con el anticodón complementario")
        l.append("  (la aminoacil-ARNt sintetasa cargó antes cada ARNt con su aminoácido):")
        l.append("")
        for i, c in enumerate(traduccion.codones):
            if c.aminoacido == "*":
                detalle = "sin ARNt: factor de liberación -> TERMINACIÓN"
            else:
                detalle = (
                    f"ARNt anticodón 3'-{c.anticodon}-5' | {c.aminoacido:<2} "
                    f"{NOMBRE_AMINOACIDO.get(c.aminoacido, c.aminoacido)}"
                )
                if i == 0:
                    detalle += " (inicio)"
            l.append(f"   pos {c.posicion:>4} | codón {c.codon} | {detalle}")
        l.append("")
        l.append("Secuencia de aminoácidos (proteína resultante):")
        l.append(f"  {_envolver(traduccion.proteina)}")
        l.append(f"  Longitud: {len(traduccion.proteina)} aminoácidos")
        if not traduccion.parada_encontrada:
            l.append("  (No se alcanzó codón de parada: el ARNm se agotó antes)")
    l.append("")
    
    # Resumen final
    l += [barra, " RESUMEN DEL FLUJO DE INFORMACION GENETICA", barra]
    l.append(
        f"  ADN ({longitud_adn} pb)  ->  ARNm ({len(transcripcion.arnm)} nt)  ->  "
        f"Proteína ({len(traduccion.proteina)} aa)"
    )

    return l


def guardar_informe_texto(ruta: str, lineas: list[str]) -> None:
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas))
