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
    l.append("")
    l.append("1. Apertura de la doble hélice (helicasa)")
    l.append(f"  Hebra codificante (5'->3'): {_envolver(replicacion.hebra_parental_codificante)}")
    l.append(f"  Hebra molde        (3'<-5'): {_envolver(replicacion.hebra_parental_molde)}")
    l.append("")
    l.append("2. Hebra líder: síntesis continua")
    l.append(f"  Cebador de ARN (primasa): 5'-{replicacion.cebador_hebra_lider}-3'")
    l.append("  Molde usado: hebra codificante (3'->5')")
    l.append(f"  Nueva hebra líder (5'->3'): {_envolver(replicacion.hebra_lider_nueva)}")
    l.append("")
    l.append(f"3. Hebra rezagada: {len(replicacion.fragmentos_okazaki)} fragmentos de Okazaki")
    l.append("  Molde usado: hebra molde (5'->3')")
    for frag in replicacion.fragmentos_okazaki:
        l.append(
            f"   - Fragmento {frag.indice:>2} | cebador 5'-{frag.cebador_arn}-3' "
            f"| ADN: {frag.fragmento_adn} ({frag.longitud} nt)"
        )
    l.append("  ADN ligasa une los fragmentos -> hebra rezagada completa:")
    l.append(f"  {_envolver(replicacion.hebra_rezagada_nueva)}")
    l.append("")
    l.append("4. Resultado: dos moléculas hijas semiconservativas")
    l.append("  Hija 1 = hebra parental codificante + hebra líder nueva")
    l.append("  Hija 2 = hebra parental molde       + hebra rezagada nueva")
    l.append("")

    # ETAPA 2 · Transcripción
    l += [barra, " ETAPA 2 - TRANSCRIPCION  (ADN -> ARN)", barra]
    l.append("Enzimas/moléculas implicadas: ARN polimerasa")
    l.append("")
    l.append("Molde utilizado (hebra de ADN, 3'->5'):")
    l.append(f"  {_envolver(transcripcion.hebra_molde_usada)}")
    l.append("")
    l.append("ARN mensajero sintetizado (5'->3'):")
    l.append(f"  {_envolver(transcripcion.arnm)}")
    l.append("")
    l.append("(Regla de complementariedad aplicada: A<->U, T<->A, C<->G, G<->C)")
    l.append("")

    # ETAPA 3 · Traducción
    l += [barra, " ETAPA 3 - TRADUCCION  (ARN -> PROTEINA)", barra]
    l.append("Enzimas/moléculas implicadas: ribosoma, ARNt, aminoacil-ARNt sintetasa")
    l.append("")
    if traduccion.posicion_inicio == -1:
        l.append("No se ha encontrado ningún codón de inicio (AUG) en el ARNm.")
    else:
        l.append(f"Codón de inicio AUG localizado en la posición {traduccion.posicion_inicio}")
        l.append("")
        l.append("Lectura del ribosoma, codón a codón:")
        for c in traduccion.codones:
            aa_nombre = NOMBRE_AMINOACIDO.get(c.aminoacido, "STOP (fin de traducción)")
            marca = "PARADA" if c.aminoacido == "*" else c.aminoacido
            l.append(f"   pos {c.posicion:>4} | codón {c.codon} -> {marca:<8} {aa_nombre}")
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
