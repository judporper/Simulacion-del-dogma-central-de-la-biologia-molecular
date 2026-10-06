from __future__ import annotations

import textwrap

from .dogma import NOMBRE_AMINOACIDO, Simulacion

ANCHO_LINEA = 70


def _envolver(secuencia: str, cada: int = 60) -> str:
    return "\n".join(textwrap.wrap(secuencia, cada)) if secuencia else "(vacía)"

def _doble_hebra(sup: str, inf: str, cada: int = 60) -> list[str]:
    out: list[str] = []
    for i in range(0, len(sup), cada):
        out.append(f"  5'-{sup[i:i + cada]}")
        out.append(f"  3'-{inf[i:i + cada]}")
        out.append("")
    return out

def construir_informe_texto(r: Simulacion) -> list[str]:
    barra = "=" * ANCHO_LINEA
    l: list[str] = []

    l += [barra, " SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR", barra, ""]
    l.append(f"Secuencia: {r.nombre}")
    if r.descripcion:
        l.append(f"Descripción: {r.descripcion}")
    longitud_adn = len(r.codificante)
    l.append(f"Longitud: {longitud_adn} pb")
    l.append("")

    l += [barra, " ETAPA 1 - REPLICACION DEL ADN  (ADN -> ADN)", barra]
    l.append("Enzimas/moléculas implicadas: helicasa, topoisomerasa, proteínas SSB, primasa,")
    l.append("  ADN polimerasa III, ADN polimerasa I, ADN ligasa")
    l.append("Modelo semiconservativo; la horquilla de replicación avanza de izquierda a derecha.")
    l.append("")
    l.append("1. Apertura de la doble hélice (helicasa rompe los puentes de hidrógeno)")
    l.append("  La topoisomerasa alivia la tensión por delante de la horquilla y las proteínas SSB")
    l.append("  mantienen las hebras separadas e impiden que vuelvan a emparejarse.")
    l.append("  Arriba: hebra codificante (5'->3')  |  Abajo: hebra molde (3'<-5')")
    l += _doble_hebra(r.codificante, r.molde)
    l.append("2. Hebra líder: síntesis continua")
    l.append("  Molde: hebra molde (se lee 3'->5', en el mismo sentido que avanza la horquilla)")
    l.append(f"  Cebador de ARN (primasa): 5'-{r.cebador_lider}-3'")
    l.append("  La ADN polimerasa III extiende el cebador sin interrupción, en sentido 5'->3'")
    l.append(f"  Nueva hebra líder (5'->3'): {_envolver(r.lider)}")
    l.append("")
    l.append(f"3. Hebra rezagada: {len(r.okazaki)} fragmentos de Okazaki")
    l.append("  Molde: hebra codificante. Cada fragmento se sintetiza en sentido 5'->3', contrario")
    l.append("  al avance de la horquilla, y empieza con su propio cebador de ARN (primasa);")
    l.append("  la ADN polimerasa III los extiende.")
    for frag in r.okazaki:
        l.append(
            f"   - Fragmento {frag.indice:>3} | molde pos. {frag.inicio_molde + 1}-{frag.fin_molde} "
            f"| cebador ARN 5'-{frag.cebador_arn}-3' + ADN {frag.fragmento_adn} ({frag.longitud} nt)"
        )
    l.append("  La ADN polimerasa I sustituye los cebadores por ADN y la ADN ligasa sella los huecos")
    l.append("  -> hebra rezagada completa (alineada con su molde, 3'<-5'):")
    l.append(f"  {_envolver(r.rezagada)}")
    l.append("")
    l.append("4. Resultado: dos moléculas hijas semiconservativas (cada una con una hebra vieja y una nueva)")
    l.append("  Hija 1: hebra líder nueva (arriba, 5'->3') + hebra molde parental (abajo, 3'<-5')")
    l += _doble_hebra(r.lider, r.molde)
    l.append("  Hija 2: hebra codificante parental (arriba, 5'->3') + hebra rezagada nueva (abajo, 3'<-5')")
    l += _doble_hebra(r.codificante, r.rezagada)

    l += [barra, " ETAPA 2 - TRANSCRIPCION  (ADN -> ARN)", barra]
    l.append("Enzimas/moléculas implicadas: ARN polimerasa, factor sigma (reconoce el promotor en")
    l.append("  bacterias), ribonucleótidos (A, U, G, C)")
    l.append("")
    l.append("Iniciación: la ARN polimerasa (guiada por el factor sigma) se une al promotor y abre")
    l.append("  una burbuja de transcripción.")
    l.append("  (Simplificación: se transcribe toda la secuencia, sin modelar el promotor.)")
    l.append("Elongación: lee la hebra molde en sentido 3'->5' y añade ribonucleótidos")
    l.append("  complementarios, sintetizando el ARNm en sentido 5'->3'.")
    l.append("")
    l.append("Molde utilizado (hebra de ADN, 3'->5'):")
    l.append(f"  {_envolver(r.molde)}")
    l.append("")
    l.append("ARN mensajero sintetizado (5'->3'):")
    l.append(f"  {_envolver(r.arnm)}")
    l.append("")
    l.append("(Regla de complementariedad aplicada: A<->U, T<->A, C<->G, G<->C)")
    coincide = r.arnm == r.codificante.replace("T", "U")
    l.append(
        "Comprobación: el ARNm coincide con la hebra codificante cambiando T por U: "
        + ("sí" if coincide else "NO")
    )
    l.append("Terminación: la polimerasa se libera y el ARNm queda listo para la traducción.")
    l.append("")

    l += [barra, " ETAPA 3 - TRADUCCION  (ARN -> PROTEINA)", barra]
    l.append("Enzimas/moléculas implicadas: ribosoma, ARNt, aminoacil-ARNt sintetasa, factor de liberación")
    l.append("")
    if r.inicio == -1:
        l.append("No se ha encontrado ningún codón de inicio (AUG) en el ARNm.")
    else:
        l.append(f"Iniciación: el ribosoma localiza el codón AUG en la posición {r.inicio}")
        l.append("Elongación: cada codón es reconocido por un ARNt con el anticodón complementario")
        l.append("  (la aminoacil-ARNt sintetasa cargó antes cada ARNt con su aminoácido):")
        l.append("")
        for i, c in enumerate(r.codones):
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
        l.append(f"  {_envolver(r.proteina)}")
        l.append(f"  Longitud: {len(r.proteina)} aminoácidos")
        if not r.parada:
            l.append("  (No se alcanzó codón de parada: el ARNm se agotó antes)")
    l.append("")
    
    l += [barra, " RESUMEN DEL FLUJO DE INFORMACION GENETICA", barra]
    l.append(
        f"  ADN ({longitud_adn} pb)  ->  ARNm ({len(r.arnm)} nt)  ->  "
        f"Proteína ({len(r.proteina)} aa)"
    )

    return l


def guardar_informe_texto(ruta: str, lineas: list[str]) -> None:
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas))
