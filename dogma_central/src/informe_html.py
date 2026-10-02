"""
informe_html.py
----------------
Genera un informe HTML autocontenido (sin dependencias externas: todo el
CSS y JS va embebido en un único fichero) que presenta de forma
interactiva y visual las tres etapas del dogma central:

    ADN -> ADN   (replicación:  helicasa, primasa, ADN polimerasa,
                  hebra líder, fragmentos de Okazaki, ADN ligasa)
    ADN -> ARN   (transcripción: ARN polimerasa, complementariedad)
    ARN -> proteína (traducción: ribosoma, ARNt, código genético)

Este informe complementa (no sustituye) la salida por consola y el
informe de texto: aquí se puede ver CADA etapa por separado (pestañas), con
las bases coloreadas, los fragmentos de Okazaki diferenciados y una
traducción codón a codón que se puede reproducir paso a paso.

No requiere matplotlib ni ninguna librería de terceros: es un único
fichero .html que se abre con cualquier navegador.
"""

from __future__ import annotations

import html
import json

from .adn import ResultadoReplicacion
from .arn import ResultadoTranscripcion
from .proteina import ResultadoTraduccion, NOMBRE_AMINOACIDO

# ---------------------------------------------------------------------- #
#  Clasificación de aminoácidos por propiedad química (para el color)
# ---------------------------------------------------------------------- #
GRUPO_AMINOACIDO = {
    **{aa: "apolar" for aa in "AVLIPFMWG"},
    **{aa: "polar" for aa in "STCYNQ"},
    **{aa: "acido" for aa in "DE"},
    **{aa: "basico" for aa in "KRH"},
}

ENZIMAS = {
    "replicacion": [
        ("Helicasa", "Abre y desenrolla la doble hélice de ADN, separando las dos hebras."),
        ("Topoisomerasa", "Alivia la tensión (superenrollamiento) que se acumula por delante de la horquilla."),
        ("Proteínas SSB", "Se unen a las hebras sencillas y evitan que se vuelvan a emparejar."),
        ("Primasa", "Sintetiza los cebadores (primers) de ARN necesarios para iniciar cada nueva hebra."),
        ("ADN polimerasa III", "Añade nucleótidos complementarios 5'→3' a partir de cada cebador (hebra líder y fragmentos de Okazaki)."),
        ("ADN polimerasa I", "Elimina los cebadores de ARN y rellena el hueco con ADN."),
        ("ADN ligasa", "Une entre sí los fragmentos de Okazaki de la hebra rezagada."),
    ],
    "transcripcion": [
        ("ARN polimerasa", "Lee la hebra molde de ADN (3'→5') y sintetiza el ARNm (5'→3')."),
        ("Factor sigma / factores de transcripción", "Reconocen el promotor y ayudan a iniciar la transcripción."),
    ],
    "traduccion": [
        ("Ribosoma", "Lee el ARNm codón a codón y cataliza la unión de aminoácidos."),
        ("ARNt (ARN de transferencia)", "Transporta cada aminoácido y lo empareja con su codón mediante el anticodón."),
        ("Aminoacil-ARNt sintetasa", "Une cada aminoácido a su ARNt correspondiente."),
        ("Factor de liberación", "Reconoce el codón de parada (UAA, UAG, UGA) y libera la proteína; no usa ARNt."),
    ],
}


def _esc(s: str) -> str:
    return html.escape(s, quote=True)


def _spans_bases(secuencia: str) -> str:
    """Envuelve cada base en un <span> coloreado por tipo de base."""
    return "".join(f'<span class="b {c}">{c}</span>' for c in secuencia)


def _doble_hebra_html(sup: str, inf: str, cada: int = 60) -> str:
    """
    Dos hebras emparejadas, troceadas en bloques de 'cada' bases para que
    sigan alineadas aunque la secuencia sea larga (cada bloque tiene una
    línea 5'->3' arriba y otra 3'<-5' debajo).
    """
    bloques = []
    for i in range(0, len(sup), cada):
        bloques.append(
            '<div class="doble">'
            f"<div>5'-{_spans_bases(sup[i:i + cada])}</div>"
            f"<div>3'-{_spans_bases(inf[i:i + cada])}</div>"
            "</div>"
        )
    return "".join(bloques)


def _bloque_okazaki(fragmentos) -> str:
    piezas = []
    for f in fragmentos:
        piezas.append(
            f'<div class="okazaki-frag">'
            f'<div class="okazaki-idx">Fragmento {f.indice}</div>'
            f'<span class="cebador" title="Cebador de ARN (primasa)">{_esc(f.cebador_arn)}</span>'
            f'<span class="unido">{_spans_bases(f.fragmento_adn)}</span>'
            f'<div class="okazaki-len">{f.longitud} nt</div>'
            f"</div>"
        )
    return "".join(piezas)


def _tabla_codones(codones) -> str:
    filas = []
    for i, c in enumerate(codones):
        grupo = "stop" if c.aminoacido == "*" else GRUPO_AMINOACIDO.get(c.aminoacido, "otro")
        nombre = "Codón STOP (fin de traducción)" if c.aminoacido == "*" else f"{NOMBRE_AMINOACIDO.get(c.aminoacido, c.aminoacido)} · anticodón ARNt 3'-{c.anticodon}-5'"
        filas.append(
            f'<div class="codon {grupo}" data-step="{i}" title="{_esc(nombre)}">'
            f'<div class="codon-triplete">{_esc(c.codon)}</div>'
            f'<div class="codon-aa">{_esc(c.aminoacido)}</div>'
            f"</div>"
        )
    return "".join(filas)


def _lista_enzimas(clave: str) -> str:
    return "".join(
        f'<span class="chip" title="{_esc(desc)}">{_esc(nombre)}</span>'
        for nombre, desc in ENZIMAS[clave]
    )


def generar_informe_html(
    ruta: str,
    nombre: str,
    descripcion: str,
    replicacion: ResultadoReplicacion,
    transcripcion: ResultadoTranscripcion,
    traduccion: ResultadoTraduccion,
) -> None:
    adn_len = len(replicacion.hebra_parental_codificante)
    arn_len = len(transcripcion.arnm)
    prot_len = len(traduccion.proteina)

    # Datos para el "reproductor" de traducción en JavaScript
    pasos_js = json.dumps(
        [
            {
                "codon": c.codon,
                "aa": c.aminoacido,
                "nombre": ("STOP" if c.aminoacido == "*" else NOMBRE_AMINOACIDO.get(c.aminoacido, c.aminoacido)),
                "pos": c.posicion,
            }
            for c in traduccion.codones
        ]
    )

    html_doc = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Dogma central — {_esc(nombre)}</title>
<style>
  :root {{
    --bg: #0f172a; --panel: #16213a; --panel2: #1c2b4a; --border: #2a3a5c;
    --text: #e6ebf5; --muted: #93a3c2; --accent: #38bdf8; --accent2: #a78bfa;
    --A: #f87171; --T: #60a5fa; --G: #4ade80; --C: #fbbf24; --U: #c084fc;
    --apolar: #60a5fa; --polar: #4ade80; --acido: #f87171; --basico: #fbbf24; --stop: #f43f5e;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; background: linear-gradient(180deg, var(--bg), #0b1223 60%);
    color: var(--text); font-family: 'Segoe UI', system-ui, sans-serif; padding: 0 0 60px;
  }}
  header {{
    padding: 28px 24px 20px; text-align: center;
    background: radial-gradient(circle at 50% -20%, rgba(56,189,248,.25), transparent 60%);
    border-bottom: 1px solid var(--border);
  }}
  header h1 {{ margin: 0 0 6px; font-size: 1.6rem; letter-spacing: .5px; }}
  header p {{ margin: 2px 0; color: var(--muted); font-size: .92rem; }}

  .flujo {{ display:flex; justify-content:center; gap:10px; margin-top:18px; flex-wrap:wrap; }}
  .flujo .paso {{ background: var(--panel2); border:1px solid var(--border); border-radius:10px;
    padding:10px 16px; font-size:.85rem; }}
  .flujo .flecha {{ color: var(--accent); font-size:1.3rem; align-self:center; }}

  nav.tabs {{ display:flex; justify-content:center; gap:8px; margin:24px auto 0; flex-wrap:wrap; max-width:1000px; }}
  nav.tabs button {{
    background: var(--panel); color: var(--text); border:1px solid var(--border);
    padding:10px 18px; border-radius:999px; cursor:pointer; font-size:.9rem; transition:.15s;
  }}
  nav.tabs button.activo {{ background: var(--accent); color:#08131f; font-weight:600; border-color: var(--accent); }}
  nav.tabs button:hover {{ border-color: var(--accent); }}

  main {{ max-width: 1000px; margin: 22px auto 0; padding: 0 20px; }}
  section.etapa {{ display:none; background: var(--panel); border:1px solid var(--border);
    border-radius:16px; padding:22px 24px; animation: fade .25s ease; }}
  section.etapa.activo {{ display:block; }}
  @keyframes fade {{ from {{opacity:0; transform: translateY(6px);}} to {{opacity:1; transform:none;}} }}

  h2 {{ margin-top:0; font-size:1.15rem; color: var(--accent); }}
  h3 {{ font-size: .95rem; color: var(--accent2); margin: 22px 0 8px; }}
  p.hint {{ color: var(--muted); font-size: .85rem; margin: 4px 0 14px; }}

  .chips {{ display:flex; flex-wrap:wrap; gap:8px; margin: 6px 0 16px; }}
  .chip {{ background: var(--panel2); border:1px solid var(--border); border-radius:999px;
    padding:5px 12px; font-size:.78rem; cursor:help; }}

  .seq {{ font-family: 'Cascadia Code', Consolas, monospace; line-height:1.9;
    word-break: break-all; background:#0b1424; border:1px solid var(--border);
    border-radius:10px; padding:12px 14px; font-size:.92rem; }}
  .b {{ padding: 1px 1px; border-radius:3px; }}
  .b.A {{ color: var(--A); }} .b.T {{ color: var(--T); }} .b.G {{ color: var(--G); }}
  .b.C {{ color: var(--C); }} .b.U {{ color: var(--U); }}
  .leyenda-bases {{ display:flex; gap:14px; margin: 8px 0 4px; font-size:.78rem; color:var(--muted); flex-wrap:wrap; }}
  .leyenda-bases span.dot {{ display:inline-block; width:9px; height:9px; border-radius:50%; margin-right:5px; }}

  .okazaki-track {{ display:flex; flex-wrap:wrap; gap:10px; margin: 6px 0 4px; }}
  .okazaki-frag {{ background:#0b1424; border:1px solid var(--border); border-radius:10px;
    padding:8px 10px; font-family: monospace; font-size:.82rem; min-width:150px; }}
  .okazaki-idx {{ color: var(--muted); font-size:.72rem; margin-bottom:4px; }}
  .okazaki-len {{ color: var(--muted); font-size:.7rem; margin-top:4px; }}
  .cebador {{ color: var(--U); border-bottom:2px dotted var(--U); margin-right:2px; }}

  .doble {{ font-family: 'Cascadia Code', Consolas, monospace; font-size:.92rem; line-height:1.5;
    white-space:nowrap; overflow-x:auto; background:#0b1424; border:1px solid var(--border);
    border-radius:10px; padding:8px 14px; margin:0 0 6px; }}
  h4 {{ font-size:.85rem; color: var(--text); margin: 14px 0 6px; }}

  .codones {{ display:flex; flex-wrap:wrap; gap:6px; margin: 6px 0 10px; }}
  .codon {{ border-radius:8px; padding:6px 8px; text-align:center; font-family:monospace;
    min-width:52px; border:1px solid var(--border); cursor:default; transition:.15s; }}
  .codon.apolar {{ background: rgba(96,165,250,.15); border-color: var(--apolar); }}
  .codon.polar {{ background: rgba(74,222,128,.15); border-color: var(--polar); }}
  .codon.acido {{ background: rgba(248,113,113,.15); border-color: var(--acido); }}
  .codon.basico {{ background: rgba(251,191,36,.15); border-color: var(--basico); }}
  .codon.stop {{ background: rgba(244,63,94,.2); border-color: var(--stop); }}
  .codon.actual {{ outline: 2px solid var(--accent); transform: scale(1.08); }}
  .codon-triplete {{ font-size:.78rem; color:var(--muted); }}
  .codon-aa {{ font-weight:700; font-size:1rem; }}

  .reproductor {{ background:#0b1424; border:1px solid var(--border); border-radius:12px;
    padding:14px 16px; margin-top:14px; }}
  .reproductor .controles {{ display:flex; gap:10px; align-items:center; margin-bottom:10px; flex-wrap:wrap; }}
  .reproductor button {{ background: var(--accent); border:none; color:#08131f; font-weight:600;
    padding:8px 14px; border-radius:8px; cursor:pointer; }}
  .reproductor button:disabled {{ opacity:.4; cursor:default; }}
  #estadoReproductor {{ font-family: monospace; font-size:.95rem; }}
  #proteinaCreciente {{ font-family: monospace; letter-spacing:2px; font-size:1.05rem; margin-top:8px; color: var(--accent2); }}

  .resumen-flujo {{ display:flex; justify-content:space-between; align-items:center;
    background:#0b1424; border:1px solid var(--border); border-radius:12px; padding:16px; margin-top:10px; flex-wrap:wrap; gap:10px; }}
  .resumen-item {{ text-align:center; }}
  .resumen-item .valor {{ font-size:1.4rem; font-weight:700; color:var(--accent); }}
  .resumen-item .etiqueta {{ font-size:.75rem; color:var(--muted); }}

  footer {{ text-align:center; color:var(--muted); font-size:.78rem; margin-top:30px; }}
</style>
</head>
<body>

<header>
  <h1>🧬 Simulación del dogma central de la biología molecular</h1>
  <p><strong>Secuencia:</strong> {_esc(nombre)} &nbsp;·&nbsp; <strong>Longitud:</strong> {adn_len} pb</p>
  {f'<p>{_esc(descripcion)}</p>' if descripcion else ""}
  <div class="flujo">
    <div class="paso">ADN ({adn_len} pb)</div>
    <div class="flecha">↺</div>
    <div class="paso">ADN replicado</div>
    <div class="flecha">→</div>
    <div class="paso">ARNm ({arn_len} nt)</div>
    <div class="flecha">→</div>
    <div class="paso">Proteína ({prot_len} aa)</div>
  </div>
</header>

<nav class="tabs">
  <button data-tab="replicacion" class="activo">1 · Replicación</button>
  <button data-tab="transcripcion">2 · Transcripción</button>
  <button data-tab="traduccion">3 · Traducción</button>
  <button data-tab="resumen">Resumen</button>
</nav>

<main>

<section class="etapa activo" id="tab-replicacion">
  <h2>Replicación del ADN &nbsp;(ADN → ADN)</h2>
  <p class="hint">La helicasa abre la doble hélice; sobre cada hebra molde actúan la primasa y la ADN polimerasa.</p>
  <div class="chips">{_lista_enzimas("replicacion")}</div>

  <h3>1 · Apertura de la doble hélice</h3>
  <p class="hint">La helicasa rompe los puentes de hidrógeno; la topoisomerasa alivia la tensión por delante de la horquilla y las proteínas SSB mantienen las hebras separadas. Arriba: hebra codificante (5'→3'). Abajo: hebra molde (3'←5').</p>
  {_doble_hebra_html(replicacion.hebra_parental_codificante, replicacion.hebra_parental_molde)}
  <div class="leyenda-bases">
    <span><span class="dot" style="background:var(--A)"></span>A</span>
    <span><span class="dot" style="background:var(--T)"></span>T</span>
    <span><span class="dot" style="background:var(--G)"></span>G</span>
    <span><span class="dot" style="background:var(--C)"></span>C</span>
  </div>

  <h3>2 · Hebra líder — síntesis continua</h3>
  <p class="hint">Molde: hebra molde (3'→5'), leída en el mismo sentido que avanza la horquilla. Cebador de ARN: <span class="cebador">5'-{_esc(replicacion.cebador_hebra_lider)}-3'</span>, extendido sin interrupción por la ADN polimerasa III.</p>
  <div class="seq">5'-{_spans_bases(replicacion.hebra_lider_nueva)}-3'</div>

  <h3>3 · Hebra rezagada — {len(replicacion.fragmentos_okazaki)} fragmentos de Okazaki</h3>
  <p class="hint">Molde: hebra codificante. Cada fragmento se sintetiza 5'→3', en sentido contrario a la horquilla, y empieza con su propio cebador de ARN (subrayado); la ADN polimerasa III lo extiende. Después, la ADN polimerasa I sustituye los cebadores por ADN y la ADN ligasa une los fragmentos.</p>
  <div class="okazaki-track">{_bloque_okazaki(replicacion.fragmentos_okazaki)}</div>

  <h3>4 · Resultado: dos moléculas hijas semiconservativas</h3>
  <p class="hint">Cada molécula hija conserva una hebra parental y lleva una hebra nueva (arriba 5'→3', abajo 3'←5').</p>
  <h4>Hija 1 = hebra líder nueva (arriba) + hebra molde parental (abajo)</h4>
  {_doble_hebra_html(*replicacion.molecula_hija_1[::-1])}
  <h4>Hija 2 = hebra codificante parental (arriba) + hebra rezagada nueva (abajo)</h4>
  {_doble_hebra_html(*replicacion.molecula_hija_2)}
</section>

<section class="etapa" id="tab-transcripcion">
  <h2>Transcripción &nbsp;(ADN → ARN)</h2>
  <p class="hint">La ARN polimerasa lee la hebra molde (3'→5') y sintetiza el ARNm (5'→3'), sustituyendo T por U. En bacterias, el factor sigma guía a la polimerasa hasta el promotor.</p>
  <div class="chips">{_lista_enzimas("transcripcion")}</div>

  <h3>Hebra molde de ADN utilizada</h3>
  <div class="seq">3'-{_spans_bases(transcripcion.hebra_molde_usada)}-5'</div>

  <h3>ARN mensajero sintetizado</h3>
  <div class="seq">5'-{_spans_bases(transcripcion.arnm)}-3'</div>
  <div class="leyenda-bases">
    <span><span class="dot" style="background:var(--A)"></span>A</span>
    <span><span class="dot" style="background:var(--U)"></span>U</span>
    <span><span class="dot" style="background:var(--G)"></span>G</span>
    <span><span class="dot" style="background:var(--C)"></span>C</span>
  </div>
  <p class="hint">Regla de complementariedad aplicada: A↔U, T↔A, C↔G, G↔C.</p>
</section>

<section class="etapa" id="tab-traduccion">
  <h2>Traducción &nbsp;(ARN → proteína)</h2>
  <p class="hint">El ribosoma lee el ARNm codón a codón desde el AUG de inicio hasta un codón de parada; el ARNt aporta cada aminoácido.</p>
  <div class="chips">{_lista_enzimas("traduccion")}</div>

  <h3>Lectura codón a codón</h3>
  <div class="codones" id="pistaCodones">{_tabla_codones(traduccion.codones)}</div>

  <div class="reproductor">
    <div class="controles">
      <button id="btnPlay">▶ Reproducir traducción</button>
      <button id="btnReset" disabled>⟲ Reiniciar</button>
      <span id="estadoReproductor">Pulsa reproducir para ver el ribosoma leyendo codón a codón</span>
    </div>
    <div id="proteinaCreciente"></div>
  </div>

  <h3 style="margin-top:22px">Proteína resultante</h3>
  <div class="seq" style="letter-spacing:3px">{_esc(traduccion.proteina) or "(no se tradujo ningún aminoácido)"}</div>
  <p class="hint">{prot_len} aminoácidos{" · se alcanzó un codón de parada" if traduccion.parada_encontrada else " · no se encontró codón de parada (el ARNm se agotó)"}</p>
</section>

<section class="etapa" id="tab-resumen">
  <h2>Resumen del flujo de información genética</h2>
  <div class="resumen-flujo">
    <div class="resumen-item"><div class="valor">{adn_len} pb</div><div class="etiqueta">ADN</div></div>
    <div class="flecha">→</div>
    <div class="resumen-item"><div class="valor">{arn_len} nt</div><div class="etiqueta">ARNm</div></div>
    <div class="flecha">→</div>
    <div class="resumen-item"><div class="valor">{prot_len} aa</div><div class="etiqueta">Proteína</div></div>
  </div>
  <h3>Enzimas y moléculas que intervienen en todo el proceso</h3>
  <div class="chips">{_lista_enzimas("replicacion")}{_lista_enzimas("transcripcion")}{_lista_enzimas("traduccion")}</div>
</section>

</main>

<footer>Informe generado automáticamente por el simulador del dogma central · Bioinformática · ULPGC</footer>

<script>
  const tabs = document.querySelectorAll("nav.tabs button");
  tabs.forEach(btn => btn.addEventListener("click", () => {{
    tabs.forEach(b => b.classList.remove("activo"));
    document.querySelectorAll("section.etapa").forEach(s => s.classList.remove("activo"));
    btn.classList.add("activo");
    document.getElementById("tab-" + btn.dataset.tab).classList.add("activo");
  }}));

  const pasos = {pasos_js};
  let i = -1, timer = null;
  const btnPlay = document.getElementById("btnPlay");
  const btnReset = document.getElementById("btnReset");
  const estado = document.getElementById("estadoReproductor");
  const proteinaEl = document.getElementById("proteinaCreciente");
  const codonEls = document.querySelectorAll("#pistaCodones .codon");

  function pintarPaso() {{
    codonEls.forEach(el => el.classList.remove("actual"));
    if (i < 0 || i >= pasos.length) return;
    codonEls[i].classList.add("actual");
    codonEls[i].scrollIntoView({{behavior:"smooth", block:"nearest", inline:"center"}});
    const p = pasos[i];
    estado.textContent = `Codón ${{p.pos}}: ${{p.codon}} → ${{p.nombre}}`;
    proteinaEl.textContent = pasos.slice(0, i+1).map(x => x.aa === "*" ? "" : x.aa).join(" - ");
  }}

  function siguientePaso() {{
    i++;
    if (i >= pasos.length) {{
      clearInterval(timer); timer = null;
      btnPlay.textContent = "▶ Reproducir traducción";
      estado.textContent = "Traducción completada.";
      return;
    }}
    pintarPaso();
  }}

  btnPlay.addEventListener("click", () => {{
    if (timer) {{
      clearInterval(timer); timer = null;
      btnPlay.textContent = "▶ Reproducir traducción";
      return;
    }}
    if (i >= pasos.length - 1) i = -1;
    btnPlay.textContent = "⏸ Pausar";
    btnReset.disabled = false;
    timer = setInterval(siguientePaso, 650);
    siguientePaso();
  }});

  btnReset.addEventListener("click", () => {{
    clearInterval(timer); timer = null;
    i = -1;
    codonEls.forEach(el => el.classList.remove("actual"));
    proteinaEl.textContent = "";
    estado.textContent = "Pulsa reproducir para ver el ribosoma leyendo codón a codón";
    btnPlay.textContent = "▶ Reproducir traducción";
  }});
</script>
</body>
</html>
"""

    with open(ruta, "w", encoding="utf-8") as f:
        f.write(html_doc)