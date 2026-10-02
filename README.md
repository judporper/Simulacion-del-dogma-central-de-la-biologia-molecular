# Simulador del Dogma Central de la Biología Molecular

**Bioinformática — Práctica 1**
ULPGC · Grado en Ciencia e Ingeniería de Datos

Simulador en **Python puro y sin dependencias externas** que representa de
forma integrada los tres procesos del dogma central y muestra cómo viaja la
información desde el ADN hasta la proteína:

```
   ADN  ──replicación──▶  ADN
    │
    └──transcripción──▶  ARNm ──traducción──▶  Proteína
```

## Objetivos cubiertos

| Objetivo de la práctica | Dónde se implementa |
|---|---|
| Replicación: apertura de la doble hélice, cebadores, enzimas, hebra líder, hebra rezagada y fragmentos de Okazaki | `src/adn.py` → `ADN.replicar()` |
| Transcripción: hebra molde → ARNm, complementariedad de bases | `src/arn.py` → `transcribir()` |
| Traducción: lectura de codones, código genético, ARNt, señal de terminación | `src/proteina.py` → `traducir()` |
| Papel de las moléculas y enzimas de cada etapa | Informe de texto/consola (`visualizador.py`) y fichas del informe HTML (`informe_html.py`) |
| Visualización de cada etapa (texto o imagen) | Consola y `.txt` (`visualizador.py`); informe HTML con pestañas (`informe_html.py`) |
| Lenguaje e interacción de libre elección | Python + menú de consola (`main.py`) |
| Datos biológicos reales (FASTA / GenBank) | `src/lector.py` + `data/ejemplos/` (gen *lacZ* de *E. coli*) |

## Qué simula cada etapa

**1. Replicación (ADN → ADN).** Modelo semiconservativo, con la horquilla
avanzando de izquierda a derecha.
- La **helicasa** abre la doble hélice y separa la hebra codificante (5'→3') de
  la hebra molde (3'←5').
- **Hebra líder:** la **primasa** pone un único cebador de ARN y la **ADN
  polimerasa** sintetiza de forma continua sobre la hebra molde, en el mismo
  sentido que avanza la horquilla.
- **Hebra rezagada:** el molde es la hebra codificante. Cada **fragmento de
  Okazaki** se sintetiza en sentido 5'→3', contrario a la horquilla, con su
  propio cebador de ARN (complementario a su molde). Después la **ADN
  polimerasa I** sustituye los cebadores por ADN y la **ADN ligasa** une los
  fragmentos.
- Resultado: dos moléculas hijas, cada una con una hebra parental y otra nueva
  (hija 1 = molde + líder; hija 2 = codificante + rezagada).

**2. Transcripción (ADN → ARNm).** La **ARN polimerasa** usa la hebra molde,
la lee en sentido 3'→5' y sintetiza el ARNm en 5'→3' (A↔U, T↔A, C↔G, G↔C). El
programa comprueba que el ARNm coincide con la hebra codificante cambiando T
por U.

**3. Traducción (ARNm → proteína).** El **ribosoma** localiza el primer AUG,
lee tripletes y, para cada codón, se indica el **ARNt** con su anticodón
(cargado antes por la aminoacil-ARNt sintetasa). Termina en el primer codón de
parada (UAA, UAG, UGA), reconocido por un factor de liberación y no por un ARNt.

## Estructura del proyecto

```
dogma_central/
├── main.py                  # Punto de entrada: menú interactivo por consola
├── src/
│   ├── adn.py               # Clase ADN + simulación de la REPLICACIÓN
│   ├── arn.py               # Simulación de la TRANSCRIPCIÓN
│   ├── proteina.py          # Código genético + simulación de la TRADUCCIÓN
│   ├── lector.py            # Lectura de FASTA / GenBank / ADN aleatorio
│   ├── visualizador.py      # Construye el informe de texto plano
│   ├── informe_html.py      # Informe HTML interactivo (pestañas, colores, reproducción)
│   └── simulador.py         # Orquesta las tres etapas
├── data/ejemplos/
│   ├── lacZ_Ecoli_NC_000913.3.fasta   # Gen lacZ real (NCBI), formato FASTA
│   └── lacZ_Ecoli_NC_000913.3.gb      # Mismo gen, formato GenBank
└── salidas/                 # Aquí se guardan los informes generados
```

## Cómo ejecutarlo en Visual Studio Code

1. Abre la carpeta `dogma_central/` en VS Code (`Archivo → Abrir carpeta…`).
2. Necesitas **Python 3.10+** y seleccionar el intérprete (`Ctrl+Shift+P` →
   *Python: Select Interpreter*). No hay que instalar nada más.
3. Abre una terminal integrada (`` Ctrl+ñ ``) y ejecuta:
```bash
   python main.py
```
4. Elige una opción del menú:
   - **1** → escribes tú mismo una secuencia de ADN (solo A, C, G, T).
   - **2** → se genera una secuencia aleatoria (rápido para probar).
   - **3** → cargas un fichero FASTA (por ejemplo el `.fasta` de ejemplo).
   - **4** → cargas un fichero GenBank (`.gb`).
   - **5** → salir.

## Qué se obtiene

Al terminar la simulación:
- **Por consola** se muestran las tres etapas completas si la secuencia tiene
  hasta 300 pb. Con secuencias más largas (como *lacZ*, 3075 pb) solo se
  muestra un resumen, y el detalle va al informe.
- Se ofrece guardar un **informe `.txt`** con las dos hebras emparejadas, los
  cebadores y fragmentos de Okazaki, la lectura codón a codón con su anticodón
  y la proteína final.
- Se ofrece generar un **informe HTML interactivo** (se abre en cualquier
  navegador): una pestaña por etapa, bases coloreadas, fragmentos de Okazaki
  con su cebador resaltado, fichas de las enzimas con su función al pasar el
  ratón y un botón para **reproducir la traducción codón a codón** viendo
  crecer la proteína.

Ambos informes se guardan en `salidas/` sin sobrescribir los anteriores
(`informe_<id>-1`, `-2`, …).

### Ejemplo (secuencia corta)

Para `ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG`:

```
ADN (39 pb)  ->  ARNm (39 nt)  ->  Proteína (7 aa)
ARNm:      AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Proteína:  MAIVMGR   (termina en el codón UGA)
```

Con el gen *lacZ* de ejemplo (3075 pb) se obtiene la β-galactosidasa de 1024
aminoácidos (MTMITDSLAVVLQRRDWENPG…), que termina con codón de parada.

## Probar con tus propios ficheros

Copia cualquier `.fasta`/`.fa`/`.fna` o `.gb`/`.gbk`/`.genbank` a
`data/ejemplos/` (o da la ruta completa cuando el programa la pida). El formato
se detecta por la extensión. Se lee solo el primer registro y se descartan las
bases ambiguas (N, etc.).

## Convenios y simplificaciones

- **Orientación:** la hebra codificante se guarda en 5'→3' y la molde es su
  complementaria alineada base a base (3'←5'). La horquilla avanza de izquierda
  a derecha.
- **Fragmentos de Okazaki:** por defecto miden 10 nt (parámetro
  `tamano_fragmento_okazaki`, siempre mayor que el cebador de 6 nt). Es un
  tamaño didáctico: en bacterias miden 1000-2000 nt.
- **Cebadores:** son ARN complementario al molde en ese punto, de 6 nt.
- **Transcripción:** se transcribe toda la secuencia a partir de la hebra
  molde, sin modelar promotor, terminador ni procesamiento (cap, poliA, splicing).
  La teoría de iniciación, elongación y terminación se indica en el informe.
- **Traducción:** empieza en el primer AUG y acaba en el primer codón de
  parada, con el código genético estándar (tabla 1 de NCBI). La tabla 11
  (bacteriana) apenas difiere de ella.
- **Sin dependencias:** los lectores de FASTA y GenBank están escritos a mano
  (sin Biopython) y el informe HTML es autocontenido, así que no hace falta
  `pip install`.