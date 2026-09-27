# Simulador del Dogma Central de la Biología Molecular

**Bioinformática — Práctica 1**
ULPGC · Grado en Ciencia e Ingeniería de Datos

Simulador en **Python puro y sin ninguna dependencia externa** que
representa de forma integrada los tres procesos del dogma central:

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
| Traducción: lectura de codones, código genético, señal de terminación | `src/proteina.py` → `traducir()` |
| Identificar el papel de moléculas y enzimas en cada proceso | Consola (`visualizador.py`) + fichas interactivas en el informe HTML (`informe_html.py`) |
| Visualización de cada etapa (texto / HTML interactivo) | `src/visualizador.py` (consola y texto) y `src/informe_html.py` (informe navegable con pestañas) |
| Lenguaje e interacción de libre elección | Python + menú de consola (`main.py`) |
| Probar con datos biológicos reales (FASTA / GenBank) | `src/lector.py` + `data/ejemplos/` (gen *lacZ* de *E. coli*) |

## Estructura del proyecto

```
dogma_central/
├── main.py                  # Punto de entrada: menú interactivo por consola
├── src/
│   ├── adn.py                # Clase ADN + simulación de la REPLICACIÓN
│   ├── arn.py                 # Simulación de la TRANSCRIPCIÓN
│   ├── proteina.py            # Código genético + simulación de la TRADUCCIÓN
│   ├── lector.py               # Lectura de FASTA / GenBank / ADN aleatorio
│   ├── visualizador.py         # Construye el informe de texto plano
│   ├── informe_html.py          # Informe HTML interactivo (pestañas, colores, reproducción)
│   └── simulador.py             # Orquesta las tres etapas
├── data/ejemplos/
│   ├── lacZ_Ecoli_NC_000913.3.fasta   # Gen lacZ real (NCBI), formato FASTA
│   └── lacZ_Ecoli_NC_000913.3.gb      # Mismo gen, formato GenBank
└── salidas/                  # Aquí se guardan los informes generados
```

## Cómo ejecutarlo en Visual Studio Code

1. Abre la carpeta `dogma_central/` en VS Code (`Archivo → Abrir carpeta…`).
2. Asegúrate de tener **Python 3.10+** instalado y selecciona el intérprete
   (esquina inferior derecha o `Ctrl+Shift+P` → *Python: Select Interpreter*).
   No hace falta instalar nada más: el proyecto no tiene dependencias.
3. Abre una terminal integrada (`` Ctrl+ñ `` o `Terminal → Nueva terminal`).
4. Ejecuta el programa:
   ```bash
   python main.py
   ```
5. Elige una opción del menú:
   - **1** → escribes tú mismo una secuencia corta de ADN.
   - **2** → se genera una secuencia aleatoria (rápido para probar).
   - **3** → cargas el fichero `.fasta` real de ejemplo (o el que tú añadas
     en `data/ejemplos/`).
   - **4** → igual que la anterior, pero con `.gb` (GenBank).

Al terminar la simulación no se imprime nada por pantalla (solo un
mensaje corto de confirmación); en su lugar puedes guardar:
- un informe en `.txt` con todas las secuencias, los fragmentos de Okazaki
  y la lectura codón a codón,
- y un **informe HTML interactivo** (recomendado): se abre en cualquier
  navegador y tiene una pestaña por etapa (replicación, transcripción,
  traducción), las bases coloreadas por tipo, los fragmentos de Okazaki
  separados y resaltados con su cebador, las enzimas de cada etapa como
  fichas con descripción al pasar el ratón, y un botón para **reproducir la
  traducción codón a codón** viendo crecer la proteína en tiempo real.

Ambos se guardan dentro de `salidas/`.

## Probar con tus propios ficheros

Puedes copiar cualquier fichero `.fasta`/`.fa` o `.gb`/`.gbk` a
`data/ejemplos/` (o darle la ruta completa cuando el programa te la pida) y
el programa detectará automáticamente el formato.

## Notas de diseño

- **Sin ninguna dependencia**: los parsers de FASTA y GenBank están
  escritos a mano (no usan Biopython) y el informe interactivo es HTML/CSS/JS
  autocontenido, por lo que el proyecto funciona con cualquier instalación
  estándar de Python 3, sin `pip install` de por medio.
- **Código genético estándar** (tabla 1 NCBI/IUPAC): válido para fines
  didácticos; se indica en `src/proteina.py` por si se quisiera adaptar a
  otra tabla (por ejemplo, la tabla 11 bacteriana, usada en el `.gb` de
  ejemplo, que apenas difiere de la estándar).
- La traducción comienza en el **primer codón AUG** encontrado y continúa
  hasta el primer codón de parada (UAA, UAG, UGA), tal y como pide el guion.
