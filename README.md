1. Descarga y Datos Brutos
descargarfasta.py: El script inicial encargado de automatizar la descarga de las secuencias de la familia de proteínas desde la base de datos.

protein-matching-IPR037093.fasta: El archivo FASTA crudo con todas las secuencias originales correspondientes a la entrada de InterPro IPR037093 que os habéis descargado.

2. Filtrado de Redundancia (Paso CD-HIT)
cdhit.py: Un script de una sola línea que lanza el comando nativo de cd-hit aplicando un umbral de identidad del 90% (-c 0.90) para limpiar el set de datos.

smaug_clean.fasta: Vuestro set de datos optimizado y libre de redundancia generado tras ejecutar CD-HIT (por cierto, ¡gran nombre en clave para el archivo!).

smaug_clean.fasta.clstr: El archivo de salida de clusters que genera CD-HIT automáticamente. Sirve para revisar cómo se han agrupado las secuencias y cuáles se han descartado por parecerse demasiado.

3. Scripts de Alineamiento por Pares (Fase 1)
needlemanwunsch.py: Una implementación propia y manual desde cero en NumPy del algoritmo clásico de Needleman-Wunsch para alineamiento global por pares. Ideal para la parte teórica de la asignatura.

ExtractSeqPairs.py: Un script que automatiza las combinaciones por parejas utilizando itertools.combinations. Coge un subconjunto de prueba de 50 secuencias de smaug_clean.fasta para empezar a calcular alineamientos masivos y exportar los resultados estructurados en un JSON.

Try_better_matrix.py: Un script hermano del anterior que también trabaja con ese bloque de 50 secuencias de prueba, diseñado para experimentar con matrices de sustitución y ajustar los parámetros antes de casaros con un resultado definitivo.

4. Control y Entorno
verificador_alineamientos.py: Un script de control de calidad. Utiliza el PairwiseAligner de Biopython en modo global con puntuaciones planas (match=1, mismatch=0) para hacer pruebas controladas y verificar que los cálculos de puntuación e índices cuadren.

mi_sesion_bioinfo.txt: Un archivo de texto con anotaciones, comandos o el registro de lo que va ocurriendo en tu terminal y tus entornos de Conda mientras trasteas con las librerías.
