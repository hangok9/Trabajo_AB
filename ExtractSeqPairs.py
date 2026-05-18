import json
from Bio import SeqIO
from Bio import Align
from Bio.Align import substitution_matrices
from itertools import combinations

# 1. Cargar secuencias (USAMOS 50 PARA LA PRUEBA)
fasta_file = "smaug_clean.fasta"
records = list(SeqIO.parse(fasta_file, "fasta"))[:50] 

# 2. Definir las matrices que queremos probar
matrices_a_probar = ["BLOSUM45", "BLOSUM62", "BLOSUM80"]

# Diccionarios para guardar datos
resultados_json = {}
sumas_identidad = {matriz: 0.0 for matriz in matrices_a_probar} # Para calcular el promedio luego

print(f"Iniciando análisis combinatorio para {len(records)} secuencias...")
print(f"Pares a evaluar: {len(list(combinations(records, 2)))} por matriz.")

# Contador de combinaciones procesadas
total_pares = 0

# 3. Bucle para cada par de secuencias
for record1, record2 in combinations(records, 2):
    total_pares += 1
    pair_name = f"{record1.id}_VS_{record2.id}"
    resultados_json[pair_name] = {}
    
    # 4. Bucle para probar cada matriz en este par específico
    for nombre_matriz in matrices_a_probar:
        # Configuramos el alineador
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        aligner.substitution_matrix = substitution_matrices.load(nombre_matriz)
        
        # Realizamos el alineamiento
        alignments = aligner.align(record1.seq, record2.seq)
        best_alignment = alignments[0]
        
        # Extraemos las dos cadenas alineadas
        seqA_aligned, seqB_aligned = best_alignment[0], best_alignment[1]
        
        # 5. Calcular la Identidad (Matches / Longitud)
        matches = sum(1 for a, b in zip(seqA_aligned, seqB_aligned) if a == b)
        length = len(seqA_aligned)
        identity = (matches / length) * 100
        
        # Sumamos la identidad al total de esta matriz para calcular el promedio final
        sumas_identidad[nombre_matriz] += identity
        
        # 6. Guardar en el diccionario JSON
        resultados_json[pair_name][nombre_matriz] = {
            "score": best_alignment.score,
            "identity_percent": round(identity, 2),
            "seqA_aligned": seqA_aligned,
            "seqB_aligned": seqB_aligned
        }

# 7. Guardar resultados detallados en JSON
with open("resultados_alineamientos_50.json", "w") as f:
    json.dump(resultados_json, f, indent=4)

print("\n=== RESULTADOS FINALES ===")
print("Promedio de identidad de la muestra (50 secuencias):\n")

# 8. Calcular promedios y encontrar la ganadora
promedios = {}
mejor_matriz = ""
mejor_promedio = 0.0

for matriz in matrices_a_probar:
    # El promedio es la suma total dividida entre el número de pares analizados
    promedio_matriz = sumas_identidad[matriz] / total_pares
    promedios[matriz] = promedio_matriz
    print(f"- {matriz}: {round(promedio_matriz, 2)}%")
    
    if promedio_matriz > mejor_promedio:
        mejor_promedio = promedio_matriz
        mejor_matriz = matriz

print(f"\n🏆 LA MATRIZ GANADORA ES: {mejor_matriz} con un {round(mejor_promedio, 2)}% de identidad promedio.")
print("Archivo JSON generado con éxito.")