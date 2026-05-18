import random
from Bio import Align

def crear_alineador(gap_open, gap_extend):
    """Crea y devuelve un alineador base sin matrices de sustitución."""
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    # Puntuación básica: match=1, mismatch=0
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    return aligner

def prueba_visual(seq1, seq2, gap_open=-1.0, gap_extend=-0.5):
    """
    Imprime el alineamiento en pantalla para buscar bloques conservados y gaps.
    """
    aligner = crear_alineador(gap_open, gap_extend)
    
    # Extraemos el texto si nos pasan objetos SeqRecord
    s1 = str(seq1.seq) if hasattr(seq1, 'seq') else str(seq1)
    s2 = str(seq2.seq) if hasattr(seq2, 'seq') else str(seq2)
    
    alineamientos = aligner.align(s1, s2)
    mejor_alineamiento = alineamientos[0] # Tomamos el mejor resultado
    
    print("\n" + "="*40)
    print("👀 PRUEBA VISUAL DE ALINEAMIENTO")
    print("="*40)
    print(f"Score obtenido: {mejor_alineamiento.score}")
    print("-" * 40)
    print(mejor_alineamiento)
    print("="*40 + "\n")
    
    return mejor_alineamiento.score

def prueba_significancia(seq1, seq2, gap_open=-1.0, gap_extend=-0.5, iteraciones=100):
    """
    Compara el score real contra el score medio de alinear secuencias barajadas.
    """
    aligner = crear_alineador(gap_open, gap_extend)
    
    s1 = str(seq1.seq) if hasattr(seq1, 'seq') else str(seq1)
    s2 = str(seq2.seq) if hasattr(seq2, 'seq') else str(seq2)
    
    # 1. Calculamos el score biológico real
    score_original = aligner.score(s1, s2)
    
    # 2. Bucle de barajado (shuffling)
    scores_aleatorios = []
    lista_s2 = list(s2) # Convertimos a lista para poder barajar las letras
    
    print(f"⏳ Calculando {iteraciones} alineamientos aleatorios...")
    for _ in range(iteraciones):
        random.shuffle(lista_s2)
        s2_barajada = "".join(lista_s2)
        
        score_falso = aligner.score(s1, s2_barajada)
        scores_aleatorios.append(score_falso)
        
    # 3. Calculamos la media de los scores al azar
    media_aleatoria = sum(scores_aleatorios) / iteraciones
    
    # 4. Mostramos resultados
    print("\n" + "="*40)
    print("🎲 PRUEBA DE SIGNIFICANCIA ESTADÍSTICA")
    print("="*40)
    print(f"Score Biológico Original : {score_original}")
    print(f"Score Aleatorio (Media)  : {media_aleatoria:.2f}")
    
    # Una regla heurística básica: si el score real es mucho mayor que el azar
    if score_original > (media_aleatoria * 2):
        print("\n✅ CONCLUSIÓN: Tus penalizaciones generan scores significativos.")
    else:
        print("\n❌ CONCLUSIÓN: El score biológico no se distingue del azar.")
        print("   -> Sugerencia: Revisa tus valores de gap_open y gap_extend.")
    print("="*40 + "\n")
    
    return score_original, media_aleatoria