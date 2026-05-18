import json

def analitzar_gaps_alineaments(nom_arxiu):
    # 1. Obrir i carregar l'arxiu JSON
    with open(nom_arxiu, 'r') as arxiu:
        dades = json.load(arxiu)
    
    alineaments = dades.get("alignments", {})
    
    total_gaps_globals = 0
    total_caracters_globals = 0
    
    print(f"{'PARELL D\'ALINEAMENT':<30} | {'% GAPS':<10} | {'DETALLS'}")
    print("-" * 75)
    
    # 2. Iterar sobre cada parell de seqüències
    for parell, dades_alineament in alineaments.items():
        
        # Obtenim les seqüències directament
        seq_a = dades_alineament["seq_A_aligned"]
        seq_b = dades_alineament["seq_B_aligned"]
        
        # Opció A: Comptar manualment els guions '-' en ambdues seqüències
        gaps_seq_a = seq_a.count("-")
        gaps_seq_b = seq_b.count("-")
        gaps_totals_parell = gaps_seq_a + gaps_seq_b
        
        # La longitud total d'elements és la suma de les longituds d'ambdues seqüències
        longitud_total = len(seq_a) + len(seq_b)
        
        # 3. Calcular el percentatge per aquest parell
        percentatge_gaps = (gaps_totals_parell / longitud_total) * 100
        
        # Mostrar el resultat
        print(f"{parell:<30} | {percentatge_gaps:>5.2f}%    | {gaps_totals_parell} gaps de {longitud_total} posicions")
        
        # 4. Sumar als comptadors globals per fer l'estadística final
        total_gaps_globals += gaps_totals_parell
        total_caracters_globals += longitud_total
        
    # 5. Calcular i mostrar el percentatge global
    if total_caracters_globals > 0:
        percentatge_global = (total_gaps_globals / total_caracters_globals) * 100
        print("-" * 75)
        print(f"PERCENTATGE DE GAPS GLOBAL: {percentatge_global:.2f}% ({total_gaps_globals} gaps en total)")

# Cridem la funció amb el nom del teu arxiu
analitzar_gaps_alineaments("gold_standard_alignments.json")