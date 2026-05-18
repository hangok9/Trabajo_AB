#!/usr/bin/env python3

import subprocess
import json
import itertools
import time
import os
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

FAMILY_NAME  = "SH3 Domain (PF00018)"
MATRIX_NAME  = "BLOSUM62"

GAP_OPEN     = -10.0
GAP_EXTEND   = -0.5

CDHIT_THRESH = 0.85

INPUT_FASTA  = "sh3_sequences.fasta"
FILTERED_FA  = "sh3_filtered.fasta"
OUTPUT_JSON  = "gold_standard_alignments.json"

SH3_SEQUENCES = {
    "SRC_HUMAN":   "MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEPKLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYIPSNYVAPVDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVCPTSKPQTQGLAKDAWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGETGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREV",

    "GRB2_HUMAN":  "MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKDGFIPKNYIEMKPHPWFFGKIPRAKAIGNDNTGWLVDEVDERHIPIYKHDITDLHNLNKGDTSRVEHIIQPRSPKQLHFLKVSKEPFYKFLHDIQMKKVLDLQNAFLNEVHKLRQEIEQDFKLKLQEIEQDFKLKLQEIEQDFKL",

    "NCK1_HUMAN":  "MASNDFYDSERDNGTYNAPGPPPYVNPFSSAGSSSGAVVGPAPVGGPGFIPLGQPPPQQHPPQQHPPQQHAPPPQHQAPPPHQAPPPHSGSGSGSGSGSGSSAAELQRQTLNMDTTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWW",

    "P85A_HUMAN":  "MSRQSTLYSFFPQTLWPVEHVNWLDELKALMKVNLPDGGSIIVAQYELDIYKNLQQELQKKLQNQQEQKIDEIASGIVQFQQGDGSGSSSGALEEDGESQPKNLQSSQGEADTQHQPPQQHIAEEVASSQNTSSNSTSGAQNPVQQTTHVPQQNPAPQQEPPPPQPQQKPRAAAHQKPQNPAAQPPQHHLHQPPQRPQRAPPQPAQPQPAQHQQPAQHHQQAQPQHPQAQPAQHQQPAQHHQQAQPQHPQAQPAQHQQPAQHHQQAQPQHPQAQPAQHQQPAQHHQQA",

    "LCK_HUMAN":   "MGCVCSSNPEDDWMENIDSMDLETNIASVGQTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYIPSNYVAPVDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVK",

    "SYK_HUMAN":   "MASPDPAAHLPFFYGSISAEEAALERIQNLTQLGSIDNILDTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYIPSNYVAPVDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFD",

    "BTK_HUMAN":   "MTMGGLKGDVQIPEEEGSGLEVLFQGPGSHRGPSAAFAPAAAEPKLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYIPSNYVAPVDSIQAEEWYFGKITRRESERLLL",

    "ITK_HUMAN":   "MNNFPKVHSVQLKFTEDEKFIFAIFQKELQNLISDGDKTLVERDVKHTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYIPSNYVAPVDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAY",
}


def write_fasta(sequences, path):
    with open(path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")

    print(f"Saved {len(sequences)} sequences to {path}")


def run_cdhit(input_fa, output_fa, threshold=0.85):
    #Run CD-HIT to remove highly similar sequences.

    # CD-HIT requires different word sizes depending on identity threshold
    word_size = 5 if threshold >= 0.7 else 4

    cmd = [
        "cd-hit",
        "-i", input_fa,
        "-o", output_fa,
        "-c", str(threshold),
        "-n", str(word_size),
        "-M", "2000",
        "-T", "2"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("CD-HIT failed.")
        print(result.stderr)
        raise RuntimeError("Error while running CD-HIT")

    print("CD-HIT completed successfully.")


def load_sequences(fasta_path):
    return list(SeqIO.parse(fasta_path, "fasta"))


def compute_alignments(records):

    aligner = Align.PairwiseAligner()

    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(MATRIX_NAME)

    aligner.open_gap_score = GAP_OPEN
    aligner.extend_gap_score = GAP_EXTEND

    n = len(records)
    total_pairs = n * (n - 1) // 2

    print(f"Computing {total_pairs} pairwise alignments...")

    results = {}
    start_time = time.time()

    # Generate all unique pairs
    for seqA, seqB in itertools.combinations(records, 2):

        pair_name = f"{seqA.id}_vs_{seqB.id}"

        alignments = aligner.align(str(seqA.seq), str(seqB.seq))
        best = alignments[0]

        # Convert alignment object into aligned strings
        lines = str(best).split("\n")

        seq1_aligned = ""
        seq2_aligned = ""

        for i, line in enumerate(lines):

            parts = line.split()

            if not parts:
                continue

            if i % 4 == 0:
                seq1_aligned += parts[-1]

            elif i % 4 == 2:
                seq2_aligned += parts[-1]

        # Basic alignment statistics
        matches = sum(
            a == b and a != "-"
            for a, b in zip(seq1_aligned, seq2_aligned)
        )

        ungapped_positions = sum(
            a != "-" and b != "-"
            for a, b in zip(seq1_aligned, seq2_aligned)
        )

        identity = (
            round(matches / ungapped_positions * 100, 2)
            if ungapped_positions > 0 else 0
        )

        gaps = (
            seq1_aligned.count("-") +
            seq2_aligned.count("-")
        )

        results[pair_name] = {
            "seq_A": seqA.id,
            "seq_B": seqB.id,

            "score": float(best.score),

            "seq_A_aligned": seq1_aligned,
            "seq_B_aligned": seq2_aligned,

            "alignment_length": len(seq1_aligned),
            "identity_pct": identity,
            "gaps": gaps,

            "algorithm": "Needleman-Wunsch",
            "matrix": MATRIX_NAME,
            "gap_open": GAP_OPEN,
            "gap_extend": GAP_EXTEND
        }

    elapsed = time.time() - start_time

    print(f"Finished {len(results)} alignments in {elapsed:.2f} seconds.")

    return results, elapsed


def save_json(records, results, elapsed, path):
    """
    Save alignment results into JSON format.
    """

    output = {
        "metadata": {
            "family": FAMILY_NAME,
            "algorithm": "Needleman-Wunsch",
            "matrix": MATRIX_NAME,

            "gap_open": GAP_OPEN,
            "gap_extend": GAP_EXTEND,

            "cdhit_threshold": CDHIT_THRESH,

            "n_sequences": len(records),
            "n_pairs": len(results),

            "time_seconds": round(elapsed, 3)
        },

        "alignments": results
    }

    with open(path, "w") as f:
        json.dump(output, f, indent=2)

    size_kb = os.path.getsize(path) / 1024

    print(f"Results saved to {path} ({size_kb:.1f} KB)")


def print_summary(results):
    """
    Print a small summary of the results.
    """

    scores = [v["score"] for v in results.values()]
    identities = [v["identity_pct"] for v in results.values()]

    print("\n Summary")
    print(f"Total alignments: {len(results)}")
    print(f"Average score:    {sum(scores)/len(scores):.2f}")
    print(f"Average identity: {sum(identities)/len(identities):.2f}%")


if __name__ == "__main__":

    print("Starting alignment pipeline...\n")

    # Step 1: Save sequences
    write_fasta(SH3_SEQUENCES, INPUT_FASTA)

    # Step 2: Reduce redundancy with CD-HIT
    run_cdhit(INPUT_FASTA, FILTERED_FA, CDHIT_THRESH)

    # Step 3: Load filtered sequences
    records = load_sequences(FILTERED_FA)

    print(f"Loaded {len(records)} filtered sequences.\n")

    # Step 4: Compute alignments
    results, elapsed = compute_alignments(records)

    # Step 5: Save results
    save_json(records, results, elapsed, OUTPUT_JSON)

    # Step 6: Print summary
    print_summary(results)


# SET 2 
 
import subprocess, json, itertools, os
from Bio import AlignIO
 
# Nombres de los archivos que usamos
INPUT = "sh3_filtered.fasta"       # secuencias de entrada (del Set 1)
MSA   = "sh3_msa.fasta"            # aqui guardaremos el alineamiento multiple
PAIRS = "msa_projected_pairs.json" # aqui guardaremos los pares extraidos
DIFFS = "discrepancy_map.json"     # aqui guardaremos las diferencias con Set 1
GOLD  = "gold_standard_alignments.json"  # resultado del Set 1
  
# Llamamos a MAFFT (programa externo) y guardamos el resultado en MSA
with open(MSA, "w") as f:
    subprocess.run([
        "mafft",
        "--globalpair",
        "--maxiterate", "1000",
        "--treeout",
        "--quiet",
        INPUT
    ], stdout=f) 
#generara el fasta de MSA y el arbol 

#2
 
seqs = list(AlignIO.read(MSA, "fasta"))
 
parells = {}
 
for A, B in itertools.combinations(seqs, 2):
 
    netA, netB = "", ""
 
    # Recorremos columna a columna
    for a, b in zip(str(A.seq), str(B.seq)):
        if a == "-" and b == "-":
            continue  # ELIMINARRRRRR si las DOS tienen gap, lo ignoramos (es un gap falso)
        netA += a
        netB += b
 
    parells[f"{A.id}_vs_{B.id}"] = {
        "seq_A_aligned": netA,
        "seq_B_aligned": netB
    }

json.dump(parells, open(PAIRS, "w"), indent=2)# lode indent=2 solo cambia el formato pa que quede bonito el archivo este json que es para que si pongo el rpint no se peierda todo asi s eguarda el alineamiento en un archivo para q lo reutilice alguien mas o lo q sea 
 
#3 
if not os.path.exists(GOLD):
    print("  No hi ha gold standard. Executa primer el Set 1!")
else:
    gold  = json.load(open(GOLD))["alignments"]
    diffs = {}
 
    for clau, proj in parells.items():
        # Buscamos el mismo par en el gold standard
        ref = gold.get(clau) or gold.get(f"{proj['seq_B']}_vs_{proj['seq_A']}")
        if not ref:
            continue
 
        # Miramos cuantas posiciones estan alineadas igual en los dos metodos
        def mapa(aA, aB):
            # Crea un diccionario: posicion en A → posicion en B
            m, pa, pb = {}, 0, 0
            for x, y in zip(aA, aB):
                if x != "-" and y != "-": m[pa] = pb; pa += 1; pb += 1
                elif x != "-":            m[pa] = None; pa += 1
                elif y != "-":            pb += 1
            return m
 
        mapa_ref  = mapa(ref["seq_A_aligned"],  ref["seq_B_aligned"])
        mapa_proj = mapa(proj["seq_A_aligned"], proj["seq_B_aligned"])
 
        # Fraccion de posiciones que coinciden
        acord = sum(mapa_proj.get(p) == v for p, v in mapa_ref.items())
        acord = round(acord / len(mapa_ref), 4) if mapa_ref else 0
 
        diffs[clau] = {
            "seq_A": proj["seq_A"], "seq_B": proj["seq_B"],
            "identitat_NW":     ref["identity_pct"],
            "identitat_MSA":    proj["identity_pct"],
            "diferencia":       round(ref["identity_pct"] - proj["identity_pct"], 2),
            "acord_columnes":   acord
        }
 
    json.dump(diffs, open(DIFFS, "w"), indent=2)
 
    vals  = list(diffs.values())
    print(f"  Diferencia mitjana d'identitat: {sum(v['diferencia'] for v in vals)/len(vals):.2f}%")
    print(f"  Acord mitjà de columnes:        {sum(v['acord_columnes'] for v in vals)/len(vals):.4f}")
 