#!/usr/bin/env python3

import argparse
import itertools
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices


def parse_args():
    p = argparse.ArgumentParser(
        description="Pairwise protein alignment pipeline with CD-HIT redundancy filtering."
    )
    p.add_argument("input_fasta",  help="Input FASTA file with sequences to align")
    p.add_argument("output_json", default=".", help="Path for the JSON results file")
    p.add_argument("--filtered-fa", default=None,
                   help="Path for the CD-HIT filtered FASTA (default: <input>.cdhit.fa)")
    p.add_argument("--family",     default="", help="Protein family label stored in metadata")
    p.add_argument("--matrix",     default="BLOSUM62", choices=["BLOSUM62","BLOSUM50","PAM250"],help="Substitution matrix (default: BLOSUM62)") # mirar si es pot posar les 3 opcions 
    p.add_argument("--gap-open",   type=float, default=-10.0, help="Gap open penalty (default: -10.0)")
    p.add_argument("--gap-extend", type=float, default=-0.5,  help="Gap extend penalty (default: -0.5)")
    p.add_argument("--cdhit-thresh", type=float, default=0.85,
                   help="CD-HIT identity threshold for clustering (default: 0.85)")
    p.add_argument("--threads",    type=int, default=2, help="Threads passed to CD-HIT (default: 2)")
    return p.parse_args()
def run_cdhit(input_fa, output_fa, threshold, threads):
    word_size = 5 if threshold >= 0.7 else 4
    cmd = [
        "cd-hit",
        "-i", input_fa,
        "-o", output_fa,
        "-c", str(threshold),
        "-n", str(word_size),
        "-M", "2000",
        "-T", str(threads),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("CD-HIT stderr:\n", result.stderr)
        raise RuntimeError("CD-HIT exited with a non-zero status")
    print("CD-HIT finished.")


def load_records(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Loaded {len(records)} sequences from {fasta_path}")
    return records


import itertools
import time
from concurrent.futures import ProcessPoolExecutor
from Bio import Align
from Bio.Align import substitution_matrices

# 1. Definimos el "trabajador" (worker) fuera de la función principal
# Esto es necesario para que el Multiprocesamiento funcione en Python.
def align_single_pair(pair_data):
    id_a, seq_a, id_b, seq_b, matrix_name, gap_open, gap_extend = pair_data
    
    # Instanciamos el alineador localmente dentro del proceso
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    
    # Alineamos usando los strings ya convertidos
    best = aligner.align(seq_a, seq_b)[0]
    
    seq_a_aln, seq_b_aln = extract_aligned_strings(seq_a, seq_b, best)

    # Optimización menor: sumar usando generadores
    matches = sum(1 for a, b in zip(seq_a_aln, seq_b_aln) if a == b and a != "-")
    ungapped = sum(1 for a, b in zip(seq_a_aln, seq_b_aln) if a != "-" and b != "-")
    
    identity = round(matches / ungapped * 100, 2) if ungapped else 0.0
    gaps = seq_a_aln.count("-") + seq_b_aln.count("-")
    
    key = f"{id_a}_vs_{id_b}"
    result_dict = {
        "seq_A": id_a,
        "seq_B": id_b,
        "score": float(best.score),
        "seq_A_aligned": seq_a_aln,
        "seq_B_aligned": seq_b_aln,
        "alignment_length": len(seq_a_aln),
        "identity_pct": identity,
        "gaps": gaps,
        "algorithm": "Needleman-Wunsch",
        "matrix": matrix_name,
        "gap_open": gap_open,
        "gap_extend": gap_extend,
    }
    return key, result_dict

# 2. La función principal refactorizada
def align_all_pairs(records, matrix_name, gap_open, gap_extend, num_cores=None):
    # Si num_cores es None, usará todos los núcleos disponibles de tu ordenador
    n_pairs = len(records) * (len(records) - 1) // 2
    print(f"Running {n_pairs} pairwise alignments using multiprocessing...")

    t0 = time.time()

    # PRE-PROCESAMIENTO: Guardamos en memoria las secuencias ya convertidas a string
    seqs = {rec.id: str(rec.seq) for rec in records}

    # Creamos un generador de tareas (muy ligero en memoria RAM)
    tasks = (
        (rec_a.id, seqs[rec_a.id], rec_b.id, seqs[rec_b.id], matrix_name, gap_open, gap_extend)
        for rec_a, rec_b in itertools.combinations(records, 2)
    )

    results = {}
    
    # MULTIPROCESAMIENTO: Repartimos la carga de trabajo
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        # El chunksize agrupa tareas, reduciendo el cuello de botella entre procesos
        for key, res in executor.map(align_single_pair, tasks, chunksize=500):
            results[key] = res

    elapsed = time.time() - t0
    print(f"Done — {len(results)} alignments in {elapsed:.2f}s")
    return results, elapsed

def extract_aligned_strings(rec_a_seq, rec_b_seq, alignment):
    seq_a_aln = ""
    seq_b_aln = ""
    
    prev_a, prev_b = 0, 0

    for (start_a, end_a), (start_b, end_b) in zip(*alignment.aligned):
        gap_a = start_a - prev_a
        gap_b = start_b - prev_b

        if gap_a > 0:
            seq_a_aln += rec_a_seq[prev_a:start_a]
            seq_b_aln += "-" * gap_a
        elif gap_b > 0:
            seq_a_aln += "-" * gap_b
            seq_b_aln += rec_b_seq[prev_b:start_b]

        seq_a_aln += rec_a_seq[start_a:end_a]
        seq_b_aln += rec_b_seq[start_b:end_b]

        prev_a, prev_b = end_a, end_b

    tail_a = len(rec_a_seq) - prev_a
    tail_b = len(rec_b_seq) - prev_b
    if tail_a > 0:
        seq_a_aln += rec_a_seq[prev_a:]
        seq_b_aln += "-" * tail_a
    elif tail_b > 0:
        seq_a_aln += "-" * tail_b
        seq_b_aln += rec_b_seq[prev_b:]

    return seq_a_aln, seq_b_aln

def save_results(records, results, elapsed, args):
    payload = {
        "metadata": {
            "family": args.family,
            "algorithm": "Needleman-Wunsch",
            "matrix": args.matrix,
            "gap_open": args.gap_open,
            "gap_extend": args.gap_extend,
            "cdhit_threshold": args.cdhit_thresh,
            "n_sequences": len(records),
            "n_pairs": len(results),
            "time_seconds": round(elapsed, 3),
        },
        "alignments": results,
    }
    with open(args.output_json, "w") as fh:
        json.dump(payload, fh, indent=2)
    size_kb = os.path.getsize(args.output_json) / 1024
    print(f"Results written to {args.output_json} ({size_kb:.1f} KB)")


def print_summary(results):
    scores     = [v["score"]        for v in results.values()]
    identities = [v["identity_pct"] for v in results.values()]
    print("\n=== Summary ===")
    print(f"  Alignments : {len(results)}")
    print(f"  Avg score  : {sum(scores)     / len(scores):.2f}")
    print(f"  Avg identity: {sum(identities) / len(identities):.2f}%")


def main():
    args = parse_args()

    filtered_fa = args.filtered_fa or os.path.splitext(args.input_fasta)[0] + ".cdhit.fa"

    print("Step 1 — filtering redundant sequences with CD-HIT")
    run_cdhit(args.input_fasta, filtered_fa, args.cdhit_thresh, args.threads)

    print("\nStep 2 — loading filtered sequences")
    records = load_records(filtered_fa)

    print("\nStep 3 — computing pairwise alignments")
    results, elapsed = align_all_pairs(
        records, args.matrix, args.gap_open, args.gap_extend
    )

    print("\nStep 4 — saving results")
    save_results(records, results, elapsed, args)

    print_summary(results)


if __name__ == "__main__":
    main()
