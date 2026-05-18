#!/usr/bin/env python3

import argparse
import itertools
import json
import os
import subprocess
import time

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices


def parse_args():
    p = argparse.ArgumentParser(
        description="Pairwise protein alignment pipeline with CD-HIT redundancy filtering."
    )
    p.add_argument("input_fasta",  help="Input FASTA file with sequences to align")
    p.add_argument("output_json",  help="Path for the JSON results file")
    p.add_argument("--filtered-fa", default=None,
                   help="Path for the CD-HIT filtered FASTA (default: <input>.cdhit.fa)")
    p.add_argument("--family",     default="", help="Protein family label stored in metadata")
    p.add_argument("--matrix",     default="BLOSUM62", help="Substitution matrix (default: BLOSUM62)") # mirar si es pot posar les 3 opcions 
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


def align_all_pairs(records, matrix_name, gap_open, gap_extend):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score   = gap_open
    aligner.extend_gap_score = gap_extend

    n_pairs = len(records) * (len(records) - 1) // 2
    print(f"Running {n_pairs} pairwise alignments...")

    results = {}
    t0 = time.time()

    for rec_a, rec_b in itertools.combinations(records, 2):
        best = aligner.align(str(rec_a.seq), str(rec_b.seq))[0]

        seq_a_aln, seq_b_aln = extract_aligned_strings(str(best))

        matches = sum(
            a == b and a != "-"
            for a, b in zip(seq_a_aln, seq_b_aln)
        )
        ungapped = sum(
            a != "-" and b != "-"
            for a, b in zip(seq_a_aln, seq_b_aln)
        )
        identity = round(matches / ungapped * 100, 2) if ungapped else 0.0
        gaps = seq_a_aln.count("-") + seq_b_aln.count("-")

        key = f"{rec_a.id}_vs_{rec_b.id}"
        results[key] = {
            "seq_A": rec_a.id,
            "seq_B": rec_b.id,
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
