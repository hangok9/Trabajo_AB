"""
Phase 2 Pipeline: Consistency of Pairwise Alignments vs. Multiple Sequence Alignments
Author: Pepe (via Expert Bioinformatics Software Engineer)
"""
import argparse
import subprocess
import os
import ijson
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Any

# ==========================================
# PHASE 1: Parameter Sweep & MSA Execution
# ==========================================
def interactive_file_selector(valid_extensions: list, file_type_label: str) -> str:
    """
    Scans the current working directory for files matching specified extensions
    and provides an interactive, numbered selection menu for the user.
    """
    # Scan current directory for files ending with any of the valid extensions
    available_files = [
        f for f in os.listdir('.') 
        if os.path.isfile(f) and any(f.endswith(ext) for ext in valid_extensions)
    ]
    
    if not available_files:
        print(f"\n[!] ERROR: No valid {file_type_label} files found in this directory.")
        print(f"    Expected extensions: {', '.join(valid_extensions)}")
        exit(1)
        
    print(f"\n📂 Multiple {file_type_label} files detected. Please select one:")
    for idx, filename in enumerate(available_files, 1):
        print(f"  [{idx}] {filename}")
        
    while True:
        try:
            choice = int(input(f"Enter the number corresponding to your file (1-{len(available_files)}): "))
            if 1 <= choice <= len(available_files):
                return available_files[choice - 1]
            else:
                print("[!] Out of range selection. Please try again.")
        except ValueError:
            print("[!] Invalid input. Please enter a valid number integer.")

def run_msa_sweep(input_fasta: str, gop: float, gep: float, output_dir: str = ".") -> Tuple[str, str]:
    """
    Executes MAFFT MSA via subprocess with specific Gap Open and Gap Extend penalties.
    
    Args:
        input_fasta: Path to the CD-hit filtered FASTA file.
        gop: Gap Open Penalty.
        gep: Gap Extension Penalty.
        output_dir: Directory to save the resulting MSA FASTA and tree.
        
    Returns:
        A tuple containing:
        - Path to the generated MSA FASTA file.
        - Path to the generated guide tree file (or "None" if not found).
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Mantenim prefixos clars: "msa_" per l'alineament i "tree_" per l'arbre
    base_suffix = f"gop{gop}_gep{gep}"
    output_fasta = os.path.join(output_dir, f"msa_{base_suffix}.fasta")
    output_tree = os.path.join(output_dir, f"tree_{base_suffix}.tree")
    
    # MAFFT command: --op (Gap Open), --ep (Gap Extend/Offset)
    cmd = [
        "mafft",
        "--auto",     # Automatically select the best strategy based on data size
        "--op", str(gop),
        "--ep", str(gep),
        "--treeout",  # Outputs the guide tree
        "--quiet",    # Suppress verbose output to terminal
        input_fasta
    ]

    print(f"[*] Running MSA (GOP: {gop}, GEP: {gep})...")
    with open(output_fasta, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)
    
    expected_mafft_tree = f"{input_fasta}.tree"
    if os.path.exists(expected_mafft_tree):
        os.rename(expected_mafft_tree, output_tree)
    else:
        print(f"[Warning] Guide tree file not found at {expected_mafft_tree}")
        output_tree = "None"

    return output_fasta, output_tree

def read_fasta(filepath: str) -> Dict[str, str]:
    """Lightweight fasta parser to avoid heavy BioPython dependencies."""
    seqs = {}
    with open(filepath, 'r') as f:
        name = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0] # Take ID before first space
                seq = []
            else:
                seq.append(line)
        if name:
            seqs[name] = "".join(seq)
    return seqs

# ==========================================
# PHASE 2: Vectorized Coordinate Indexing
# ==========================================

def build_msa_coordinate_maps(msa_fasta: str) -> Dict[str, np.ndarray]:
    """
    Reads the MSA and creates a vectorized coordinate mapping using np.cumsum().
    Maps absolute biological coordinate (0-indexed) -> MSA column index.
    
    Args:
        msa_fasta: Path to the generated MSA FASTA.
        
    Returns:
        Dictionary mapping sequence ID to a 1D numpy array of absolute MSA columns.
    """
    msa_seqs = read_fasta(msa_fasta)
    msa_maps = {}
    
    for seq_id, seq_str in msa_seqs.items():
        # Convert sequence string to a numpy array of characters
        seq_arr = np.array(list(seq_str))
        
        # Boolean array: True if amino acid, False if gap '-'
        is_aa = (seq_arr != '-')
        
        # Apply np.cumsum() to map biological coordinates.
        # Subtracting 1 makes the coordinates 0-indexed.
        # e.g. is_aa: [True, False, True] -> cumsum - 1: [0, 0, 1]
        cum_aas = np.cumsum(is_aa) - 1
        
        # Initialize the mapping array of size equal to the biological sequence length
        num_aa = np.sum(is_aa)
        mapping = np.zeros(num_aa, dtype=int)
        
        # Map: Biological coordinate -> Current MSA Column Index
        # np.arange(len(seq_str))[is_aa] gives the absolute MSA column indices
        mapping[cum_aas[is_aa]] = np.arange(len(seq_str))[is_aa]
        
        msa_maps[seq_id] = mapping
        
    return msa_maps

# ==========================================
# PHASE 3 & 4: Mapping Engine, RLE, & Stats
# ==========================================

def calculate_rle(delta_array: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Applies Run-Length Encoding (RLE) to an array of error signals (Delta).
    Groups continuous blocks of the same gap shift value.
    
    Args:
        delta_array: 1D numpy array of Delta(i) positional errors.
        
    Returns:
        Tuple of (rle_values, rle_lengths)
    """
    # Find boolean mask of where the array values change
    changes = np.concatenate(([True], delta_array[1:] != delta_array[:-1]))
    
    rle_values = delta_array[changes]
    rle_lengths = np.diff(np.where(np.concatenate((changes, [True])))[0])
    
    return rle_values, rle_lengths


def process_gold_standard_json(
    json_file: str, 
    msa_maps: Dict[str, np.ndarray], 
    original_fasta_dict: Dict[str, str], 
    gop: float, 
    gep: float
) -> List[Dict[str, Any]]:
    """
    Streams the JSON file iteratively (O(1) memory), parses pairwise alignments 
    (supporting both global Needleman-Wunsch and local Smith-Waterman via dynamic offset detection),
    calculates the Delta(i) error signal, compresses via RLE, and aggregates stats on-the-fly.
    """
    results = []
    
    print(f"[*] Streaming JSON and applying Mapping/RLE logic...")
    # 'rb' mode is preferred for ijson high-performance backend parsing
    with open(json_file, 'rb') as f:
        alignments = ijson.kvitems(f, 'alignments')
        
        for pair_id, data in alignments:
            seq_A_id = data['seq_A']
            seq_B_id = data['seq_B']
            nw_seq_A = data['seq_A_aligned']
            nw_seq_B = data['seq_B_aligned']
            identity = data['identity_pct']
            
            # Extract the pure, gapless biological sequence fragments from the JSON alignment
            pure_align_A = nw_seq_A.replace('-', '')
            pure_align_B = nw_seq_B.replace('-', '')
            
            # Dynamic Offset Detection:
            # Returns 0 for global alignments (NW) and the true structural start index for local alignments (SW).
            offset_A = original_fasta_dict[seq_A_id].find(pure_align_A)
            offset_B = original_fasta_dict[seq_B_id].find(pure_align_B)
            
            # Integrity check to ensure alignment segments exist within the reference FASTA sequences
            if offset_A == -1 or offset_B == -1:
                print(f"Warning: Alignment segment not found in original FASTA for pair {pair_id}. Skipping.")
                continue
            
            # Vectorized character analysis of the pairwise alignment strings
            arr_A = np.array(list(nw_seq_A))
            arr_B = np.array(list(nw_seq_B))
            
            is_aa_A = (arr_A != '-')
            is_aa_B = (arr_B != '-')
            
            # Logical mask highlighting positions where BOTH sequences possess aligned amino acids
            match_mask = is_aa_A & is_aa_B
            
            if not np.any(match_mask):
                continue  # Skip processing if no structurally aligned residues exist
            
            # Map local pairwise indices to global coordinates using cumulative sums and the dynamic offsets
            bio_A = (np.cumsum(is_aa_A) - 1) + offset_A
            bio_B = (np.cumsum(is_aa_B) - 1) + offset_B
            
            # Extract global biological coordinates for matching residue pairs
            match_bio_A = bio_A[match_mask]
            match_bio_B = bio_B[match_mask]
            
            # Project biological coordinates directly into MSA column coordinates using Phase 2 vector maps
            try:
                msa_col_A = msa_maps[seq_A_id][match_bio_A]
                msa_col_B = msa_maps[seq_B_id][match_bio_B]
            except KeyError as e:
                print(f"Warning: Sequence {e} not found in MSA mapping indices. Skipping pair {pair_id}.")
                continue
                
            # --- Phase 3: Positional Error Signal & Run-Length Encoding ---
            # Delta(i) measures the column displacement shift between the two sequences inside the MSA matrix
            delta = msa_col_A - msa_col_B
            
            # Compress continuous streams of identical displacement values to maintain O(1) memory efficiency
            rle_values, rle_lengths = calculate_rle(delta)
            
            # --- Phase 4: On-the-fly Statistical Aggregation ---
            total_aligned_residues = len(delta)
            
            # Count how many residue alignments were completely consistent (Delta == 0) using compressed RLE tracking
            consistent_residues = np.sum(rle_lengths[rle_values == 0])
            
            # Calculate final performance summary scores for Phase 3 analytics
            cr_score = consistent_residues / total_aligned_residues
            max_shift = np.max(np.abs(rle_values))
            
            results.append({
                'Pair_ID': pair_id,
                'Seq_A': seq_A_id,
                'Seq_B': seq_B_id,
                'GOP': gop,
                'GEP': gep,
                'Identity_Pct': identity,
                'Total_Aligned': total_aligned_residues,
                'Consistent_Pairs': consistent_residues,
                'CR_Score': cr_score,
                'Max_Shift_Magnitude': max_shift
            })
            
    return results

# ==========================================
# MASTER PIPELINE EXECUTION
# ==========================================

import pandas as pd
# Assegura't d'importar també les teves funcions: run_msa_sweep, build_msa_coordinate_maps, process_gold_standard_json

def main():
    # 1. Initialize Argparse CLI Options
    parser = argparse.ArgumentParser(description="Phase 2: MSA Consistency Engine & Parameter Sweep Pipeline")
    parser.add_argument("-f", "--fasta", type=str, help="Path to input FASTA file (e.g., sequences.cdhit.fa)")
    parser.add_argument("-j", "--json", type=str, help="Path to Gold Standard JSON file (e.g., gold_standard.json)")
    parser.add_argument("-o", "--output", type=str, default="msa_consistency_summary.csv", help="Name of final summary CSV output")
    parser.add_argument("-f_orig", "--fasta_original", type=str, help="Path to original UNFILTERED FASTA file")
    args = parser.parse_args()

    # 2. Resolve Input FASTA (Explicit Argument vs Interactive Fallback)
    fasta_extensions = ['.fasta', '.fa', '.cdhit.fa', '.cdhit.fasta']
    if args.fasta and os.path.exists(args.fasta):
        input_fasta = args.fasta
    else:
        if args.fasta:
            print(f"[!] Target path '{args.fasta}' does not exist. Switching to interactive mode...")
        input_fasta = interactive_file_selector(fasta_extensions, "FASTA")

    # 3. Resolve Input JSON (Explicit Argument vs Interactive Fallback)
    json_extensions = ['.json']
    if args.json and os.path.exists(args.json):
        json_file = args.json
    else:
        if args.json:
            print(f"[!] Target path '{args.json}' does not exist. Switching to interactive mode...")
        json_file = interactive_file_selector(json_extensions, "JSON")
    orig_fasta = args.fasta_original if args.fasta_original and os.path.exists(args.fasta_original) else input_fasta

    output_csv = args.output

    # 4. Display Execution Banner
    print("\n=========================================")
    print(f"🚀 RUNNING PIPELINE: PHASE 2 ENGINE")
    print(f"🔹 FASTA Target: {input_fasta}")
    print(f"🔹 JSON Target:  {json_file}")
    print(f"🔹 CSV Destination: {output_csv}")
    print("=========================================\n")

    # Define Parameter Sweeps
    gop_values = [1.53, 2.0]  # Gap Open Penalties (1.53 is MAFFT default)
    gep_values = [0.0, 0.123] # Gap Extension Penalties
    
    all_metrics = []
    print("[*] Charging FASTA original in memory...")
    original_fasta_dict = read_fasta(orig_fasta)
    for gop in gop_values:
        for gep in gep_values:
            # Execute Phase 1: Unpack MSA FASTA and Newick Guide Tree paths
            msa_fasta, msa_tree = run_msa_sweep(input_fasta, gop, gep, output_dir="msa_outputs")
            
            # Execute Phase 2: Compute Vector Maps from current MSA configuration
            msa_maps = build_msa_coordinate_maps(msa_fasta)
            
            # Execute Phase 3 & 4: Stream and evaluate against Gold Standard JSON
            regime_metrics = process_gold_standard_json(json_file, msa_maps,original_fasta_dict, gop, gep)
            all_metrics.extend(regime_metrics)
            
    # Compile final metrics into a Pandas DataFrame and hand off to Phase 3 (Resi)
    df_results = pd.DataFrame(all_metrics)
    df_results.to_csv(output_csv, index=False)
    print(f"\n[*] Master Pipeline Complete. Summary metrics successfully saved to: {output_csv}")
if __name__ == "__main__":
    main()