import os
import sys
import subprocess
import shutil

def check_dependencies(scripts):
    missing = []
    
    # 1. Check if the scripts for all 3 phases exist in the current folder
    for script in scripts:
        if not os.path.isfile(script):
            missing.append(f"Script not found: {script}")

    # 2. Check if cd-hit is installed (typical in Linux/Ubuntu environments)
    if shutil.which("cd-hit") is None and shutil.which("cdhit") is None:
        missing.append("cd-hit not found in PATH. Install with 'sudo apt install cd-hit' or via conda.")

    if missing:
        print("Error: Missing required dependencies:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)

def select_fasta():
    print("\n--- FASTA File Selection ---")
    
    # Look for files with common FASTA extensions
    fastas = [f for f in os.listdir('.') if f.endswith(('.fasta', '.fa', '.fna', '.faa')) and os.path.isfile(f)]
    
    if not fastas:
        print("Error: No FASTA files found in the current directory.")
        print("Please ensure your FASTA file is in the same folder as this script.")
        sys.exit(1)
    
    for i, f in enumerate(fastas, 1):
        print(f"[{i}] {f}")
    
    while True:
        try:
            choice = int(input("\nSelect the number of the FASTA file you want to use: "))
            if 1 <= choice <= len(fastas):
                selected = fastas[choice - 1]
                print(f"Selected input file: {selected}")
                return selected
            else:
                print("Invalid selection. Please choose a valid number.")
        except ValueError:
            print("Please enter a valid numeric input.")

def configure_phase1():
    # Default parameters for Phase 1
    params = {
        "family": "",
        "matrix": "BLOSUM62",
        "gap_open": -10.0,
        "gap_extend": -0.5,
        "cdhit_thresh": 0.85,
        "threads": 2
    }
    
    while True:
        print("\n--- Phase 1: Parameter Configuration ---")
        print(f"[1] Protein Family Label (--family)       : '{params['family']}'")
        print(f"[2] Substitution Matrix (--matrix)        : {params['matrix']}")
        print(f"[3] Gap Open Penalty (--gap-open)         : {params['gap_open']}")
        print(f"[4] Gap Extend Penalty (--gap-extend)     : {params['gap_extend']}")
        print(f"[5] CD-HIT Threshold (--cdhit-thresh)     : {params['cdhit_thresh']}")
        print(f"[6] CPU Threads (--threads)               : {params['threads']}")
        print("\n[C] Continue to next phase (Accept current values)")
        
        choice = input("\nSelect a parameter to change (or 'C' to continue): ").strip().upper()
        
        if choice == 'C':
            break
        elif choice == '1':
            params['family'] = input("Enter new protein family label: ").strip()
        elif choice == '2':
            mat = input("Enter new matrix (Options: BLOSUM62, BLOSUM50, PAM250): ").strip().upper()
            if mat in ["BLOSUM62", "BLOSUM50", "PAM250"]:
                params['matrix'] = mat
            else:
                print("Invalid matrix. Must be BLOSUM62, BLOSUM50, or PAM250.")
        elif choice == '3':
            try: params['gap_open'] = float(input("Enter new gap open penalty (e.g., -10.0): "))
            except ValueError: print("Invalid numeric input.")
        elif choice == '4':
            try: params['gap_extend'] = float(input("Enter new gap extend penalty (e.g., -0.5): "))
            except ValueError: print("Invalid numeric input.")
        elif choice == '5':
            try: params['cdhit_thresh'] = float(input("Enter new CD-HIT identity threshold (e.g., 0.85): "))
            except ValueError: print("Invalid numeric input.")
        elif choice == '6':
            try: params['threads'] = int(input("Enter number of threads to use (e.g., 2): "))
            except ValueError: print("Invalid numeric input.")
        else:
            print("Invalid choice. Try again.")
            
    return params

def main():
    print("Starting Protein Alignment Pipeline...")
    
    # EXACT FILENAMES ASSIGNED HERE
    script1 = "fase1align_pipeline.py"
    script2 = "phase2_project.py"
    script3 = "fase3codi.py"
    
    # 1. Dependency Check
    scripts = [script1, script2, script3]
    check_dependencies(scripts)
    
    # 2. Interactive Setup
    input_fasta = select_fasta()
    p1_params = configure_phase1()
    
    print("\n--- Phase 2: Parameter Configuration ---")
    print("Note: Phase 2 only requires input/output file paths. These paths are locked for pipeline integrity.")
    input("Press Enter to continue...")

    print("\n--- Phase 3: Parameter Configuration ---")
    print("Note: Phase 3 only requires input/output file paths. These paths are locked for pipeline integrity.")
    input("Press Enter to start the pipeline execution...")

    # ==========================================
    # AUTOMATED DIRECTORY SETUP
    # ==========================================
    out_dir = "pipeline_results"
    os.makedirs(out_dir, exist_ok=True)
    
    msa_out_dir = os.path.join(os.getcwd(), "msa_outputs")
    os.makedirs(msa_out_dir, exist_ok=True)
    
    phase3_out_dir = os.path.join(out_dir, "results_phase3")
    os.makedirs(phase3_out_dir, exist_ok=True)

    output_json = os.path.join(out_dir, "gold_standard.json")
    filtered_fa = os.path.join(out_dir, "filtered.cdhit.fa")
    summary_csv = os.path.join(out_dir, "msa_consistency_summary.csv")

    # ==========================================
    # PHASE 1 EXECUTION
    # ==========================================
    print(f"\n[INFO] Starting Phase 1: CD-HIT Redundancy Filtering ({script1})...")
    cmd_phase1 = [
        sys.executable, script1,
        input_fasta,
        output_json,
        "--filtered-fa", filtered_fa,
        "--family", p1_params["family"],
        "--matrix", p1_params["matrix"],
        "--gap-open", str(p1_params["gap_open"]),
        "--gap-extend", str(p1_params["gap_extend"]),
        "--cdhit-thresh", str(p1_params["cdhit_thresh"]),
        "--threads", str(p1_params["threads"])
    ]
    subprocess.run(cmd_phase1, check=True)

    # ==========================================
    # PHASE 2 EXECUTION
    # ==========================================
    print(f"\n[INFO] Starting Phase 2: MSA Consistency Engine ({script2})...")
    cmd_phase2 = [
        sys.executable, script2,
        "-f", filtered_fa,
        "-j", output_json,
        "-o", summary_csv,
        "-f_orig", input_fasta
    ]
    subprocess.run(cmd_phase2, check=True)

    # ==========================================
    # PHASE 3 EXECUTION
    # ==========================================
    print(f"\n[INFO] Starting Phase 3: Consistency Analysis & Plotting ({script3})...")
    cmd_phase3 = [
        sys.executable, script3,
        "--csv", summary_csv,
        "--tree-dir", "msa_outputs",
        "--fasta-dir", "msa_outputs",
        "--out-dir", phase3_out_dir
    ]
    subprocess.run(cmd_phase3, check=True)

    print("\n[SUCCESS] Pipeline completed.")
    print(f"  - Intermediate results saved in: {out_dir}/")
    print(f"  - MSAs and Trees generated in  : {msa_out_dir}/")
    print(f"  - Final plots and tables in    : {phase3_out_dir}/")

if __name__ == "__main__":
    main()