#!/usr/bin/env python3
"""
Phase 3 — MSA Consistency Analysis
===================================
Project: Algorithms in Biology — Evaluating NW vs. MAFFT consistency
         across parameter regimes.

Usage:
    python phase3_analysis.py \
        --csv  msa_consistency_summary.csv \
        --tree-dir msa_outputs \
        --fasta-dir msa_outputs \
        --out-dir  results_phase3

Outputs (saved to --out-dir):
    block1_length_bias.png
    block2a_twilight_zone.png
    block2b_anova.txt
    block3_guide_tree_effect.png
    block4_error_typology.png
    block5_pareto_gappiness.png
    phase3_summary_table.csv
"""

import argparse
import glob
import os
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")      
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import pearsonr
from statsmodels.formula.api import ols
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.anova import anova_lm
from Bio import Phylo

def extract_base_id(raw_name):
    '''
    we clean the strings with strips. then we search for the pure ID, the regex ignores numeric prefixes
    for ex. from  "12_PROT_A|" we want to extract "PROT"
    '''
    clean_name = str(raw_name).strip(" '\"")
    match = re.search(r'^(?:\d+_)?([A-Z0-9]+)[_|]', clean_name)
    if match:
        return match.group(1)
    return clean_name

def parse_args():
    #here we are setting up CLI args so we can run this from the terminal.
    parser = argparse.ArgumentParser(
        description="Phase 3: MSA Consistency Analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # main input file (the outputs we got from Phase 2)
    parser.add_argument("--csv", default="msa_consistency_summary.csv",
                        help="Path to the Phase 2 summary CSV")    
    # directories where we'll read the trees and alignments from           
    parser.add_argument("--tree-dir", default="msa_outputs",
                        help="Folder containing the guide trees (.tree files)")
    parser.add_argument("--fasta-dir", default="msa_outputs",
                        help="Folder containing the .fasta alignments")
    #output of the figures and tables
    parser.add_argument("--out-dir", default="results_phase3",
                        help="Output folder for the figures and tables")     
    return parser.parse_args()

def setup_plots():
    #here we use the seaborn library
    sns.set_theme(style="whitegrid", font_scale=1.1)
    # we keep the configuration in a quick dictionary for later
    return {
        "dpi": 150,       
        "palette": "Set2" 
    }

def save_plot(fig, path, dpi=150):
    # we use tight layout so labels don't get cut off at the edges
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[Saved] {path}")

def load_fasta(filepath):
    # we make a custom parser to avoid loading heavy libraries for this
    seqs = {}
    current_name = None
    sequence_chunks = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # we use this to save the previous sequence if we have one before moving to the next
                if current_name:
                    seqs[current_name] = "".join(sequence_chunks)
                # we get the new sequence header (just the first token, drop the rest)
                current_name = line[1:].split()[0]
                sequence_chunks = []
            else:
                sequence_chunks.append(line)
        # here we grab the final sequence at the end of the file
        if current_name:
            seqs[current_name] = "".join(sequence_chunks)
    return seqs


def load_csv_data(csv_path):
    df = pd.read_csv(csv_path)
    # we create a nice label combining GOP and GEP for grouping the plots later
    # we use vectorized string concatenation, because it is way faster than using .apply(axis=1)
    df["Regime"] = "GOP=" + df["GOP"].astype(str) + " / GEP=" + df["GEP"].astype(str)
    print(f"Loaded {len(df)} records across {df['Regime'].nunique()} regimes.")
    return df

def check_length_bias(df, out_dir, cfg):
    # we check if longer alignments just score differently by pure chance.
    # we decided that if r < 0.15 and p > 0.05, we're mostly safe to ignore length as a confounder
    print("\n--- Block 1: Length Bias Check ---")
    #we calculate correlation coeff and p-value
    r, p = pearsonr(df["Total_Aligned"], df["CR_Score"])
    print(f"Pearson r = {r:.3f}, p = {p:.3e}")
    #as we said earlier if abs(r) < 0.15 and p > 0.05 the correlation is weak  
    #and if else then it is significant
    if abs(r) < 0.15 and p > 0.05:
        print("[OK] Correlation is weak/non-significant. Length isn't a major confounder.")
    else:
        print("[WARNING] Significant correlation found. Might need to watch out for length bias later.")
    #we make the plot
    fig, ax = plt.subplots(figsize=(7, 5))
    # the base scatter plot
    sns.scatterplot(
        data=df, x="Total_Aligned", y="CR_Score",
        alpha=0.3, s=20, color="steelblue", ax=ax, label="Pair"
    )
    # we make the trendline overlayed
    sns.regplot(
        data=df, x="Total_Aligned", y="CR_Score",
        scatter=False, color="firebrick", line_kws={"linewidth": 2},
        ax=ax, label=f"OLS fit (r={r:.2f}, p={p:.2e})"
    )
    #we set the axis titles, the tile and the legend
    ax.set_xlabel("Alignment Length (Total Aligned Residue Pairs)")
    ax.set_ylabel("Consistency Ratio Score (CR_Score)")
    ax.set_title("Length Bias Check: Total_Aligned vs CR_Score")
    ax.legend()
    
    # we save it using the helper that we made earlier
    save_plot(fig, os.path.join(out_dir, "block1_length_bias.png"), cfg["dpi"])

def analyze_twilight_zone(df, out_dir, cfg):
    # here we analyze how consistency drops when the sequences hit the "twilight zone", 
    # that occurs when there is a 20-35% identity.
    # we will also run an ANOVA to check if our GOP/GEP tweaks actually make a statistical difference.
    print("\n--- Block 2: Twilight Zone & Parameter Sweep ---")
    #scatter plot with the LOESS smoothing
    regimes = df["Regime"].unique()
    palette = sns.color_palette(cfg["palette"], n_colors=len(regimes))
    regime_colors = dict(zip(regimes, palette))

    fig, ax = plt.subplots(figsize=(9, 6))
    
    # we highlight the twilight zone so it's obvious on the chart
    ax.axvspan(20, 35, alpha=0.15, color="gold", label="Twilight Zone (20-35%)")

    for regime, color in regime_colors.items():
        # here we drop NaNs just in case, otherwise the lowess smoother might crash
        sub = df[df["Regime"] == regime].dropna(subset=["Identity_Pct", "CR_Score"])
        
        ax.scatter(sub["Identity_Pct"], sub["CR_Score"], color=color, alpha=0.15, s=15)
        
        # we add a trendline (we found that frac=0.3 is normally a decent baseline for smoothing)
        sm_xy = lowess(sub["CR_Score"].values, sub["Identity_Pct"].values, frac=0.3, return_sorted=True)
        ax.plot(sm_xy[:, 0], sm_xy[:, 1], color=color, linewidth=2, label=regime)

    ax.set_xlabel("Pairwise Sequence Identity (%)")
    ax.set_ylabel("Consistency Ratio Score (CR_Score)")
    ax.set_title("Identity vs Consistency across MAFFT Regimes")
    ax.legend(title="Regime (LOESS trend)", framealpha=0.9)
    
    save_plot(fig, os.path.join(out_dir, "block2_twilight_zone.png"), cfg["dpi"])

    # here we are going to do the two-way anova
    # we treat the gap penalties as categories to see their main and interaction effects
    df_anova = df.copy()
    df_anova["GOP_cat"] = df_anova["GOP"].astype("category")
    df_anova["GEP_cat"] = df_anova["GEP"].astype("category")

    model = ols("CR_Score ~ C(GOP_cat) + C(GEP_cat) + C(GOP_cat):C(GEP_cat)", data=df_anova).fit()
    anova_tbl = anova_lm(model, typ=2)

    # we format the output string 
    anova_out = (
        "--- Two-Way ANOVA: GOP & GEP effects on CR_Score ---\n"
        f"{anova_tbl.to_string()}\n\n"
        "Quick interpretation guide:\n"
        "- p < 0.05 for GOP_cat: GOP changes consistency\n"
        "- p < 0.05 for GEP_cat: GEP changes consistency\n"
        "- p < 0.05 for interaction: The effect of GOP depends on the GEP used\n"
    )
    
    print(anova_out)

    anova_path = os.path.join(out_dir, "block2_anova.txt")
    with open(anova_path, "w") as f:
        f.write(anova_out)
        
    print(f"[Saved] {anova_path}")

def parse_tree_distances(tree_path):
    # we read the newick guide tree using Biopython
    tree = Phylo.read(tree_path, "newick")
    leaves = tree.get_terminals()
    distances = {}
    # we use a loop through every possible pair of leaves to get their branch distance
    for i, leaf_a in enumerate(leaves):
        for leaf_b in leaves[i + 1:]:
            dist = tree.distance(leaf_a, leaf_b)
            # we clean up the names using our earlier function
            id_a = extract_base_id(leaf_a.name)
            id_b = extract_base_id(leaf_b.name)
            # we use a frozenset for the key so we don't have to care about the order
            # (for examples looking up A-B will give the same result as B-A)
            pair_key = frozenset([id_a, id_b])
            distances[pair_key] = dist
    return distances


def add_cophenetic_distance(df, tree_dir):
    # we process the dataframe by pieces (every combination of GOP and GEP) 
    # to be able to assign its corresponding tree
    chunks = []
    for (gop, gep), group in df.groupby(["GOP", "GEP"]):
        # we make a copy because we got the  "SettingWithCopyWarning" error of Pandas
        group = group.copy()
        tree_file = os.path.join(tree_dir, f"tree_gop{gop}_gep{gep}.tree")
        # if the tree is missing, me fill the Nans and we go to the next
        if not os.path.exists(tree_file):
            print(f"[Warning] The tree is not found: {tree_file} - filling with NaNs.")
            group["Cophenetic_Dist"] = np.nan
            chunks.append(group)
            continue
            
        print(f"Reading the tree: {os.path.basename(tree_file)}")
        dist_map = parse_tree_distances(tree_file) # we use the function from above
        
        # this function is aplied to every road. 
        def fetch_dist(row):
            id_a = extract_base_id(row["Seq_A"])
            id_b = extract_base_id(row["Seq_B"])
            
            #we search the ppair in the dictionary, if its not found we return NaN
            return dist_map.get(frozenset([id_a, id_b]), np.nan)
            
        group["Cophenetic_Dist"] = group.apply(fetch_dist, axis=1)
        chunks.append(group)
        
    # we paste again every piece into a single dataframe
    return pd.concat(chunks, ignore_index=True)

def analyze_guide_tree_effect(df, tree_dir, out_dir, cfg):
    # we want to check the "once a gap, always a gap" hypotesis
    # the more distant sequences are aligned later and they drag previous mistakes
    # if the theory is true, we should see a negative slope
    print("\n--- Block 3: Guide Tree Effect ---")
    # we add the distances usigng the function that we made earlier
    df = add_cophenetic_distance(df, tree_dir)
    # we erase Nans so that the OLS stats model can work
    df_clean = df.dropna(subset=["Cophenetic_Dist", "CR_Score"])
    # we make a lineal model (OLS) to extract p-values r^2
    X = sm.add_constant(df_clean["Cophenetic_Dist"])
    model = sm.OLS(df_clean["CR_Score"], X).fit()
    print(model.summary())

    # we assign the colors again
    regimes = df_clean["Regime"].unique()  
    palette = sns.color_palette(cfg["palette"], n_colors=len(regimes))
    regime_colors = dict(zip(regimes, palette))

    fig, ax = plt.subplots(figsize=(8, 5))

    for regime, color in regime_colors.items():
        sub = df_clean[df_clean["Regime"] == regime]
        ax.scatter(sub["Cophenetic_Dist"], sub["CR_Score"],
                   color=color, alpha=0.3, s=15, label=regime)

    # here we create, the regression line
    x_min = df_clean["Cophenetic_Dist"].min()
    x_max = df_clean["Cophenetic_Dist"].max()
    x_rng = np.linspace(x_min, x_max, 100)
    
    # the equation y = a + bx
    y_hat = model.params["const"] + model.params["Cophenetic_Dist"] * x_rng
    
    label_text = f"OLS R²={model.rsquared:.2f}, p={model.pvalues['Cophenetic_Dist']:.2e}"
    ax.plot(x_rng, y_hat, color="black", linewidth=2, label=label_text)

    ax.set_xlabel("Cophenetic Distance in Guide Tree")
    ax.set_ylabel("Consistency Ratio Score (CR_Score)")
    ax.set_title("Guide Tree Distance vs CR_Score\n('Once a gap, always a gap' effect)")
    ax.legend(framealpha=0.9, fontsize=9)
    
    save_plot(fig, os.path.join(out_dir, "block3_guide_tree_effect.png"), cfg["dpi"])

    return df 


def plot_error_typology(df, out_dir, cfg):
    # we want to see if the alignment errors are just micro jitter (<= 2 columns)
    # or massive structural macro misalignments (> 2 columns)
    # we saw that the KDE plot works best here to visualize the distribution
    print("\n--- Block 4: Error Typology ---")

    regimes = df["Regime"].unique()
    palette = sns.color_palette(cfg["palette"], n_colors=len(regimes))
    regime_colors = dict(zip(regimes, palette))

    fig, ax = plt.subplots(figsize=(9, 5))

    for regime, color in regime_colors.items():
        # dropna is important here, otherwise sns.kdeplot will throw a fit
        data = df[df["Regime"] == regime]["Max_Shift_Magnitude"].dropna()
        
        sns.kdeplot(
            data, ax=ax, color=color, fill=True, alpha=0.25, 
            linewidth=2, label=regime, 
            bw_adjust=0.8 # we lower the bandwidth slightly to spot bimodality better
        )

    # draw a line to separate micro from macro errors
    ax.axvline(x=2, color="dimgray", linestyle="--", linewidth=1.5, label="Micro/Macro cutoff (2 cols)")
    
    # here we add a text label right on the plot
    # we grab the y-limit so it stays near the top no matter the data
    y_top = ax.get_ylim()[1]
    ax.text(2.2, y_top * 0.85, "Micro -> Macro", color="dimgray", fontsize=10)

    ax.set_xlabel("Max Shift Magnitude (MSA columns)")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of Max Shift Magnitude by Regime")
    ax.legend(title="Regime", framealpha=0.9)
    
    save_plot(fig, os.path.join(out_dir, "block4_error_typology.png"), cfg["dpi"])

    # we make summary stats to print in the console
    summary = df.groupby("Regime")["Max_Shift_Magnitude"].describe(
        percentiles=[0.25, 0.5, 0.75, 0.90]
    ).round(3)
    
    print("\nMax Shift Magnitude summary:")
    print(summary.to_string())


def compute_gappiness(msa_fasta_path):
    # we need to figure out if the alignment is cheating by over-gapping
    # we have to check two things:
    # 1. the length ratio: hoow much longer is the MSA compared to the raw unaligned sequences?
    # 2. High-gap columns: what fraction of the columns are mostly (>50%) gaps?
    
    seqs = load_fasta(msa_fasta_path) # we use our custom parser from earlier
    ids = list(seqs.keys())
    
    if not ids:
        return None

    # since it's an alignment, the length of the first sequence is the length of the whole MSA
    msa_len = len(seqs[ids[0]])
    
    # we calculate the original lengths by stripping out the gap characters
    raw_lengths = [len(s.replace("-", "")) for s in seqs.values()]
    mean_raw = np.mean(raw_lengths)

    # here we build a boolean matrix where True means the character is a gap.
    # we realized that doing this with numpy is way faster for large alignments
    gap_matrix = np.array([[char == "-" for char in seq] for seq in seqs.values()])
    
    # we get the fraction of gaps for each column (axis=0)
    col_gap_frac = gap_matrix.mean(axis=0)
    
    # we find what percentage of columns consist of more than 50% gaps
    high_gap_frac = float(np.mean(col_gap_frac > 0.5))

    # we get the filename and extract GOP and GEP using a quick regex
    fname = os.path.basename(msa_fasta_path)
    match = re.search(r"gop([\d.]+)_gep([\d.]+)\.fasta", fname)

    gop = float(match.group(1)) if match else float("nan")
    gep = float(match.group(2)) if match else float("nan")

    # here we group everything into a dictionary and return
    return {
        "GOP": gop,
        "GEP": gep,
        "Regime": f"GOP={gop} / GEP={gep}",
        "Length_Ratio": msa_len / mean_raw,
        "High_Gap_Frac": high_gap_frac,
        "N_Seqs": len(ids),
        "MSA_Length": msa_len,
        "Mean_Raw_Length": mean_raw
    }

def analyze_gappiness_tradeoff(df, fasta_dir, out_dir, cfg):
    # we need to see if a high CR_Score is actually good, or if the alignment
    # just "cheated" by inserting a million gaps (over-gapping).
    # top-left quadrant = true winner (high consistency, low gaps).
    # top-right quadrant = suspect (high consistency, but too many gaps).
    print("\n--- Block 5: Gappy-ness vs Consistency Trade-off ---")

    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "msa_gop*.fasta")))
    print(f"Found {len(fasta_files)} MSA FASTA file(s).")

    # we compute the metrics for each file using our helper function
    records = [compute_gappiness(f) for f in fasta_files]
    records = [r for r in records if r is not None]

    if not records:
        print("[Warning] No MSA FASTA files found - skipping Block 5.")
        return

    df_gappy = pd.DataFrame(records)

    # we calculate the mean CR score for each regime
    mean_cr = df.groupby(["GOP", "GEP"])["CR_Score"].mean().reset_index()
    mean_cr.rename(columns={"CR_Score": "Mean_CR_Score"}, inplace=True)
    
    # we merge the gappiness data with the CR scores
    df_pareto = df_gappy.merge(mean_cr, on=["GOP", "GEP"], how="left")

    print("\nGappy-ness + Consistency per regime:")
    print(df_pareto[["Regime", "Length_Ratio", "High_Gap_Frac", "Mean_CR_Score"]].to_string(index=False))

    # save the raw data for the pareto chart
    pareto_path = os.path.join(out_dir, "block5_pareto_table.csv")
    df_pareto.to_csv(pareto_path, index=False)
    print(f"[Saved] {pareto_path}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    pareto_palette = sns.color_palette(cfg["palette"], n_colors=len(df_pareto))

    x_specs = [
        ("Length_Ratio", r"Length Ratio $L_{msa}\,/\,\bar{L}_{seq}$"),
        ("High_Gap_Frac", "Fraction of Columns with >50% Gaps"),
    ]

    for ax, (x_col, x_label) in zip(axes, x_specs):
        for idx, row in df_pareto.iterrows():
            ax.scatter(
                row[x_col], row["Mean_CR_Score"],
                color=pareto_palette[int(idx)],
                s=160, zorder=5, edgecolors="white", linewidth=0.8,
            )
            
            # we have to label each dot with its regime
            ax.annotate(
                row["Regime"],
                xy=(row[x_col], row["Mean_CR_Score"]),
                xytext=(6, 4), textcoords="offset points",
                fontsize=8.5, color="dimgray",
            )

        # we draw dashed lines at the medians to create the 4 quadrants
        ax.axhline(df_pareto["Mean_CR_Score"].median(), color="steelblue", linestyle=":", label="Median CR_Score")
        ax.axvline(df_pareto[x_col].median(), color="tomato", linestyle=":", label=f"Median {x_col}")

        ax.set_xlabel(x_label)
        ax.set_ylabel("Mean CR_Score")
        ax.set_title(f"Gappy-ness vs Consistency ({x_col})")
        ax.legend(fontsize=8, framealpha=0.9)

    fig.suptitle("Pareto Trade-off: Gappy-ness vs Mean Consistency by Regime", fontsize=13, y=1.02)
    save_plot(fig, os.path.join(out_dir, "block5_pareto_gappiness.png"), cfg["dpi"])

    # we made a quick cheat sheet for when we run the archive
    print("\nQuick interpretation guide:")
    print("- Top-left (high CR, low gappy)     -> Genuinely superior regime [WINNER]")
    print("- Top-right (high CR, high gappy)   -> Artificially inflated by over-gapping [SUSPECT]")
    print("- Bottom-left (low CR, low gappy)   -> Clean but inconsistent")
    print("- Bottom-right (low CR, high gappy) -> Messy AND inconsistent [WORST]")



def export_summary_table(df, out_dir):
    # we aggregate all the key metrics per regime into a final CSV
    # so we have all the numbers in one place without having to re-run everything.
    
    summary = df.groupby("Regime").agg(
        N_Pairs=("CR_Score", "count"),
        Mean_CR_Score=("CR_Score", "mean"),
        Median_CR_Score=("CR_Score", "median"),
        Std_CR_Score=("CR_Score", "std"),
        Mean_Max_Shift=("Max_Shift_Magnitude", "mean"),
        Mean_Cophenetic=("Cophenetic_Dist", "mean")
    ).round(4).reset_index()
    
    path = os.path.join(out_dir, "phase3_summary_table.csv")
    summary.to_csv(path, index=False)
    
    # we print it to the console
    print("\n--- Final Summary Table ---")
    print(summary.to_string(index=False))
    print(f"\n[Saved] {path}")

def main():
    # we load the parameters and we make sure that the folder for the output does exist 
    args = parse_args()
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    # we prepare the plots
    cfg = setup_plots()

    # initial log to make sure that everything is fine
    print("=======================================================")
    print(" Phase 3: MSA Consistency Analysis")
    print(f" Input CSV  : {args.csv}")
    print(f" Tree dir   : {args.tree_dir}")
    print(f" FASTA dir  : {args.fasta_dir}")
    print(f" Output dir : {args.out_dir}")
    print("=======================================================")

    # we load the base data
    df = load_csv_data(args.csv)

    # we execute the analysis 
    check_length_bias(df, args.out_dir, cfg)
    analyze_twilight_zone(df, args.out_dir, cfg)
    
    df = analyze_guide_tree_effect(df, args.tree_dir, args.out_dir, cfg)
    
    plot_error_typology(df, args.out_dir, cfg)
    analyze_gappiness_tradeoff(df, args.fasta_dir, args.out_dir, cfg)
    
    # we save the final summary
    export_summary_table(df, args.out_dir)

    print("\n[Done] Phase 3 complete. All outputs saved to:", args.out_dir)


if __name__ == "__main__":
    main()