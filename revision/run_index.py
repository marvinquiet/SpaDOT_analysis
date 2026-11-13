"""List the SpaDOT analysis additions made for the manuscript revision.

Run from any directory with Python 3. Use each script's help or header for
required packages and input schemas. Inputs and generated outputs live outside
the repository. Workflows are grouped by manuscript analysis, with one folder
for each set of additions.
"""
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
ANALYSES = {
    "ablation_studies/encoder_ablations": "Note 7, S39–S40: depth and replacement models, Moran's I and Jacobian sensitivity",
    "ablation_studies/loss_ablations": "Table S8 and related panels: architecture and loss controls",
    "sensitivity_analysis/kernel_sensitivity": "S47: kernel-scale fits, paper-setting WOT scoring, aggregation and figures",
    "sensitivity_analysis/threshold_sensitivity": "Note 8 and S41: transition-score threshold sensitivity",
    "benchmark_study/single_slice_alignment": "Note 9, S42–S46 and Table S9: seven independent-slice workflows and comparisons",
}
METHODS = {
    "Method-specific benchmark runners": "benchmark_study/explore_*/",
}
EXISTING = {
    "Harmony comparisons (S30, S34–S37)": "summarize_benchmark/scripts/analyze_Harmony_integrated_latent_*.py",
    "Native vs K-means (S31)": "summarize_benchmark/scripts/*/*_clustering.R",
    "Real-data metrics and figures (Figs. 3–5)": "summarize_benchmark/scripts/",
    "Spot OT and transition benchmarks (S18 and Figs. 3c–5c)": "summarize_benchmark/scripts/*/benchmark_OT_transition_table.py",
    "Transition rendering and STIFT comparison": "summarize_benchmark/scripts/{render_transitions.R,merge_transition_benchmark.py}",
}


def main():
    print("Manuscript revision analyses")
    for directory, description in ANALYSES.items():
        path = REPOSITORY / directory
        print(f"\n{description}\n  {directory}/")
        for script in sorted(path.rglob("*")):
            if script.suffix in {".py", ".R", ".sh"}:
                print(f"  {script.relative_to(REPOSITORY)}")
    print("\nMethod-specific workflows")
    for description, directory in METHODS.items():
        print(f"  {description}\n    {directory}")
    print("\nExisting workflows")
    for label, pattern in EXISTING.items():
        print(f"  {label}")
        print(f"    {pattern}")


if __name__ == "__main__":
    main()
