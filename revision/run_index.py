"""List the SpaDOT analysis additions made for the manuscript revision.

Run from any directory with Python 3. Use each script's help or header for
required packages and input schemas. Inputs and generated outputs live outside
the repository. Workflows are grouped by manuscript analysis, with one folder
for each set of additions.
"""
from pathlib import Path

REVISION = Path(__file__).resolve().parent
ANALYSES = {
    "encoder_ablations": "Note 7, S39–S40: depth and replacement models, Moran's I and Jacobian sensitivity",
    "loss_ablations": "Table S8 and related panels: architecture and loss controls",
    "kernel_sensitivity": "S47: kernel-scale fits, paper-setting WOT scoring, aggregation and figures",
    "transition_analysis": "Note 8, S41, S18 and Figs. 3c–5c: spot OT, transition evaluation and rendering",
    "single_slice": "Note 9, S42–S46 and Table S9: seven independent-slice workflows and comparisons",
    "stift": "STIFT on the three real datasets and two simulation strategies",
    "real_data_figures": "Figs. 3–5 and related supplement panels: matched predictions, metrics and visualizations",
}
EXISTING = {
    "Harmony comparisons (S30, S34–S37)": "summarize_benchmark/scripts/analyze_Harmony_integrated_latent_*.py",
    "Native vs K-means (S31)": "summarize_benchmark/scripts/*/*_clustering.R",
    "Original transition comparators": "summarize_benchmark/scripts/*/benchmark_OT_transition_table.py",
}


def main():
    print("Manuscript revision analyses")
    for directory, description in ANALYSES.items():
        path = REVISION / directory
        print(f"\n{description}\n  revision/{directory}/")
        for script in sorted(path.rglob("*")):
            if script.suffix in {".py", ".R", ".sh"}:
                print(f"  {script.relative_to(REVISION.parent)}")
    print("\nExisting workflows")
    for label, pattern in EXISTING.items():
        print(f"  {label}")
        print(f"    {pattern}")


if __name__ == "__main__":
    main()
