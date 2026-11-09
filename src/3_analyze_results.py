"""
3_analyze_results.py -- compare AD vs Control iMAT models and plot the threshold sweep.

Reads : iMAT_parameter_comparison.csv, reaction_results_<GROUP>_<CONFIG>.csv, Human-GEM.xml
Writes: subsystem_changes_<CONFIG>.csv, iMAT_sweep_summary.png

Design choice (be ready to explain it): iMAT only constrains reactions that got a HIGH or LOW
label. A reaction with a MEDIUM label can carry flux or not in the solver's solution for arbitrary
reasons, so we compare ONLY reactions that were constrained (label != 0) in at least one group.
"""
import numpy as np
import pandas as pd
import cobra
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

CONFIGS = ["Baseline_25_75", "Strict_15_85"]

model = cobra.io.read_sbml_model("Human-GEM.xml")
subsystem = {r.id: (r.subsystem if r.subsystem else "Unassigned") for r in model.reactions}


def bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1)
    return out


for cfg in CONFIGS:
    try:
        ad = pd.read_csv(f"reaction_results_AD_{cfg}.csv")
        ct = pd.read_csv(f"reaction_results_Control_{cfg}.csv")
    except FileNotFoundError:
        print(f"skip {cfg}: result files not found")
        continue

    d = ad.merge(ct, on="Reaction_ID", suffixes=("_AD", "_Ctrl"))
    d["Subsystem"] = d["Reaction_ID"].map(subsystem)
    d = d[(d.Label_AD != 0) | (d.Label_Ctrl != 0)].copy()            # constrained in >= 1 group
    d["Gained_in_AD"] = d.Active_AD & ~d.Active_Ctrl                 # ON in AD, OFF in control
    d["Lost_in_AD"] = ~d.Active_AD & d.Active_Ctrl
    d["Changed"] = d.Gained_in_AD | d.Lost_in_AD
    print(f"\n=== {cfg}: {len(d)} constrained reactions, {int(d.Changed.sum())} differ between AD and Control ===")

    rows = []
    n_ch, n_all = int(d.Changed.sum()), len(d)
    for sub, g in d.groupby("Subsystem"):
        a, b = int(g.Changed.sum()), len(g) - int(g.Changed.sum())
        c, e = n_ch - a, (n_all - len(g)) - (n_ch - a)
        _, p = fisher_exact([[a, b], [c, e]])                       # subsystem vs rest
        rows.append({"Subsystem": sub, "Constrained": len(g), "Changed": a,
                     "Gained_in_AD": int(g.Gained_in_AD.sum()), "Lost_in_AD": int(g.Lost_in_AD.sum()), "p": p})
    res = pd.DataFrame(rows)
    res["q_BH"] = bh(res["p"])
    res = res.sort_values("q_BH")
    res.to_csv(f"subsystem_changes_{cfg}.csv", index=False)
    print(res.head(10).to_string(index=False))

# ------------------ sweep figure (uses columns that 2_run_imat.py really writes) ------------------
df = pd.read_csv("iMAT_parameter_comparison.csv")
labels = df["Group"] + "\n" + df["Config"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
ax[0].bar(labels, df["High_rxns"] + df["Low_rxns"], color="#7570b3")
ax[0].set_ylabel("Reactions constrained by expression (HIGH + LOW)")
ax[0].set_title("How many reactions each threshold constrains")
x = np.arange(len(df))
ax[1].bar(x - 0.2, df["High_satisfied"] / df["High_rxns"], 0.4, label="HIGH carrying flux", color="#d95f02")
ax[1].bar(x + 0.2, df["Low_satisfied"] / df["Low_rxns"], 0.4, label="LOW at zero flux", color="#1b9e77")
ax[1].set_xticks(x); ax[1].set_xticklabels(labels, fontsize=8)
ax[1].set_ylim(0, 1); ax[1].set_ylabel("Fraction satisfied"); ax[1].legend()
ax[1].set_title("How well the model fits the expression calls")
for a in ax:
    a.grid(axis="y", linestyle="--", alpha=0.5)
    plt.setp(a.get_xticklabels(), fontsize=8)
plt.tight_layout()
plt.savefig("iMAT_sweep_summary.png", dpi=300)
print("\nSaved iMAT_sweep_summary.png")