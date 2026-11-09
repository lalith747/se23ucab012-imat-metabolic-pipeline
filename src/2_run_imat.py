"""
2_run_imat.py  --  iMAT (Zur et al. 2010) on Human-GEM, with a threshold sweep.

Inputs : Human-GEM.xml, processed_<GROUP>_expression.csv  (columns: Gene_ID, Mean_Expression)
Outputs: reaction_results_<GROUP>_<CONFIG>.csv, iMAT_parameter_comparison.csv

Usage  : python 2_run_imat.py            (AD and Control, every config)
         python 2_run_imat.py AD         (one group only)

iMAT in one paragraph
---------------------
Genes are labelled HIGH / LOW / MEDIUM by expression percentile. Each reaction gets a
label from its gene rule (GPR): AND -> min of the gene labels, OR -> max. A MILP then
looks for a flux vector v with  S v = 0  (steady state) that maximises
    (# HIGH reactions with |v| >= EPS)  +  (# LOW reactions with v = 0).
MEDIUM reactions are left unconstrained.

The MILP is solved with HiGHS through scipy.optimize.milp (no extra solver install).
"""
import ast
import sys
import numpy as np
import pandas as pd
import cobra
from cobra.util.array import create_stoichiometric_matrix
from scipy import sparse
from scipy.optimize import milp, LinearConstraint, Bounds

# ------------------------- settings -------------------------
MODEL_FILE = "Human-GEM.xml"
GROUPS = sys.argv[1:] or ["AD", "Control"]
CONFIGS = [                # (name, low quantile, high quantile)
    ("Baseline_25_75", 0.25, 0.75),
    ("Strict_15_85",   0.15, 0.85),
]
EPS = 1.0                  # minimum |flux| for a HIGH reaction to count as active
TIME_LIMIT = 900           # seconds per MILP (HiGHS returns the best solution found so far)
MIP_GAP = 0.02             # stop when within 2 % of the proven upper bound

# ------------------------- load model -------------------------
print("[1] Loading Human-GEM ...")
model = cobra.io.read_sbml_model(MODEL_FILE)
n_rxn = len(model.reactions)
S = create_stoichiometric_matrix(model, array_type="dok").tocsr()   # metabolites x reactions
LB = np.array([r.lower_bound for r in model.reactions], dtype=float)
UB = np.array([r.upper_bound for r in model.reactions], dtype=float)
print(f"  {n_rxn} reactions, {S.shape[0]} metabolites, {len(model.genes)} genes")


# ------------------------- gene -> reaction labels -------------------------
def gene_states_from_expression(expr, low_q, high_q):
    """+1 = high (>= high_q quantile), -1 = low (<= low_q quantile), 0 = medium."""
    lo, hi = expr.quantile(low_q), expr.quantile(high_q)
    states = {g: (1 if v >= hi else -1 if v <= lo else 0) for g, v in expr.items()}
    return states, lo, hi


def reaction_label(node, gene_state):
    """Evaluate a GPR syntax tree on labels: AND = min, OR = max.
    Genes with no measurement count as 0 (medium), so they never force HIGH or LOW."""
    if isinstance(node, ast.Name):
        return gene_state.get(node.id, 0)
    if isinstance(node, ast.BoolOp):
        vals = [reaction_label(v, gene_state) for v in node.values]
        return min(vals) if isinstance(node.op, ast.And) else max(vals)
    if hasattr(node, "body"):              # cobra's GPR wrapper object
        return reaction_label(node.body, gene_state)
    return 0


def classify_reactions(model, gene_state):
    return np.array([
        reaction_label(r.gpr, gene_state) if r.gene_reaction_rule else 0
        for r in model.reactions
    ])


# ------------------------- the iMAT MILP -------------------------
def run_imat(labels):
    """Decision vector x = [v (n_rxn fluxes), y (binary indicators)].

    HIGH, forward : y=1 -> v >= EPS      written as   v - (EPS - lb) * y >= lb
    HIGH, reverse : y=1 -> v <= -EPS     written as   v + (ub + EPS) * y <= ub
    LOW           : y=1 -> v = 0         written as   v + ub * y <= ub   and   v + lb * y >= lb
    (when y=0 each row just reduces to the reaction's own bounds, so it never blocks flux)
    """
    rows, cols, vals = [], [], []         # extra constraint rows (sparse triplets)
    row_lo, row_hi = [], []
    y_info = []                           # (reaction index, kind) for every binary
    r_i = 0

    def add_row(entries, lo, hi):
        nonlocal r_i
        for c, v in entries:
            rows.append(r_i); cols.append(c); vals.append(v)
        row_lo.append(lo); row_hi.append(hi)
        r_i += 1

    for i in np.where(labels == 1)[0]:
        if UB[i] >= EPS:
            y = n_rxn + len(y_info); y_info.append((i, "high_fwd"))
            add_row([(i, 1.0), (y, -(EPS - LB[i]))], LB[i], np.inf)
        if LB[i] <= -EPS:
            y = n_rxn + len(y_info); y_info.append((i, "high_rev"))
            add_row([(i, 1.0), (y, UB[i] + EPS)], -np.inf, UB[i])
    for i in np.where(labels == -1)[0]:
        y = n_rxn + len(y_info); y_info.append((i, "low"))
        if UB[i] > 0:
            add_row([(i, 1.0), (y, UB[i])], -np.inf, UB[i])
        if LB[i] < 0:
            add_row([(i, 1.0), (y, LB[i])], LB[i], np.inf)

    n_y = len(y_info)
    print(f"  HIGH rxns: {(labels == 1).sum()} | LOW rxns: {(labels == -1).sum()} | binaries: {n_y}")

    A_extra = sparse.csr_matrix((vals, (rows, cols)), shape=(r_i, n_rxn + n_y))
    A_ss = sparse.hstack([S, sparse.csr_matrix((S.shape[0], n_y))]).tocsr()   # S v = 0
    A = sparse.vstack([A_ss, A_extra]).tocsr()
    lo = np.concatenate([np.zeros(S.shape[0]), row_lo])
    hi = np.concatenate([np.zeros(S.shape[0]), row_hi])

    c = np.concatenate([np.zeros(n_rxn), -np.ones(n_y)])            # maximise sum(y)
    integrality = np.concatenate([np.zeros(n_rxn), np.ones(n_y)])
    bounds = Bounds(np.concatenate([LB, np.zeros(n_y)]), np.concatenate([UB, np.ones(n_y)]))

    res = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=integrality,
               bounds=bounds,
               options={"time_limit": TIME_LIMIT, "mip_rel_gap": MIP_GAP, "disp": False})
    if res.x is None:
        raise RuntimeError(f"No feasible solution found: {res.message}")
    x = res.x
    y = np.round(x[n_rxn:]).astype(int)
    return x[:n_rxn], y, y_info, res


# ------------------------- run everything -------------------------
rows_out = []
for group in GROUPS:
    expr_df = pd.read_csv(f"processed_{group}_expression.csv", index_col=0)
    expr = pd.to_numeric(expr_df.iloc[:, 0], errors="coerce").dropna()
    expr.index = expr.index.astype(str)
    print(f"\n[2] Group {group}: {len(expr)} genes with expression")

    for name, lq, hq in CONFIGS:
        print(f"\n=== {group} / {name} ===")
        gstate, lo_cut, hi_cut = gene_states_from_expression(expr, lq, hq)
        labels = classify_reactions(model, gstate)
        v, y, y_info, res = run_imat(labels)

        out = pd.DataFrame({
            "Reaction_ID": [r.id for r in model.reactions],
            "Label": labels,
            "Flux": v,
        })
        out["Active"] = out["Flux"].abs() >= EPS * 0.99
        out.to_csv(f"reaction_results_{group}_{name}.csv", index=False)

        high_sat = int(((out.Label == 1) & out.Active).sum())
        low_sat = int(((out.Label == -1) & (out.Flux.abs() < 1e-6)).sum())
        rows_out.append({
            "Group": group, "Config": name,
            "Low_cut": round(lo_cut, 4), "High_cut": round(hi_cut, 4),
            "High_genes": sum(s == 1 for s in gstate.values()),
            "Low_genes": sum(s == -1 for s in gstate.values()),
            "High_rxns": int((labels == 1).sum()), "Low_rxns": int((labels == -1).sum()),
            "High_satisfied": high_sat, "Low_satisfied": low_sat,
            "Objective": int(y.sum()),
            "MIP_gap": None if res.mip_gap is None else round(res.mip_gap, 4),
            "Status": "optimal" if res.status == 0 else "time_limit/best_found",
        })
        print(f"  objective {int(y.sum())} | HIGH satisfied {high_sat} | LOW satisfied {low_sat} "
              f"| status {rows_out[-1]['Status']}")
        pd.DataFrame(rows_out).to_csv("iMAT_parameter_comparison.csv", index=False)

print("\nDone.\n", pd.DataFrame(rows_out).to_string(index=False))