"""
1_prepare_data.py -- GSE33000 (Agilent two-colour, platform GPL4372) -> per-group Human-GEM expression profiles.

Output files (one per diagnosis group):
    processed_AD_expression.csv       310 Alzheimer's samples
    processed_Control_expression.csv  157 non-demented controls
    processed_HD_expression.csv       157 Huntington's samples (not used by iMAT here; kept for completeness)
Each has columns Gene_ID (Ensembl, matches Human-GEM) and Mean_Expression.

IMPORTANT about the values: this is a two-colour array. Channel 1 is a common reference pool,
channel 2 is the individual brain sample, and the matrix holds log10(sample / pool). So values are
RELATIVE to the pool (centred near 0), not absolute expression. That is fine for iMAT because it only
uses percentile cut-offs, but say so in your README.

The disease label lives in !Sample_characteristics_ch2 ("disease status: ..."), NOT in channel 1.
"""
import gzip
import io
import os
import urllib.request
import pandas as pd

MATRIX = "GSE33000_series_matrix.txt.gz"
MODEL_FILE = "Human-GEM.xml"
PROBE_MAP_CACHE = "probe_to_ensembl.csv"       # built once, reused on later runs

print("=" * 60)
print("GSE33000 -> per-group Human-GEM expression profiles")
print("=" * 60)

# ---------------------------------------------------------
# 1. SAMPLE LABELS (from channel 2 characteristics)
# ---------------------------------------------------------
print("\n[1] Reading sample labels...")
gsm_ids, status = [], []
with gzip.open(MATRIX, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            gsm_ids = [x.strip('" \t\r\n') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch2"):
            vals = [x.strip('" \t\r\n') for x in line.split("\t")[1:]]
            if vals and vals[0].lower().startswith("disease status"):
                status = [v.split(":", 1)[1].strip().lower() for v in vals]
        elif line.startswith("!series_matrix_table_begin"):
            break

if not gsm_ids or not status:
    raise SystemExit("Could not find 'disease status' in !Sample_characteristics_ch2.")

def to_group(s):
    if "alzheimer" in s:
        return "AD"
    if "huntington" in s:
        return "HD"
    if "non-demented" in s or "control" in s:
        return "Control"
    return "Other"

groups = {}
for gsm, s in zip(gsm_ids, status):
    groups.setdefault(to_group(s), []).append(gsm)
print({g: len(v) for g, v in groups.items()})
assert len(groups.get("AD", [])) > 0 and len(groups.get("Control", [])) > 0

# ---------------------------------------------------------
# 2. EXPRESSION MATRIX
# ---------------------------------------------------------
print("\n[2] Parsing expression matrix...")
lines, in_table = [], False
with gzip.open(MATRIX, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            lines.append(line)
df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", quotechar='"')
df.rename(columns={df.columns[0]: "ID_REF"}, inplace=True)
df["ID_REF"] = df["ID_REF"].astype(str).str.strip()
print(f"  {len(df)} probes x {df.shape[1] - 1} samples")

# ---------------------------------------------------------
# 3. PROBE -> HUMAN-GEM ENSEMBL ID MAP (cached)
# ---------------------------------------------------------
def build_probe_map():
    import cobra
    print("\n[3] Building probe -> Ensembl map (first run only)...")
    model = cobra.io.read_sbml_model(MODEL_FILE)
    lookup = {}
    for g in model.genes:
        ens = g.id.split(".")[0]
        lookup[ens.upper()] = ens
        if getattr(g, "name", None):
            lookup[g.name.upper()] = ens

    BEGIN = ("!annotation_table_begin", "!platform_table_begin")
    END = ("!annotation_table_end", "!platform_table_end")
    URLS = [
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/annot/GPL4372.annot.gz",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL4372&targ=self&form=text&view=full",
    ]

    def table_from(text):
        out, on = [], False
        for ln in text.splitlines(keepends=True):
            if ln.startswith(BEGIN):
                on = True
                continue
            if ln.startswith(END):
                break
            if on:
                out.append(ln)
        return out

    def to_text(raw):
        return (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode("utf-8", "ignore")

    gpl = None
    for p in ("GPL4372.annot.gz", "GPL4372.annot", "GPL4372.txt"):
        if os.path.exists(p):
            with open(p, "rb") as fh:
                t = table_from(to_text(fh.read()))
            if t:
                gpl = pd.read_csv(io.StringIO("".join(t)), sep="\t", low_memory=False, on_bad_lines="skip")
                break
    if gpl is None:
        for url in URLS:
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=90) as r:
                    t = table_from(to_text(r.read()))
                if t:
                    gpl = pd.read_csv(io.StringIO("".join(t)), sep="\t", low_memory=False, on_bad_lines="skip")
                    break
            except Exception as e:
                print("  download failed:", e)
    if gpl is None:
        raise SystemExit("Could not load GPL4372 table. Save it as GPL4372.txt next to this script.")
    print("  platform columns:", list(gpl.columns))

    def tok(x):
        s = str(x).split("///")[0].strip()
        return "" if s.lower() == "nan" else s

    id_col = gpl.columns[0]
    sym_cols = [c for c in gpl.columns if "symbol" in c.lower() or c == "ORF"]
    ens_cols = [c for c in gpl.columns if "ensembl" in c.lower()]
    key_cols = [c for c in gpl.columns if "entrez" in c.lower()] + \
               [c for c in gpl.columns if any(k in c.lower() for k in ("gb_acc", "genbank", "refseq"))]

    probe_map, need_api = {}, {}
    for _, row in gpl.iterrows():
        pid = str(row[id_col]).strip()
        hit = None
        for c in ens_cols:
            v = tok(row[c]).split(".")[0].upper()
            if v in lookup:
                hit = lookup[v]; break
        if not hit:
            for c in sym_cols:
                v = tok(row[c]).upper()
                if v in lookup:
                    hit = lookup[v]; break
        if hit:
            probe_map[pid] = hit
        else:
            for c in key_cols:
                k = tok(row[c]).split(".")[0]
                if k:
                    need_api[pid] = k; break

    if need_api:
        import httpx
        client = httpx.Client(timeout=30.0)
        key_to_probes = {}
        for pid, k in need_api.items():
            key_to_probes.setdefault(k, []).append(pid)
        keys = list(key_to_probes)
        for i in range(0, len(keys), 500):
            batch = keys[i:i + 500]
            try:
                r = client.post("https://mygene.info/v3/query", data={
                    "q": ",".join(batch), "scopes": "entrezgene,accession,refseq",
                    "fields": "symbol,ensembl.gene", "species": "human"})
                if r.status_code != 200:
                    continue
                for item in r.json():
                    q, target = str(item.get("query", "")), None
                    e_info = item.get("ensembl")
                    for e in (e_info if isinstance(e_info, list) else [e_info] if e_info else []):
                        eid = str(e.get("gene", "")).split(".")[0].upper()
                        if eid in lookup:
                            target = lookup[eid]; break
                    if not target:
                        sym = str(item.get("symbol", "")).upper()
                        target = lookup.get(sym)
                    if target:
                        for pid in key_to_probes.get(q, []):
                            probe_map.setdefault(pid, target)
            except Exception:
                pass

    pd.Series(probe_map, name="Gene_ID").rename_axis("ID_REF").reset_index().to_csv(PROBE_MAP_CACHE, index=False)
    print(f"  mapped {len(probe_map)} probes; saved {PROBE_MAP_CACHE}")
    return probe_map


if os.path.exists(PROBE_MAP_CACHE):
    pm = pd.read_csv(PROBE_MAP_CACHE, dtype=str)
    probe_map = dict(zip(pm["ID_REF"], pm["Gene_ID"]))
    print(f"\n[3] Loaded cached probe map: {len(probe_map)} probes ({PROBE_MAP_CACHE})")
else:
    probe_map = build_probe_map()

# ---------------------------------------------------------
# 4. ONE PROFILE PER GROUP
# ---------------------------------------------------------
print("\n[4] Writing per-group profiles...")
df["Gene_ID"] = df["ID_REF"].map(probe_map)
mapped = df.dropna(subset=["Gene_ID"])
if mapped.empty:
    raise SystemExit("No probes mapped to Human-GEM genes.")

for grp in ("AD", "Control", "HD"):
    cols = [c for c in groups.get(grp, []) if c in mapped.columns]
    if not cols:
        continue
    vals = mapped[cols].apply(pd.to_numeric, errors="coerce")
    per_probe = vals.mean(axis=1)                                   # mean over samples in the group
    prof = (pd.DataFrame({"Gene_ID": mapped["Gene_ID"], "Mean_Expression": per_probe})
            .dropna().groupby("Gene_ID")["Mean_Expression"].mean().reset_index())   # mean over probes per gene
    prof.to_csv(f"processed_{grp}_expression.csv", index=False)
    print(f"  {grp}: {len(cols)} samples -> {len(prof)} genes -> processed_{grp}_expression.csv")

print("\nDone.")