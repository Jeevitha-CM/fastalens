"""
segmentation.py — sorts sequences into groups (segments) for the segregation charts.

Every function here is plain pandas: no Streamlit. That keeps it testable and reusable.
"""

import pandas as pd

INF = float("inf")

# Each class definition: (bin edges, labels). Bins include the lower edge and exclude the upper edge.
LENGTH_CLASSES = ([0, 200, 500, 1000, 2000, 5000, INF],
                  ["< 200", "200–499", "500–999", "1,000–1,999", "2,000–4,999", "≥ 5,000"])
GC_CLASSES = ([0, 35, 45, 55, 65, INF],
              ["AT-rich (< 35%)", "35–44.9%", "45–54.9%", "55–64.9%", "GC-rich (≥ 65%)"])
ORF_LENGTH_CLASSES = ([0, 100, 300, 500, 1000, INF],
                      ["< 100 aa", "100–299 aa", "300–499 aa", "500–999 aa", "≥ 1,000 aa"])
PI_CLASSES = ([0, 6, 8, INF],
              ["Acidic (pI < 6)", "Neutral (pI 6–8)", "Basic (pI ≥ 8)"])
MW_CLASSES = ([0, 10, 30, 60, 100, INF],
              ["< 10 kDa", "10–29.9 kDa", "30–59.9 kDa", "60–99.9 kDa", "≥ 100 kDa"])
GRAVY_CLASSES = ([-INF, 0, INF],
                 ["Hydrophilic (GRAVY < 0)", "Hydrophobic (GRAVY ≥ 0)"])
STABILITY_CLASSES = ([-INF, 40, INF],
                     ["Stable (index < 40)", "Unstable (index ≥ 40)"])
ISOFORM_CLASSES = ([1, 2, 3, 4, 6, INF],
                   ["1 isoform", "2 isoforms", "3 isoforms", "4–5 isoforms", "≥ 6 isoforms"])

# Chart name -> short explanation (shown under each chart)
CHART_DESCRIPTIONS = {
    "Length classes": "Sequences grouped by length (bp for nucleotides, aa for proteins).",
    "GC content classes": "Nucleotide sequences grouped by GC content.",
    "Molecule type": "Nucleotide versus protein sequences.",
    "Header format": "Source format recognised from the sequence headers.",
    "Source file": "Number of selected sequences contributed by each input file.",
    "ORF status": "Coding potential of nucleotide sequences (longest ORF of at least 50 codons).",
    "ORF length classes": "Length of the longest ORF in nucleotide sequences.",
    "Quality flags": "Sequences without flags (Clean) and with each quality flag. A sequence can carry several flags.",
    "Organism (top 15)": "The 15 most frequent organisms in the headers.",
    "Isoforms per gene (Trinity)": "Trinity transcripts grouped by the number of isoforms of their gene.",
    "Protein charge (pI)": "Proteins grouped by theoretical isoelectric point.",
    "Protein size (MW)": "Proteins grouped by theoretical molecular weight.",
    "Hydropathy (GRAVY)": "Proteins grouped by grand average of hydropathy.",
    "Stability (instability index)": "Proteins grouped by predicted in vitro stability.",
}


def _cut(values, classes):
    edges, labels = classes
    return pd.cut(values, bins=edges, labels=labels, right=False)


def orf_status(table):
    if "ORF complete" not in table.columns:
        return None
    status = pd.Series("Not predicted", index=table.index)
    status[table["Flags"].str.contains("no ORF", regex=False)] = "No ORF (≥ 50 aa)"
    status[table["ORF complete"] == True] = "Complete ORF"      # noqa: E712  (column also holds empty values)
    status[table["ORF complete"] == False] = "Incomplete ORF (no stop codon)"  # noqa: E712
    status[table["Type"] != "DNA/RNA"] = None
    return status


def classify(table):
    """Return one column per chart, each holding the group of every sequence."""
    segments = pd.DataFrame(index=table.index)
    segments["Length classes"] = _cut(table["Length"], LENGTH_CLASSES)
    if "GC%" in table.columns:
        segments["GC content classes"] = _cut(table["GC%"], GC_CLASSES)
    segments["Molecule type"] = table["Type"]
    segments["Header format"] = table["Format"]
    segments["Source file"] = table["Source file"]

    status = orf_status(table)
    if status is not None:
        segments["ORF status"] = status
        segments["ORF length classes"] = _cut(table["ORF length (aa)"], ORF_LENGTH_CLASSES)

    if "Organism" in table.columns:
        segments["Organism (top 15)"] = table["Organism"]

    is_trinity = table["Format"] == "Trinity"
    if is_trinity.any():
        isoform_counts = table[is_trinity].groupby(["Source file", "Gene"])["ID"].transform("count")
        segments["Isoforms per gene (Trinity)"] = _cut(isoform_counts, ISOFORM_CLASSES)

    if "pI" in table.columns:
        segments["Protein charge (pI)"] = _cut(table["pI"], PI_CLASSES)
        segments["Protein size (MW)"] = _cut(table["MW (kDa)"], MW_CLASSES)
        segments["Hydropathy (GRAVY)"] = _cut(table["GRAVY"], GRAVY_CLASSES)
        segments["Stability (instability index)"] = _cut(table["Instability index"], STABILITY_CLASSES)

    return segments


def flag_lists(flags):
    """'short, no ORF' -> ['short', 'no ORF'];  '' -> ['Clean']"""
    return flags.where(flags != "", "Clean").str.split(", ")


def available_charts(segments):
    charts = [name for name in CHART_DESCRIPTIONS if name in segments.columns and segments[name].notna().any()]
    charts.insert(min(len(charts), 7), "Quality flags")
    return charts


def count_groups(segments, flags, chart):
    """Counts and percentages for one chart, as a small table: Group | Sequences | Percent | Label."""
    if chart == "Quality flags":
        counts = flag_lists(flags).explode().value_counts()
        total = len(flags)
    else:
        values = segments[chart]
        if isinstance(values.dtype, pd.CategoricalDtype):
            counts = values.value_counts(sort=False)
        else:
            counts = values.dropna().astype(str).value_counts()
            if chart == "Organism (top 15)":
                counts = counts.head(15)
        total = int(values.notna().sum())

    result = counts.rename_axis("Group").reset_index(name="Sequences")
    result["Group"] = result["Group"].astype(str)
    result["Percent"] = (result["Sequences"] / max(total, 1) * 100).round(1)
    result["Label"] = [f"{n:,} ({p:.1f}%)" for n, p in zip(result["Sequences"], result["Percent"])]
    return result


def selection_mask(segments, flags, chart, groups):
    """True for the sequences that belong to any of the selected groups of one chart."""
    if chart == "Quality flags":
        lists = flag_lists(flags)
        return lists.apply(lambda items: any(group in items for group in groups))
    return segments[chart].astype(str).isin(groups)


# ---------- custom segregation (user chooses the feature, range and number of classes) ----------

NUMERIC_FEATURES = ["Length", "GC%", "N%", "ORF length (aa)", "MW (kDa)", "pI", "GRAVY",
                    "Instability index", "Ambiguous bases"]
CATEGORY_FEATURES = ["Type", "Format", "Source file", "Organism", "Gene", "Database", "Molecule",
                     "ORF strand", "Isoform"]


def custom_features(table):
    numeric = [f for f in NUMERIC_FEATURES if f in table.columns and table[f].notna().any()]
    categories = [f for f in CATEGORY_FEATURES if f in table.columns and table[f].notna().any()]
    return numeric, categories


def custom_classes(values, low, high, n_classes, whole_numbers=False):
    """Split values between low and high into n equal-width classes. Values outside the range are left out."""
    width = (high - low) / n_classes
    edges = [low + width * i for i in range(n_classes + 1)]
    labels = []
    for start, end in zip(edges[:-1], edges[1:]):
        if whole_numbers:
            labels.append(f"{start:,.0f}–{end:,.0f}")
        else:
            labels.append(f"{start:,.2f}–{end:,.2f}")
    edges[-1] = edges[-1] + 1e-9            # so that the maximum value itself is included in the last class
    inside = values.where(values.between(low, high))
    return pd.cut(inside, bins=edges, labels=labels, right=False)


def count_values(values, top=None):
    """Counts and percentages for any Series of groups (ordered classes keep their order)."""
    if isinstance(values.dtype, pd.CategoricalDtype):
        counts = values.value_counts(sort=False)
    else:
        counts = values.dropna().astype(str).value_counts()
        if top:
            counts = counts.head(top)
    total = int(values.notna().sum())
    result = counts.rename_axis("Group").reset_index(name="Sequences")
    result["Group"] = result["Group"].astype(str)
    result["Percent"] = (result["Sequences"] / max(total, 1) * 100).round(1)
    result["Label"] = [f"{n:,} ({p:.1f}%)" for n, p in zip(result["Sequences"], result["Percent"])]
    return result


def is_ordered(segments, chart):
    return chart in segments.columns and isinstance(segments[chart].dtype, pd.CategoricalDtype)