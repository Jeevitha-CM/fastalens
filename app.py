import gzip
import hashlib
from datetime import datetime
from io import StringIO
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

from header_parser import parse_header
from segmentation import (CHART_DESCRIPTIONS, available_charts, classify, count_groups, count_values,
                          custom_classes, custom_features, is_ordered, selection_mask)
from sequence_tools import analyse_sequence, calculate_n50

VERSION = "1.0.0"
EXAMPLE_FILES = ["data/human_mrna.fasta", "data/hbb_proteins.fasta"]

st.set_page_config(page_title="FastaLens", layout="wide")

st.markdown(
    """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Source+Serif+4:wght@600;700&display=swap');

p, li, label, input, textarea {font-family: 'Inter', sans-serif;}
h2, h3 {font-family: 'Source Serif 4', Georgia, serif; color: #1f4e79;}
.block-container {padding-top: 1.5rem; max-width: 1350px;}

.app-header {background: linear-gradient(135deg, #1f4e79 0%, #2e6ca4 100%); color: #ffffff;
             padding: 1.4rem 1.8rem; border-radius: 10px; margin-bottom: 1.2rem;}
.app-title {font-family: 'Source Serif 4', Georgia, serif; font-size: 2.2rem; font-weight: 700; color: #ffffff;}
.app-version {font-family: 'Inter', sans-serif; font-size: 0.8rem; background: rgba(255, 255, 255, 0.18);
              padding: 2px 8px; border-radius: 10px; vertical-align: middle; margin-left: 6px;}
.app-subtitle {font-size: 1.05rem; opacity: 0.92; margin-top: 0.2rem;}
.app-tags span {display: inline-block; font-size: 0.75rem; border: 1px solid rgba(255, 255, 255, 0.5);
                border-radius: 12px; padding: 2px 10px; margin: 0.7rem 6px 0 0;}

[data-testid="stMetric"] {background: #ffffff; border: 1px solid #dde3ea; border-radius: 8px;
                          padding: 0.8rem 1rem; box-shadow: 0 1px 2px rgba(16, 24, 40, 0.05);}
[data-testid="stMetricLabel"] p {font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.04em; color: #5b6b7c;}
[data-testid="stMetricValue"] {font-size: 1.5rem; color: #1a2b3c;}

.stTabs [data-baseweb="tab-list"] {gap: 4px; border-bottom: 2px solid #dde3ea;}
.stTabs [data-baseweb="tab"] {padding: 0.5rem 1rem; font-weight: 500;}

[data-testid="stSidebar"] {border-right: 1px solid #dde3ea;}
[data-testid="stSidebar"] h2 {font-family: 'Inter', sans-serif; font-size: 0.8rem; text-transform: uppercase;
                              letter-spacing: 0.06em; color: #5b6b7c;}

[data-testid="stFileUploader"] section {border: 1.5px dashed #9fb3c8; border-radius: 8px; background: #ffffff;}
.stButton button, .stDownloadButton button {border-radius: 6px; font-weight: 500;}

.figure-caption {font-size: 0.88rem; color: #5b6b7c; margin-top: -0.4rem; margin-bottom: 1.2rem;}
.app-footer {text-align: center; font-size: 0.8rem; color: #7a8a9a; margin-top: 2rem;}
</style>
""",
    unsafe_allow_html=True,
)

CHART_STYLE = dict(template="simple_white",
                color_discrete_sequence=["#1f4e79", "#c55a11", "#548235", "#7030a0", "#7f6000"])


def caption(text):
    st.markdown(f'<p class="figure-caption">{text}</p>', unsafe_allow_html=True)


REPORT_STYLE = """
<style>
body {font-family: Georgia, serif; max-width: 1000px; margin: 40px auto; color: #1a1a1a;}
h1 {border-bottom: 2px solid #1f4e79; padding-bottom: 6px;}
h2 {color: #1f4e79; margin-top: 28px;}
table {border-collapse: collapse; margin-bottom: 16px;}
th, td {border: 1px solid #cccccc; padding: 4px 8px; text-align: left; font-size: 13px;}
th {background: #f4f6f8;}
</style>
"""


# ---------- chart helpers ----------

PALETTES = {
    "Automatic": None,
    "Academic": ["#1f4e79", "#c55a11", "#548235", "#7030a0", "#bf9000", "#2e75b6", "#843c0c", "#385723",
                 "#7f7f7f", "#5b9bd5"],
    "Colour-blind safe (Okabe–Ito)": ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7",
                                      "#F0E442", "#000000"],
    "Blues (graded)": ["#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"],
    "Greens (graded)": ["#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b"],
}


def pick_colours(n, palette_name, ordered):
    if palette_name == "Automatic":
        palette_name = "Blues (graded)" if ordered else "Academic"
    shades = PALETTES[palette_name]
    if "graded" in palette_name:
        if n == 1:
            return [shades[-2]]
        return [shades[round(i * (len(shades) - 1) / (n - 1))] for i in range(n)]
    return [shades[i % len(shades)] for i in range(n)]


def draw_segregation_chart(counts, title, chart_type, palette_name, ordered, image_format, key=None):
    groups = list(counts["Group"])
    colour_map = dict(zip(groups, pick_colours(len(groups), palette_name, ordered)))
    file_stem = "fastalens_" + "".join(c if c.isalnum() else "_" for c in title.lower()).strip("_")[:60]
    config = {"displaylogo": False,
              "toImageButtonOptions": {"format": image_format.lower(), "filename": file_stem, "scale": 3}}

    if chart_type == "Pie":
        figure_chart = px.pie(counts, names="Group", values="Sequences", color="Group",
                              color_discrete_map=colour_map, hole=0.45, template="simple_white")
        figure_chart.update_traces(sort=False, direction="clockwise", textinfo="percent",
                                   hovertemplate="%{label}<br>%{value:,} sequences (%{percent})<extra></extra>")
        figure_chart.update_layout(title=title, height=380, margin=dict(l=10, r=10, t=50, b=10))
        st.plotly_chart(figure_chart, config=config, key=key)
        return

    if chart_type == "Horizontal bar":
        figure_chart = px.bar(counts, x="Sequences", y="Group", orientation="h", text="Label", color="Group",
                              color_discrete_map=colour_map, custom_data=["Group"], template="simple_white")
        figure_chart.update_layout(yaxis=dict(categoryorder="array", categoryarray=groups[::-1]),
                                   xaxis_title="Number of sequences", yaxis_title=None,
                                   height=110 + 34 * len(groups))
    else:
        figure_chart = px.bar(counts, x="Group", y="Sequences", text="Label", color="Group",
                              color_discrete_map=colour_map, custom_data=["Group"], template="simple_white")
        figure_chart.update_layout(xaxis=dict(categoryorder="array", categoryarray=groups),
                                   yaxis_title="Number of sequences", xaxis_title=None, height=420)
    figure_chart.update_traces(textposition="outside", cliponaxis=False,
                               hovertemplate="%{customdata[0]}<br>%{text}<extra></extra>")
    figure_chart.update_layout(title=title, showlegend=False, margin=dict(l=10, r=40, t=50, b=40))
    if key:
        st.plotly_chart(figure_chart, config=config, on_select="rerun", selection_mode="points", key=key)
    else:
        st.plotly_chart(figure_chart, config=config)


# ---------- PART 1: read each file into a table ----------

def load_fasta(file_name, file_bytes, find_orfs):
    if file_name.endswith(".gz"):
        try:
            file_bytes = gzip.decompress(file_bytes)
        except (OSError, EOFError):
            raise ValueError("the name ends in .gz but the file is not a valid gzip archive.")
    text = file_bytes.decode("utf-8", errors="replace")

    beginning = text.lstrip()[:10]
    if beginning.startswith("@"):
        raise ValueError("this looks like a FASTQ file (the first line starts with '@'). Convert it to FASTA first.")
    if beginning.startswith("LOCUS"):
        raise ValueError("this looks like a GenBank file. Export it in FASTA format first.")
    if not beginning.startswith(">"):
        raise ValueError("this is not a FASTA file (the first line must start with '>').")

    rows = []
    for record in SeqIO.parse(StringIO(text), "fasta"):
        sequence = str(record.seq)
        if len(sequence) == 0:
            continue
        row = {"Source file": file_name}
        row.update(parse_header(record.description))
        row.update(analyse_sequence(sequence, find_orf=find_orfs))
        row["Header"] = record.description
        row["Sequence"] = sequence
        rows.append(row)
    return pd.DataFrame(rows)


def checksum(file_bytes):
    return hashlib.sha256(file_bytes).hexdigest()


class ExampleFile:
    def __init__(self, path):
        self.name = Path(path).name
        self.data = Path(path).read_bytes()
        self.size = len(self.data)

    def getvalue(self):
        return self.data


def qc_summary(group):
    return {
        "Sequences": len(group),
        "Duplicate IDs": int(group["Duplicate ID"].sum()),
        "Duplicate sequences": int(group["Duplicate"].sum()),
        "With ambiguous bases": int((group["Ambiguous bases"] > 0).sum()),
        "With invalid characters": int((group["Invalid characters"] > 0).sum()),
        "Flagged": int((group["Flags"] != "").sum()),
    }


st.markdown(
    f"""
<div class="app-header">
<div class="app-title">FastaLens <span class="app-version">v{VERSION}</span></div>
<div class="app-subtitle">Parse, characterise, quality-check and export multi-FASTA files</div>
<div class="app-tags"><span>Trinity</span><span>UniProt</span><span>NCBI RefSeq</span><span>Nucleotide and protein</span><span>Runs locally</span></div>
</div>
""",
    unsafe_allow_html=True,
)

uploaded_files = st.file_uploader("Input files (FASTA, optionally gzip-compressed)",
                                type=["fasta", "fa", "fna", "faa", "fas", "txt", "gz"],
                                accept_multiple_files=True)

if not uploaded_files:
    if st.session_state.get("use_example"):
        uploaded_files = [ExampleFile(path) for path in EXAMPLE_FILES if Path(path).exists()]
        st.session_state["analysed_files"] = [file.name for file in uploaded_files]
        if st.button("Clear example data"):
            st.session_state["use_example"] = False
            st.rerun()
    elif st.button("Load example data"):
        st.session_state["use_example"] = True
        st.rerun()

if not uploaded_files:
    st.info("Upload one or more FASTA files to begin the analysis, or click Load example data.")
    st.stop()

current_files = [file.name for file in uploaded_files]
if st.session_state.get("analysed_files") != current_files:
    st.info(f"{len(current_files)} file(s) ready: {', '.join(current_files)}. "
            "Add more files if needed, then click Start analysis.")
    if st.button("Start analysis", type="primary"):
        st.session_state["analysed_files"] = current_files
        st.rerun()
    st.stop()

total_mb = sum(file.size for file in uploaded_files) / 1_000_000
find_orfs = st.sidebar.checkbox("Predict ORFs (slower for large files)", value=total_mb < 20,
                                help="Six-frame ORF search. Switched off by default when the input is larger than 20 MB.")

# Read the files only when the file list or the ORF setting changes, and keep the result in memory.
table_key = (tuple((file.name, file.size) for file in uploaded_files), find_orfs)
if st.session_state.get("table_key") != table_key:
    tables = []
    errors = []
    file_rows = []
    with st.spinner("Reading sequences... large files can take a few minutes."):
        for file in uploaded_files:
            data = file.getvalue()
            try:
                tables.append(load_fasta(file.name, data, find_orfs))
                file_rows.append({"File": file.name, "Size (bytes)": f"{len(data):,}",
                                  "SHA-256 checksum": checksum(data)})
            except ValueError as error:
                errors.append(f"{file.name}: {error}")

        if tables:
            table = pd.concat(tables, ignore_index=True)
        else:
            table = pd.DataFrame()
        if not table.empty:
            table["Duplicate"] = table["Sequence"].str.upper().duplicated(keep=False)
            table["Duplicate ID"] = table.duplicated(subset=["Source file", "ID"], keep=False)

    st.session_state["table"] = table
    st.session_state["errors"] = errors
    st.session_state["file_info"] = pd.DataFrame(file_rows)
    st.session_state["table_key"] = table_key

table = st.session_state["table"]
for message in st.session_state["errors"]:
    st.error(message)

if table.empty:
    st.error("No sequences were found. Please check that the files are in FASTA format.")
    st.stop()

multiple_files = table["Source file"].nunique() > 1


# ---------- PART 2: sidebar controls ----------

st.sidebar.header("Search and filter")

search = st.sidebar.text_input("Search", placeholder="ID, gene, organism, description")

all_types = sorted(table["Type"].unique())
chosen_types = st.sidebar.multiselect("Molecule type", all_types, default=all_types)

all_formats = sorted(table["Format"].unique())
chosen_formats = st.sidebar.multiselect("Header format", all_formats, default=all_formats)

all_files = sorted(table["Source file"].unique())
chosen_files = st.sidebar.multiselect("Source file", all_files, default=all_files)

shortest = int(table["Length"].min())
longest = int(table["Length"].max())
if shortest < longest:
    length_range = st.sidebar.slider("Length range", shortest, longest, (shortest, longest))
else:
    length_range = (shortest, longest)

unigenes_only = False
if (table["Format"] == "Trinity").any():
    unigenes_only = st.sidebar.checkbox("Longest isoform per gene (Trinity)")

exclude_flagged = st.sidebar.checkbox("Exclude flagged sequences")
remove_duplicates = st.sidebar.checkbox("Remove duplicate sequences (keep one copy)")

st.sidebar.header("Sort")
sort_options = [column for column in table.columns if column not in ("Header", "Sequence")]
sort_column = st.sidebar.selectbox("Sort by", sort_options, index=sort_options.index("Length"))
largest_first = st.sidebar.toggle("Descending order", value=True)


# ---------- PART 3: apply the sieves ----------

view = table[
    table["Type"].isin(chosen_types)
    & table["Format"].isin(chosen_formats)
    & table["Source file"].isin(chosen_files)
    & table["Length"].between(length_range[0], length_range[1])
]

if search:
    found = pd.Series(False, index=view.index)
    for column in ["ID", "Description", "Gene", "Organism"]:
        if column in view.columns:
            matches = view[column].fillna("").astype(str).str.contains(search, case=False, regex=False)
            found = found | matches
    view = view[found]

if unigenes_only:
    view = view.sort_values("Length", ascending=False).drop_duplicates(subset="Gene")

if exclude_flagged:
    view = view[view["Flags"] == ""]

if remove_duplicates:
    view = view.drop_duplicates(subset="Sequence")

view = view.sort_values(sort_column, ascending=not largest_first)

# ---------- chart selections (clicking a bar in the segregation charts filters the data) ----------

base_view = view
segments = classify(base_view)
chart_list = available_charts(segments)
reset_number = st.session_state.get("chart_reset", 0)

chart_filters = {}
for chart in chart_list:
    state = st.session_state.get(f"chart_{chart}_{reset_number}")
    if state and state["selection"]["points"]:
        groups = []
        for point in state["selection"]["points"]:
            group = point.get("customdata", [point.get("y")])[0]
            if group not in groups:
                groups.append(group)
        chart_filters[chart] = groups

for chart, groups in chart_filters.items():
    view = view[selection_mask(segments, base_view["Flags"], chart, groups).loc[view.index]]

if view.empty:
    st.warning("No sequences match the current filters.")
    st.stop()

dna = view[view["Type"] == "DNA/RNA"]
protein = view[view["Type"] == "Protein"]


# ---------- STATUS LINE AND TABS ----------

st.info(f"Files analysed: **{table['Source file'].nunique()}**  ·  "
        f"Sequences loaded: **{len(table):,}**  ·  Selected by filters: **{len(view):,}**")

tab_overview, tab_qc, tab_table, tab_viewer, tab_export, tab_methods = st.tabs(
    ["Overview", "Quality control", "Sequence table", "Sequence viewer", "Export and report", "Methods"]
)


# ---------- SECTION 1: overview ----------

with tab_overview:
    st.subheader("Dataset overview")

    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Sequences", f"{len(view):,}")
    col2.metric("Total length (residues)", f"{int(view['Length'].sum()):,}")
    col3.metric("N50", f"{calculate_n50(view['Length']):,}")
    col4.metric("Mean length", f"{view['Length'].mean():,.1f}")
    col5.metric("Median length", f"{int(view['Length'].median()):,}")

    col6, col7, col8, col9 = st.columns(4)
    col6.metric("Nucleotide sequences", f"{len(dna):,}")
    col7.metric("Protein sequences", f"{len(protein):,}")
    col8.metric("Mean GC content (%)", f"{dna['GC%'].mean():.2f}" if not dna.empty else "n/a")
    col9.metric("Mean pI", f"{protein['pI'].mean():.2f}" if not protein.empty else "n/a")

    if multiple_files:
        st.subheader("Per-file comparison")
        groups = view.groupby("Source file")
        per_file = groups.agg(Sequences=("ID", "count"),
                            Total=("Length", "sum"),
                            Mean=("Length", "mean"),
                            Median=("Length", "median"))
        per_file.columns = ["Sequences", "Total length", "Mean length", "Median length"]
        per_file["N50"] = groups["Length"].apply(calculate_n50)
        if "GC%" in view.columns:
            per_file["Mean GC (%)"] = groups["GC%"].mean()
        st.dataframe(per_file.round(2))

    st.subheader("Length and composition")
    figure = 1

    left, right = st.columns(2)
    with left:
        colour_by = "Source file" if multiple_files else None
        chart = px.histogram(view, x="Length", nbins=60, color=colour_by,
                            barmode="overlay", opacity=0.75, **CHART_STYLE)
        chart.update_layout(xaxis_title="Length (bp or aa)", yaxis_title="Number of sequences")
        st.plotly_chart(chart)
        caption(f"Figure {figure}. Length distribution of the selected sequences.")
        figure = figure + 1

    with right:
        if not dna.empty:
            chart = px.histogram(dna, x="GC%", nbins=50, **CHART_STYLE)
            chart.update_layout(xaxis_title="GC content (%)", yaxis_title="Number of sequences")
            st.plotly_chart(chart)
            caption(f"Figure {figure}. GC content distribution of nucleotide sequences.")
            figure = figure + 1

    if not protein.empty:
        chart = px.scatter(protein, x="pI", y="MW (kDa)", hover_name="ID", log_y=True, **CHART_STYLE)
        chart.update_layout(xaxis_title="Isoelectric point (pI)", yaxis_title="Molecular weight (kDa, log scale)")
        st.plotly_chart(chart)
        caption(f"Figure {figure}. Virtual two-dimensional gel: theoretical isoelectric point versus molecular weight.")


# ---------- SECTION 2 + 3: quality control and table ----------

with tab_qc:
    st.subheader("Input quality control")

    qc_table = pd.DataFrame([{"File": name, **qc_summary(group)} for name, group in table.groupby("Source file")])
    st.dataframe(qc_table, hide_index=True)

    if qc_table["Duplicate IDs"].sum() > 0:
        st.warning("Some sequence IDs occur more than once. Many tools (for example makeblastdb) "
                "require unique IDs. Sort the table by 'Duplicate ID' to inspect them.")
    if qc_table["With invalid characters"].sum() > 0:
        st.warning("Some sequences contain characters that are not valid nucleotide or amino-acid codes "
                "(for example alignment gaps, digits or punctuation). ORF prediction was skipped for these.")
    caption("Quality control is computed on all input sequences, before filtering.")


with tab_table:
    st.subheader("Sequence table")

    if chart_filters:
        described = "; ".join(f"{chart}: {', '.join(groups)}" for chart, groups in chart_filters.items())
        info_column, button_column = st.columns([4, 1])
        info_column.info(f"Chart selection active. {described}")
        if button_column.button("Clear chart selection"):
            st.session_state["chart_reset"] = reset_number + 1
            st.rerun()

    st.dataframe(view.drop(columns=["Header", "Sequence"]).dropna(axis=1, how="all"), hide_index=True, height=420)
    caption(f"{len(view):,} sequences shown.")

    # ---------- segregation charts ----------
    st.divider()
    st.subheader("Segregation charts")
    st.markdown("Choose the charts to display. In bar charts, **click a bar** to show only those sequences in the "
                "table above (hold **Shift** to select several bars). To save a chart, use the **camera icon** that "
                "appears at its top right when you point at it.")

    default_charts = [chart for chart in ["Length classes", "GC content classes", "ORF status", "Quality flags",
                                          "Protein charge (pI)", "Protein size (MW)"] if chart in chart_list]
    chosen_charts = st.multiselect("Charts", chart_list, default=default_charts)

    style_1, style_2, style_3 = st.columns(3)
    chart_type = style_1.radio("Chart type", ["Horizontal bar", "Vertical bar", "Pie"], horizontal=True)
    palette_name = style_2.selectbox("Colours", list(PALETTES))
    image_format = style_3.radio("Image format", ["PNG", "JPEG", "SVG"], horizontal=True)

    segregation_tables = {}
    chart_columns = st.columns(2)
    for number, chart in enumerate(chosen_charts):
        counts = count_groups(segments, base_view["Flags"], chart)
        segregation_tables[chart] = counts
        with chart_columns[number % 2].container(border=True):
            draw_segregation_chart(counts, chart, chart_type, palette_name, is_ordered(segments, chart),
                                   image_format, key=f"chart_{chart}_{reset_number}")
            caption(CHART_DESCRIPTIONS[chart])

    # ---------- custom segregation ----------
    st.markdown("#### Custom segregation")
    st.markdown("Build your own chart: choose a feature, and for numeric features the range and number of classes.")

    numeric_features, category_features = custom_features(base_view)
    custom_1, custom_2, custom_3 = st.columns([2, 3, 1])
    feature = custom_1.selectbox("Feature", numeric_features + category_features)

    if feature in numeric_features:
        values = base_view[feature]
        whole_numbers = feature in ("Length", "ORF length (aa)", "Ambiguous bases")
        lowest = float(values.min())
        highest = float(values.max())
        if lowest < highest:
            low, high = custom_2.slider("Range", lowest, highest, (lowest, highest))
        else:
            low, high = lowest, highest + 1
        n_classes = custom_3.number_input("Classes", min_value=2, max_value=30, value=5)
        custom_groups = custom_classes(values, low, high, int(n_classes), whole_numbers)
        custom_title = f"{feature}: {n_classes} classes from {low:,.2f} to {high:,.2f}"
        ordered = True
    else:
        top_n = custom_3.number_input("Top groups", min_value=2, max_value=50, value=15)
        custom_groups = base_view[feature]
        custom_title = f"{feature} (top {top_n})"
        ordered = False

    custom_counts = count_values(custom_groups, top=int(top_n) if not ordered else None)
    if custom_counts["Sequences"].sum() > 0:
        with st.container(border=True):
            draw_segregation_chart(custom_counts, custom_title, chart_type, palette_name, ordered, image_format)
            caption("Custom chart. Lower bounds are included and upper bounds excluded, except that the last class "
                    "also includes the maximum. Sequences outside the range are not counted.")
        segregation_tables[custom_title] = custom_counts
    else:
        st.info("No sequences fall in this range.")

    if segregation_tables:
        all_counts = pd.concat([counts.assign(Chart=name)[["Chart", "Group", "Sequences", "Percent"]]
                                for name, counts in segregation_tables.items()], ignore_index=True)
        st.download_button("Download segregation counts (CSV)", all_counts.to_csv(index=False),
                           file_name="fastalens_segregation_counts.csv", mime="text/csv")


# ---------- SECTION 4: sequence viewer ----------

with tab_viewer:
    st.subheader("Sequence viewer")

    candidate_ids = view["ID"].head(1000).tolist()
    if len(view) > 1000:
        st.caption("The list shows the first 1,000 sequences of the current view. Use the search box to narrow it.")
    selected_id = st.selectbox("Select a sequence", candidate_ids)
    selected = view[view["ID"] == selected_id].iloc[0]
    sequence = selected["Sequence"].upper()

    st.markdown(f"**Header:** `{selected['Header']}`")

    options = ["Original sequence"]
    if selected["Type"] == "DNA/RNA":
        options = options + ["Reverse complement", "Translation"]
        if pd.notna(selected.get("ORF protein")):
            options.append("Longest ORF protein")
    display_mode = st.radio("Display", options, horizontal=True)

    viewer_dna = sequence.replace("U", "T")
    if display_mode == "Reverse complement":
        shown = str(Seq(viewer_dna).reverse_complement())
    elif display_mode == "Translation":
        frame_choice = st.radio("Reading frame", ["+1", "+2", "+3", "-1", "-2", "-3"], horizontal=True)
        offset = int(frame_choice[1]) - 1
        if frame_choice.startswith("-"):
            frame_dna = str(Seq(viewer_dna).reverse_complement())
        else:
            frame_dna = viewer_dna
        frame_dna = frame_dna[offset:]
        usable = len(frame_dna) // 3 * 3
        try:
            shown = str(Seq(frame_dna[:usable]).translate())
        except TranslationError:
            st.warning("This sequence contains characters that cannot be translated. Showing the original sequence.")
            shown = sequence
    elif display_mode == "Longest ORF protein":
        shown = selected["ORF protein"]
    else:
        shown = sequence

    viewer_left, viewer_right = st.columns([3, 2])

    with viewer_left:
        formatted = "\n".join(shown[i: i + 60] for i in range(0, len(shown), 60))
        st.code(formatted, language=None)
        annotation = selected.drop(labels=["Header", "Sequence"]).dropna().astype(str)
        st.dataframe(annotation.rename("Value"))

    with viewer_right:
        counts = pd.Series(list(shown)).value_counts().sort_index()
        composition = pd.DataFrame({"Residue": counts.index,
                                    "Frequency (%)": (counts.values / len(shown) * 100).round(2)})
        chart = px.bar(composition, x="Residue", y="Frequency (%)", **CHART_STYLE)
        st.plotly_chart(chart)
        caption(f"Residue composition of the displayed sequence ({display_mode.lower()}, {len(shown):,} characters).")


# ---------- SECTION 5: export and report ----------

with tab_export:
    st.subheader("Export and report")
    st.markdown("Download the filtered sequences and their annotations for downstream analysis "
                "(for example BLAST, primer design or functional annotation).")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    file_info = st.session_state["file_info"]

    settings = pd.DataFrame({
        "Setting": ["Search text", "Molecule types", "Header formats", "Source files", "Length range",
                    "ORF prediction", "Longest isoform per gene", "Exclude flagged", "Remove duplicate sequences",
                    "Sort", "Chart selection"],
        "Value": [search or "(none)", ", ".join(chosen_types), ", ".join(chosen_formats), ", ".join(chosen_files),
                f"{length_range[0]:,} to {length_range[1]:,}", find_orfs, unigenes_only, exclude_flagged,
                remove_duplicates, f"{sort_column} ({'descending' if largest_first else 'ascending'})",
                "; ".join(f"{chart}: {', '.join(groups)}" for chart, groups in chart_filters.items()) or "(none)"],
    })

    summary = pd.DataFrame({
        "Statistic": ["Sequences selected", "Total length", "N50", "Mean length", "Median length",
                    "Nucleotide sequences", "Protein sequences"],
        "Value": [f"{len(view):,}", f"{int(view['Length'].sum()):,}", f"{calculate_n50(view['Length']):,}",
                f"{view['Length'].mean():,.1f}", f"{int(view['Length'].median()):,}",
                f"{len(dna):,}", f"{len(protein):,}"],
    })

    segregation_html = "".join(f"<h3>{name}</h3>" + counts[["Group", "Sequences", "Percent"]].to_html(index=False)
                               for name, counts in segregation_tables.items())
    if not segregation_html:
        segregation_html = "<p>No segregation charts were selected.</p>"

    report = f"""<!DOCTYPE html>
    <html><head><meta charset="utf-8"><title>FastaLens report</title>{REPORT_STYLE}</head><body>
    <h1>FastaLens analysis report</h1>
    <p>Generated: {timestamp}<br>Software: FastaLens version {VERSION}</p>
    <h2>1. Input files</h2>
    {file_info.to_html(index=False)}
    <p>Checksums are computed on the uploaded bytes (the compressed file for .gz input).</p>
    <h2>2. Analysis settings</h2>
    {settings.to_html(index=False)}
    <h2>3. Input quality control</h2>
    {qc_table.to_html(index=False)}
    <h2>4. Summary of selected sequences</h2>
    {summary.to_html(index=False)}
    <h2>5. Segregation of selected sequences</h2>
    <p>Counts and percentages for the charts selected in the Sequence table tab (Percent = share of sequences
    with a value for that chart).</p>
    {segregation_html}
    <h2>6. Note</h2>
    <p>Results are computational predictions for exploratory analysis and require experimental or
    database validation before clinical or diagnostic interpretation.</p>
    </body></html>"""

    st.download_button("Download analysis report (HTML)", report,
                       file_name=f"fasta_report_{datetime.now():%Y%m%d_%H%M}.html", mime="text/html")

    st.markdown("**Sequence files**")
    prepare = st.toggle("Prepare FASTA and CSV downloads",
                        help="The files are built only when requested, which keeps the app fast for large datasets.")
    if prepare:
        fasta_text = "".join(f">{header}\n{sequence}\n" for header, sequence in zip(view["Header"], view["Sequence"]))

        include_sequences = st.checkbox("Include sequences in the CSV table", value=True)
        if include_sequences:
            columns_to_drop = ["Header"]
        else:
            columns_to_drop = ["Header", "Sequence"]
        csv_text = view.drop(columns=columns_to_drop).dropna(axis=1, how="all").to_csv(index=False)

        if include_sequences and view["Length"].max() > 32767:
            st.warning("Some sequences are longer than 32,767 characters, the maximum Excel can display in one cell. "
                       "Use the FASTA download for complete sequences.")

        export_1, export_2 = st.columns(2)
        export_1.download_button("Download sequences (FASTA)", fasta_text,
                                 file_name="filtered_sequences.fasta", mime="text/plain")
        export_2.download_button("Download annotation table (CSV)", csv_text,
                                 file_name="sequence_annotations.csv", mime="text/csv")


# ---------- SECTION 6: methods ----------

with tab_methods:
    st.subheader("Methods")
    with st.expander("Definitions, assumptions and references"):
        st.markdown("""
**Molecule type.** A sequence is classified as nucleotide when at least 90% of its characters are
A, C, G, T, U or N, and as protein otherwise. Short peptides composed mainly of A, C, G and T may be misclassified.

**GC content.** (G + C) / sequence length × 100.

**N50.** Lengths are sorted in descending order; N50 is the length at which the cumulative sum
first reaches 50% of the total length.

**Open reading frames.** All six reading frames are translated with the standard genetic code. The longest
ORF beginning with ATG and of at least 50 codons is reported, with 1-based coordinates on the input sequence.
ORFs without a stop codon are reported as incomplete. The longest ORF is not necessarily the true coding sequence.
In the sequence viewer, frames +1 to +3 start at positions 1 to 3 of the input sequence, and frames -1 to -3 start
at positions 1 to 3 of its reverse complement.

**Protein properties.** Molecular weight, isoelectric point, GRAVY and instability index were computed
with Biopython ProtParam after removal of non-standard residues (X, B, Z, *).

**Input validation.** Files are rejected when they are not FASTA (for example FASTQ or GenBank) or when a .gz file
is not a valid gzip archive. IUPAC ambiguity codes (R, Y, K, M, S, W, B, D, H, V) are counted as ambiguous bases;
any other non-standard character (for example alignment gaps, digits or punctuation) is counted as invalid, and ORF
prediction is skipped for such sequences.

**Quality flags.** short (nucleotide < 200 bp; protein < 50 aa), high N (N > 5%), ambiguous bases, invalid
characters, no ORF, incomplete ORF, internal stop, unknown residues. Duplicate marks sequences that occur more than
once across all input files; Duplicate ID marks repeated identifiers.
These thresholds are conventions and should be adapted to the application.

**Header parsing.** Key–value pairs (key=value) are extracted from all headers. Format-specific fields are
extracted for Trinity (cluster, gene, isoform), UniProt (accession, entry name, organism, gene) and
NCBI RefSeq (accession, version, molecule type, gene symbol, organism). The organism in NCBI headers is
inferred heuristically and should be verified.

**Segregation charts.** Sequences are grouped into classes (lower bound included, upper bound excluded):
length < 200, 200–499, 500–999, 1,000–1,999, 2,000–4,999 and ≥ 5,000; GC content < 35% (AT-rich), 35–45%, 45–55%,
55–65% and ≥ 65% (GC-rich); ORF length < 100, 100–299, 300–499, 500–999 and ≥ 1,000 aa; pI < 6 (acidic), 6–8
(neutral) and ≥ 8 (basic); molecular weight < 10, 10–30, 30–60, 60–100 and ≥ 100 kDa; GRAVY < 0 (hydrophilic) or
≥ 0 (hydrophobic); instability index < 40 (stable) or ≥ 40 (unstable; Guruprasad et al., 1990). Selecting bars in
a chart restricts all results to the selected groups. Custom segregation divides a user-defined range of a numeric feature into
equal-width classes, or counts the most frequent values of a categorical feature.

**Longest isoform per gene.** For each Trinity gene identifier, only the longest transcript is retained.

**Analysis report.** The HTML report records the date and time, software version, SHA-256 checksum of each input
file, all filter settings, quality-control results and summary statistics, to support reproducibility and audit.

**References**
- Cock PJA *et al.* (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25:1422–1423.
- Grabherr MG *et al.* (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nature Biotechnology* 29:644–652.
- Kyte J, Doolittle RF (1982) A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157:105–132.
- Guruprasad K, Reddy BV, Pandit MW (1990) Correlation between stability of a protein and its dipeptide composition. *Protein Engineering* 4:155–161.
- Bjellqvist B *et al.* (1993) The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. *Electrophoresis* 14:1023–1031.
""")

st.markdown(f'<div class="app-footer">FastaLens v{VERSION} · Jeevitha C M ·  · All processing runs locally; sequences are not sent to any external server.</div>', unsafe_allow_html=True)