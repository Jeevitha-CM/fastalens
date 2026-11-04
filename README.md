# FastaLens

**Parse, characterise, quality-check and export multi-FASTA files in an interactive web application.**

FastaLens reads any FASTA file, separates the information hidden in its headers into columns, detects whether each
sequence is nucleotide or protein, computes the appropriate sequence features, flags quality problems, and exports
clean datasets together with an auditable analysis report.


## Why FastaLens?

Almost every sequence analysis starts with a FASTA file: a de novo transcriptome assembly (e.g. `Trinity.fasta`), a set of reference mRNAs, or a protein list from a proteomics experiment. Before running BLAST, designing primers or
annotating genes, an analyst needs to know what is actually in the file, whether it can be trusted, and how to extract the right subset. This normally requires several command-line tools, web servers and spreadsheets.

FastaLens brings these steps together in one application that can run entirely on the user's own computer.

## Features

| Feature | Description |
|---|---|
| **Multi-format header parsing** | Extracts annotations from Trinity (cluster, gene, isoform), UniProt (accession, organism, gene, evidence level) and NCBI RefSeq (accession, version, molecule type, gene symbol, organism) headers, plus any `key=value` pairs in other formats |
| **Molecule-aware analysis** | Automatic nucleotide/protein detection; GC and N content for nucleotides; molecular weight, isoelectric point, GRAVY and instability index for proteins |
| **Six-frame ORF prediction** | Longest open reading frame per nucleotide sequence on both strands, with 1-based coordinates, completeness and translated protein |
| **Input quality control** | Rejects FASTQ, GenBank and corrupt gzip files with clear messages; reports duplicate IDs, duplicate sequences, IUPAC ambiguity codes, invalid characters, fragments, missing or incomplete ORFs and internal stop codons |
| **Search, filter and sort** | Search by ID, gene, organism or description; filter by molecule type, header format, source file, length and quality; sort by any column |
| **Trinity support** | Gene and isoform counts and extraction of the longest isoform per gene |
| **Multi-file comparison** | Upload several files (plain or `.gz`), then start the analysis with one click and compare them side by side |
| **Example data** | One click loads human mRNAs and haemoglobin proteins to explore the application |
| **Sequence viewer** | Original sequence, reverse complement, translation in any of six frames, longest ORF protein, and residue composition |
| **Overview figures** | Length and GC distributions, and a virtual two-dimensional gel (theoretical pI versus molecular weight) for proteins |
| **Segregation charts** | 14 charts that group sequences into classes (length, GC content, ORF status and length, quality flags, organism, isoforms per gene, protein charge, size, hydropathy and stability), with counts and percentages |
| **Custom segregation** | Any numeric feature split into a user-defined range and number of classes, or any categorical feature by its most frequent values |
| **Click-to-filter** | Clicking a bar restricts the table, statistics, exports and report to that group |
| **Publication-ready charts** | Horizontal bar, vertical bar or pie; academic, colour-blind-safe (Okabe–Ito) or graded palettes; download as PNG, JPEG or SVG at 3× resolution |
| **Export** | Filtered FASTA (original headers), annotation table (CSV), segregation counts (CSV) and an HTML analysis report |
| **Auditable report** | Date and time, software version, SHA-256 checksum of every input file, all analysis settings, QC results, summary statistics and segregation tables |

## Installation

Requires Python 3.10 or later.

```bash
git clone https://github.com/<your-username>/fastalens.git
cd fastalens
python -m venv .venv

# Windows (PowerShell):      .venv\Scripts\Activate.ps1
# Windows (Command Prompt):  .venv\Scripts\activate.bat
# macOS / Linux:             source .venv/bin/activate

pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

The application opens at `http://localhost:8501`. Upload one or more FASTA files and click **Start analysis**, or click **Load example data**. The results are organised in six tabs: Overview, Quality control, Sequence table (with
segregation charts), Sequence viewer, Export and report, and Methods.

### Step-by-step guide

| Step | Where | What to do |
|---|---|---|
| 1 | Start page | Upload one or more FASTA files (`.fasta`, `.fa`, `.fna`, `.faa`, `.fas`, `.txt` or `.gz`) and click **Start analysis**, or click **Load example data** |
| 2 | Sidebar | Switch ORF prediction on or off (off by default for inputs above 20 MB); search; filter by molecule type, header format, file, length and quality; keep the longest isoform per gene for Trinity assemblies; sort by any column |
| 3 | **Overview** | Summary statistics (sequences, total length, N50, mean and median length, GC content, pI, flagged sequences), per-file comparison, length and GC distributions, and the virtual 2D gel for proteins |
| 4 | **Quality control** | Per-file QC table (duplicate IDs and sequences, ambiguous bases, invalid characters, flagged sequences) with warnings |
| 5 | **Sequence table** | All parsed and computed columns, followed by the **segregation charts**: choose charts, chart type and colours; click bars to filter; build a **custom segregation**; download images (camera icon) and counts (CSV) |
| 6 | **Sequence viewer** | Inspect any sequence: original, reverse complement, translation in six frames, longest ORF protein, residue composition and annotation |
| 7 | **Export and report** | Download the HTML analysis report, and prepare filtered FASTA and CSV files on request |
| 8 | **Methods** | Definitions, thresholds, assumptions and references for every calculation |

All results (statistics, charts, exports and report) always refer to the sequences currently selected by the filters and chart selections.

### Example data

| File | Content | Source |
|---|---|---|
| `data/human_mrna.fasta` | Human mRNAs of *HBB*, *HBA1*, *TP53* and *BRCA1* (NM_000518.5, NM_000558.5, NM_000546.6, NM_007294.4) | NCBI RefSeq, downloaded with `scripts/fetch_data.py` |
| `data/hbb_proteins.fasta` | Reviewed haemoglobin beta proteins from multiple species | UniProtKB/Swiss-Prot (CC BY 4.0) |

## Project structure

```
fastalens/
├── app.py                 # Streamlit user interface (display only)
├── header_parser.py       # Header parsing: key=value pairs, Trinity, UniProt, NCBI RefSeq
├── sequence_tools.py      # Molecule detection, GC/N content, protein properties, ORFs, N50, QC flags
├── segmentation.py        # Grouping of sequences into classes for the segregation charts
├── tests/                 # Automated tests (pytest)
├── scripts/fetch_data.py  # Downloads the example mRNAs from NCBI (Entrez)
├── data/                  # Example data
├── learning/              # Practice scripts written while learning
├── docs/                  # Screenshots
├── .streamlit/config.toml # Theme and upload limit (500 MB)
└── requirements.txt
```

The scientific logic (`header_parser.py`, `sequence_tools.py`, `segmentation.py`) is separated from the interface (`app.py`), so every
calculation can be tested without a browser and reused in other projects.

## Testing

```bash
python -m pytest -v
```

26 tests cover header parsing for all supported formats, molecule-type detection, GC content, N50, character
validation, ORF coordinates on both strands, quality flags, segregation class boundaries, and agreement with
reference annotations (below).

## Validation

ORF predictions were compared with the annotated coding sequences of four NCBI RefSeq mRNAs, and protein molecular
weight with UniProt.

| Gene | RefSeq | Expected protein length | FastaLens longest ORF |
|---|---|---|---|
| *HBB* | NM_000518.5 | 147 aa (CDS 51–494) | 147 aa (51–494) |
| *HBA1* | NM_000558.5 | 142 aa | 142 aa |
| *TP53* | NM_000546.6 | 393 aa | 393 aa |
| *BRCA1* | NM_007294.4 | 1,863 aa | 1,863 aa |

The molecular weight computed for β-globin (16.0 kDa) agrees with UniProt P68871 (15,998 Da).

This confirms correct behaviour on well-annotated reference transcripts. It does not establish accuracy on every
input type (see Limitations).

## Methods (summary)

- **Molecule type:** nucleotide if at least 90% of characters (excluding alignment gaps) are A, C, G, T, U or N; otherwise protein.
- **GC content:** (G + C) / length × 100. **N50:** the length at which the cumulative sum of lengths sorted in descending order first reaches 50% of the total.
- **ORFs:** six-frame translation with the standard genetic code; the longest ORF starting with ATG and of at least 50 codons is reported.
- **Protein properties:** Biopython ProtParam (molecular weight, isoelectric point, GRAVY, instability index).
- **Quality flags:** short (nucleotide < 200 bp, protein < 50 aa), high N (> 5%), ambiguous bases, invalid characters, no ORF, incomplete ORF, internal stop, unknown residues; duplicate sequences and duplicate IDs.

Full definitions are available in the **Methods** tab of the application.

## Deployment

FastaLens is deployed on [Streamlit Community Cloud](https://streamlit.io/cloud) directly from this repository (main file `app.py`, dependencies from `requirements.txt`). Every push to the `main` branch updates the web demo
automatically.

## Privacy

When FastaLens is run locally (`streamlit run app.py`), all processing happens on your own computer and sequences are never sent to an external server. On the public web demo, uploaded files are processed on Streamlit Community
Cloud: use it only with public data, never with patient or unpublished sequences.

## Limitations

- Molecule-type detection is heuristic; very short peptides rich in A, C, G and T may be classified as nucleotide.
- The longest ORF is not necessarily the true coding sequence (e.g. for short proteins, non-coding RNAs or transcripts with upstream ORFs). Only ATG start codons and the standard genetic code are used.
- Organism names in NCBI RefSeq headers are inferred heuristically.
- QC thresholds are conventions and may need adjustment for specific applications.
- Files are loaded into memory; inputs larger than about 500 MB may require a command-line workflow. The web demo has
  limited memory and is intended for small files.
- Results are computational predictions for exploratory analysis and require validation before clinical or diagnostic use.

## Roadmap

- [ ] Command-line interface for batch processing in pipelines
- [ ] Motif and primer search on both strands (IUPAC-aware)
- [ ] Codon usage and GC3 analysis
- [ ] GC-outlier detection for possible contamination
- [ ] Figures embedded in the HTML report; PDF export
- [ ] Integration with BLAST annotation

## References

- Cock PJA *et al.* (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25:1422–1423.
- Grabherr MG *et al.* (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nature Biotechnology* 29:644–652.
- Kyte J, Doolittle RF (1982) A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157:105–132.
- Guruprasad K, Reddy BV, Pandit MW (1990) Correlation between stability of a protein and its dipeptide composition. *Protein Engineering* 4:155–161.
- Bjellqvist B *et al.* (1993) The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. *Electrophoresis* 14:1023–1031.
- O'Leary NA *et al.* (2016) Reference sequence (RefSeq) database at NCBI. *Nucleic Acids Research* 44:D733–D745.
- The UniProt Consortium (2023) UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51:D523–D531.
- Streamlit (https://streamlit.io), pandas (https://pandas.pydata.org) and Plotly (https://plotly.com/python/).

## Version history

| Version | Date | Changes |
|---|---|---|
| 1.0.0 | November 2026 | First release: multi-format header parsing, sequence and protein analysis, six-frame ORF prediction, input QC, segregation charts with click-to-filter and custom segregation, sequence viewer, exports and auditable HTML report |

## How to cite

If you use FastaLens in your work, please cite the repository:

> Jeevitha C M (2026). *FastaLens: an interactive application for parsing, characterising, quality-checking and
> exporting multi-FASTA files.* Version 1.0.0. https://github.com/Jeevitha-CM/fastalens

## Contributing and feedback

Bug reports and suggestions are welcome through GitHub Issues. Please include the FastaLens version, a description of the input file (format and size, never confidential data) and the steps to reproduce the problem.

## Author

**Jeevitha C M** · PhD Research Scholar, Genomics and Bioinformatics, SRM Institute of Science and Technology, Chennai · Founder, [Biovagon](https://biovagon.org)

## License

Released under the MIT License. See [LICENSE](LICENSE).
