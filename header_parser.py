import re

# ---------- PART A: lookup tables ----------

REFSEQ_TYPES = {
    "NM": "mRNA", "XM": "mRNA (predicted)",
    "NR": "non-coding RNA", "XR": "non-coding RNA (predicted)",
    "NP": "protein", "XP": "protein (predicted)",
    "NC": "genomic (chromosome)", "NG": "genomic (gene region)",
}

UNIPROT_LABELS = {"OS": "Organism", "OX": "Taxon ID", "GN": "Gene", "PE": "Evidence level", "SV": "Version"}


# ---------- PART B: the detective ----------

def detect_format(seq_id):
    if seq_id.startswith("TRINITY_"):
        return "Trinity"
    if seq_id.startswith("sp|") or seq_id.startswith("tr|"):
        return "UniProt"
    if re.match(r"^[A-Z]{2}_\d+\.\d+$", seq_id):
        return "NCBI RefSeq"
    return "Generic"


# ---------- PART C: one specialist reader per database ----------

def read_trinity(seq_id, info):
    parts = seq_id.split("_")
    if len(parts) == 5:
        info["Cluster"] = parts[1]
        info["Gene"] = "_".join(parts[0:4])
        info["Isoform"] = parts[4]


def read_uniprot(seq_id, info):
    parts = seq_id.split("|")
    info["Database"] = "Swiss-Prot (reviewed)" if parts[0] == "sp" else "TrEMBL (unreviewed)"
    info["Accession"] = parts[1]
    info["Entry name"] = parts[2]
    for short, long in UNIPROT_LABELS.items():
        if short in info:
            info[long] = info.pop(short)


def read_ncbi(seq_id, info):
    accession, version = seq_id.split(".")
    info["Accession"] = accession
    info["Version"] = version
    info["Molecule"] = REFSEQ_TYPES.get(accession[:2], "unknown")
    description = info.get("Description", "")

    symbols = re.findall(r"\(([^()]+)\)", description)
    if symbols:
        info["Gene"] = symbols[-1]

    square = re.findall(r"\[([^\[\]]+)\]", description)
    if square:
        info["Organism"] = square[-1]
    else:
        words = description.split()
        if len(words) >= 2:
            info["Organism"] = words[0] + " " + words[1]


# ---------- PART D: the main recipe ----------

def parse_header(header):
    header = header.strip().lstrip(">")

    pieces = header.split(" ", 1)
    seq_id = pieces[0]
    rest = pieces[1] if len(pieces) == 2 else ""
    info = {"ID": seq_id}

    pairs = re.findall(r"(\w+)=(.+?)(?=\s\w+=|$)", rest)
    for key, value in pairs:
        info[key] = value

    if pairs:
        description = rest[: rest.find(pairs[0][0] + "=")].strip()
    else:
        description = rest.strip()
    if description:
        info["Description"] = description

    header_format = detect_format(seq_id)
    info["Format"] = header_format
    if header_format == "Trinity":
        read_trinity(seq_id, info)
    elif header_format == "UniProt":
        read_uniprot(seq_id, info)
    elif header_format == "NCBI RefSeq":
        read_ncbi(seq_id, info)

    return info