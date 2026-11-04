from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

DNA_LETTERS = set("ACGTUN")
PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWY")
IUPAC_AMBIGUOUS = set("RYKMSWBDHV")
VALID_NUCLEOTIDE = DNA_LETTERS | IUPAC_AMBIGUOUS
VALID_PROTEIN = PROTEIN_LETTERS | set("BXZJUO*")


def detect_type(sequence):
    letters = sequence.upper().replace("-", "").replace(".", "")
    if len(letters) == 0:
        return "Empty"
    dna_count = sum(letters.count(letter) for letter in "ACGTUN")
    if dna_count / len(letters) >= 0.9:
        return "DNA/RNA"
    return "Protein"


def character_check(sequence, sequence_type):
    sequence = sequence.upper()
    if sequence_type == "DNA/RNA":
        ambiguous = sum(sequence.count(letter) for letter in IUPAC_AMBIGUOUS)
        valid = sum(sequence.count(letter) for letter in VALID_NUCLEOTIDE)
    else:
        ambiguous = 0
        valid = sum(sequence.count(letter) for letter in VALID_PROTEIN)
    invalid = len(sequence) - valid
    return ambiguous, invalid


def dna_features(sequence):
    sequence = sequence.upper()
    length = len(sequence)
    gc = (sequence.count("G") + sequence.count("C")) / length * 100
    n = sequence.count("N") / length * 100
    return {"GC%": round(gc, 2), "N%": round(n, 2)}


def protein_features(sequence):
    clean = "".join(letter for letter in sequence.upper() if letter in PROTEIN_LETTERS)
    if len(clean) == 0:
        return {}
    analysis = ProteinAnalysis(clean)
    return {
        "MW (kDa)": round(analysis.molecular_weight() / 1000, 2),
        "pI": round(analysis.isoelectric_point(), 2),
        "GRAVY": round(analysis.gravy(), 3),
        "Instability index": round(analysis.instability_index(), 1),
    }


def analyse_sequence(sequence, find_orf=True):
    sequence_type = detect_type(sequence)
    result = {"Type": sequence_type, "Length": len(sequence)}
    flags = []

    ambiguous, invalid = character_check(sequence, sequence_type)
    result["Ambiguous bases"] = ambiguous
    result["Invalid characters"] = invalid
    if ambiguous > 0:
        flags.append("ambiguous bases")
    if invalid > 0:
        flags.append("invalid characters")

    if sequence_type == "DNA/RNA":
        result.update(dna_features(sequence))
        if len(sequence) < 200:
            flags.append("short")
        if result["N%"] > 5:
            flags.append("high N")
        if find_orf and invalid == 0:
            orf = longest_orf(sequence)
            if orf is None:
                flags.append("no ORF")
            else:
                result.update(orf)
                if not orf["ORF complete"]:
                    flags.append("incomplete ORF")

    elif sequence_type == "Protein":
        result.update(protein_features(sequence))
        if len(sequence) < 50:
            flags.append("short")
        if "*" in sequence.rstrip("*"):
            flags.append("internal stop")
        if "X" in sequence.upper():
            flags.append("unknown residues")

    result["Flags"] = ", ".join(flags)
    return result


def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    half_total = sum(sorted_lengths) / 2
    running_total = 0
    for length in sorted_lengths:
        running_total = running_total + length
        if running_total >= half_total:
            return int(length)
    return 0


def longest_orf(sequence, min_protein_length=50):
    dna = sequence.upper().replace("U", "T")
    length = len(dna)
    strands = [("+", dna), ("-", str(Seq(dna).reverse_complement()))]
    best = None

    for strand, strand_dna in strands:
        for frame in range(3):
            usable = (length - frame) // 3 * 3
            protein = str(Seq(strand_dna[frame: frame + usable]).translate())

            aa_position = 0
            pieces = protein.split("*")
            for number, piece in enumerate(pieces):
                start_in_piece = piece.find("M")
                if start_in_piece != -1:
                    orf_protein = piece[start_in_piece:]
                    has_stop = number < len(pieces) - 1
                    if len(orf_protein) >= min_protein_length:
                        if best is None or len(orf_protein) > best["ORF length (aa)"]:
                            nt_start = frame + (aa_position + start_in_piece) * 3
                            nt_length = (len(orf_protein) + (1 if has_stop else 0)) * 3
                            if strand == "+":
                                start, end = nt_start + 1, nt_start + nt_length
                            else:
                                start, end = length - (nt_start + nt_length) + 1, length - nt_start
                            best = {
                                "ORF strand": strand,
                                "ORF start": start,
                                "ORF end": end,
                                "ORF length (aa)": len(orf_protein),
                                "ORF complete": has_stop,
                                "ORF protein": orf_protein,
                            }
                aa_position = aa_position + len(piece) + 1
    return best