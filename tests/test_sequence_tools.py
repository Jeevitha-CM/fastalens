from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from sequence_tools import (analyse_sequence, calculate_n50, character_check,
                            detect_type, dna_features, longest_orf, protein_features)

EXAMPLE = Path("data/human_mrna.fasta")


# ---------- molecule type ----------

def test_detect_dna():
    assert detect_type("ATGCGTACGTTAGC") == "DNA/RNA"


def test_detect_protein():
    assert detect_type("MVHLTPEEKSAVTALWGKV") == "Protein"


def test_detect_aligned_dna_with_gaps():
    assert detect_type("ATG---CGTACG--TTAGC") == "DNA/RNA"


# ---------- simple calculations ----------

def test_gc_content():
    assert dna_features("GGCCAATT")["GC%"] == 50.0


def test_gc_content_lowercase():
    assert dna_features("ggcc")["GC%"] == 100.0


def test_n50_worked_example():
    # total 1000, half 500: 400 + 300 = 700 >= 500, so N50 = 300
    assert calculate_n50([100, 200, 300, 400]) == 300


def test_character_check():
    assert character_check("ACGTRY", "DNA/RNA") == (2, 0)     # R and Y are ambiguous
    assert character_check("ACG-T1", "DNA/RNA") == (0, 2)     # '-' and '1' are invalid


# ---------- ORF finder ----------

def test_orf_on_plus_strand():
    sequence = "CCCCC" + "ATG" + "GCT" * 59 + "TAA" + "GGGG"
    orf = longest_orf(sequence)
    assert (orf["ORF strand"], orf["ORF start"], orf["ORF end"]) == ("+", 6, 188)
    assert orf["ORF length (aa)"] == 60
    assert orf["ORF complete"] is True


def test_orf_on_minus_strand():
    sequence = "CCCCC" + "ATG" + "GCT" * 59 + "TAA" + "GGGG"
    reverse = str(Seq(sequence).reverse_complement())
    orf = longest_orf(reverse)
    assert (orf["ORF strand"], orf["ORF start"], orf["ORF end"]) == ("-", 5, 187)


def test_no_orf_in_short_sequence():
    assert longest_orf("ATGAAATAA") is None


@pytest.mark.skipif(not EXAMPLE.exists(), reason="example data not downloaded")
def test_orfs_match_ncbi_annotation():
    expected = {"NM_000518.5": 147, "NM_000558.5": 142, "NM_000546.6": 393, "NM_007294.4": 1863}
    for record in SeqIO.parse(EXAMPLE, "fasta"):
        orf = longest_orf(str(record.seq))
        assert orf["ORF length (aa)"] == expected[record.id]


@pytest.mark.skipif(not EXAMPLE.exists(), reason="example data not downloaded")
def test_hbb_cds_coordinates():
    record = next(SeqIO.parse(EXAMPLE, "fasta"))
    orf = longest_orf(str(record.seq))
    assert (orf["ORF start"], orf["ORF end"]) == (51, 494)
    assert orf["ORF protein"].startswith("MVHLTPEEK")


# ---------- proteins ----------

@pytest.mark.skipif(not EXAMPLE.exists(), reason="example data not downloaded")
def test_hbb_molecular_weight_matches_uniprot():
    record = next(SeqIO.parse(EXAMPLE, "fasta"))
    hbb_protein = longest_orf(str(record.seq))["ORF protein"]
    # UniProt P68871 lists 15,998 Da for the 147-residue beta-globin
    assert protein_features(hbb_protein)["MW (kDa)"] == pytest.approx(16.0, abs=0.01)


# ---------- quality flags ----------

def test_flags_short_and_invalid():
    result = analyse_sequence("ACGT-ACGTACGT")
    assert "short" in result["Flags"]
    assert "invalid characters" in result["Flags"]


def test_internal_stop_flag():
    protein = "MVHLTPEEKSAVTALW*GKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAV*"
    assert "internal stop" in analyse_sequence(protein)["Flags"]