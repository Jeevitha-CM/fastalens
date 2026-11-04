from header_parser import parse_header


def test_trinity_header():
    info = parse_header(">TRINITY_DN1000_c115_g5_i1 len=247 path=[0:0-246]")
    assert info["Format"] == "Trinity"
    assert info["Gene"] == "TRINITY_DN1000_c115_g5"
    assert info["Isoform"] == "i1"
    assert info["len"] == "247"


def test_uniprot_header():
    info = parse_header(">sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2")
    assert info["Format"] == "UniProt"
    assert info["Accession"] == "P68871"
    assert info["Organism"] == "Homo sapiens"
    assert info["Gene"] == "HBB"
    assert info["Description"] == "Hemoglobin subunit beta"


def test_ncbi_mrna_header():
    info = parse_header(">NM_000518.5 Homo sapiens hemoglobin subunit beta (HBB), mRNA")
    assert info["Molecule"] == "mRNA"
    assert info["Gene"] == "HBB"
    assert info["Organism"] == "Homo sapiens"


def test_ncbi_protein_header_organism_in_brackets():
    # this was a real bug we found on Day 1
    info = parse_header(">NP_000509.1 hemoglobin subunit beta [Homo sapiens]")
    assert info["Molecule"] == "protein"
    assert info["Organism"] == "Homo sapiens"


def test_generic_header():
    info = parse_header(">my_sequence_1")
    assert info == {"ID": "my_sequence_1", "Format": "Generic"}