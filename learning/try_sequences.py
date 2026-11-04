from Bio import Entrez, SeqIO
from header_parser import parse_header
from sequence_tools import analyse_sequence

Entrez.email = "jeevithacm21@gmail.com"

answer = Entrez.efetch(db="protein", id="NP_000509", rettype="fasta", retmode="text")
protein_record = SeqIO.read(answer, "fasta")
answer.close()

print(parse_header(protein_record.description))
print(analyse_sequence(str(protein_record.seq)))

for record in SeqIO.parse("data/human_mrna.fasta", "fasta"):
    print(record.id, analyse_sequence(str(record.seq)))