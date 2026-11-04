from Bio import SeqIO
from sequence_tools import longest_orf

for record in SeqIO.parse("data/human_mrna.fasta", "fasta"):
    orf = longest_orf(str(record.seq))
    print(record.id, orf["ORF strand"], orf["ORF start"], orf["ORF end"],
          orf["ORF length (aa)"], orf["ORF complete"])
    print("   protein starts:", orf["ORF protein"][:20])