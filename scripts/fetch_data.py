# Step 1: bring in the "Entrez" tool from Biopython
from Bio import Entrez

# Step 2: tell NCBI who we are (they ask everyone for an email)
Entrez.email = "your.email@gmail.com"

# Step 3: the list of genes we want (their NCBI ID numbers)
genes = ["NM_000518", "NM_000558", "NM_000546", "NM_007294"]

# Step 4: ask NCBI for them, in FASTA format
answer = Entrez.efetch(db="nucleotide", id=genes, rettype="fasta", retmode="text")
fasta_text = answer.read()
answer.close()

# Step 5: save what we got into a file
file = open("data/human_mrna.fasta", "w")
file.write(fasta_text)
file.close()

# Step 6: count the ">" signs to check we got 4 genes
print("Saved", fasta_text.count(">"), "sequences")