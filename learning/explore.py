from Bio import SeqIO

for record in SeqIO.parse("data/human_mrna.fasta", "fasta"):
    name = record.id
    length = len(record.seq)
    g = record.seq.count("G")
    c = record.seq.count("C")
    gc_percent = (g + c) / length * 100
    print(name, length, round(gc_percent, 2))