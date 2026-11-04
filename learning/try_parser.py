from header_parser import parse_header

examples = [
    ">TRINITY_DN1000_c115_g5_i1 len=247 path=[0:0-246]",
    ">sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2",
    ">NM_000518.5 Homo sapiens hemoglobin subunit beta (HBB), mRNA",
    ">my_sequence_1",
]

for header in examples:
    print(parse_header(header))