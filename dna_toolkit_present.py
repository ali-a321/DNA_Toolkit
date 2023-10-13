from dna_toolkit import *
from Utility import DNA_Sequence
import random

randomStr = "".join([random.choice(NucleoTides)
                for nuc in range(20)])

# print(makeComplementDNARNA(DNA_Sequence))
print(f"The Nucleotide count is {countNucFreq(DNA_Sequence)}")
print(f"The GC percentage is {findGC(DNA_Sequence)}%")
print(findGCsubset(DNA_Sequence, 100))
amino_acid_seq = translate_dna_seq(DNA_Sequence)
print(f"The amino acid sequence is: {amino_acid_seq}")
print(f"The first possible generated protein is: {generateProtein(amino_acid_seq)}")
print(f"The frequency of the specified codon is: {codon_usage(DNA_Sequence, 'A')}")
print(f"All possible proteins made are: ")
for prot in allProteins(DNA_Sequence, 0,0,True):
    print(f'{prot}')
