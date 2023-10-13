from Utility import NucleoTides
from Utility import DNA_Codons
from Utility import RNA_Codons
from collections import Counter

#Makes sure its valid DNA string
def validateSeq(dna_seq):
    tmp_string = dna_seq.upper()
    for nuc in tmp_string:
        if nuc not in NucleoTides:
            return False
        return tmp_string

def makeComplementDNARNA(seq):
    dnaComp = ""
    rnaSeq = ""
    for nuc in seq:
        if nuc == "A":
            dnaComp += "T"
            rnaSeq += "U"
        elif nuc == "T":
            dnaComp += "A"
            rnaSeq += "A"
        elif nuc == "G":
            dnaComp += "C"
            rnaSeq += "C"
        elif nuc == "C":
            dnaComp += "G"
            rnaSeq += "G"
    formatted_string = f"DNA String: 3' {dnaComp} 5', mRNA String: 3' {rnaSeq} 5'"

    return formatted_string

def countNucFreq(seq):
    tmpfreqDic = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpfreqDic[nuc] += 1

    final = f"A: {tmpfreqDic['A']}, C: {tmpfreqDic['C']}, G: {tmpfreqDic['G']}, T: {tmpfreqDic['T']}"
    return final

def findGC(seq):
    totalLength = len(seq)
    count = 0
    for nuc in seq:
        if nuc == "G" or nuc == "C":
            count += 1
    gcPercentage = round((count / totalLength) * 100, 2)
    return gcPercentage

def findGCsubset(seq, subset):
    result = []

    for i in range(0, len(seq), subset):
        section = seq[i:i+subset]
        gc_percentage = findGC(section)
        result.append(gc_percentage)
    return f"k={subset}: {result}"


def translate_dna_seq(seq, startPosition=0):
    protein_sequence = ""
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    
    # Iterate through the DNA sequence in groups of 3 (codons)
    for i in range(startPosition, len(seq), 3):
        codon = seq[i:i+3]
        
        # Check if it's a start codon
        if codon == start_codon:
            protein_sequence += "M"
        # Check if it's a stop codon
        elif codon in stop_codons:
            protein_sequence += "_"
        else:
            amino_acid = DNA_Codons.get(codon, "")
            protein_sequence += amino_acid
    return protein_sequence


def contains_start_codon(aminoAcidSeq):
    return "M" in aminoAcidSeq

def contains_stop_codon(aminoAcidSeq):
    return "_" in aminoAcidSeq

def generateProtein(aminoAcidSeq):
    sliced_sequence = ""
    if contains_start_codon(aminoAcidSeq):
        sliced_sequence = "M"
    else:
        return print("The amino acid sequence does not contain a start, therefore no protein is made.")

    position_of_M = aminoAcidSeq.find("M")
    sliced_sequence = aminoAcidSeq[position_of_M:]

    if contains_stop_codon(sliced_sequence):
        position_of_stop = sliced_sequence.find("_")
        sliced_sequence = sliced_sequence[:position_of_stop]
        return sliced_sequence
    else:
        return sliced_sequence


def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


def codon_usage(dna_seq, aminoacid_letter):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmpList = []
 
    for i in range(0, len(dna_seq) - 2, 3):
        if DNA_Codons[dna_seq[i:i + 3]] == aminoacid_letter:
            tmpList.append(dna_seq[i:i + 3])
    freqDict = dict(Counter(tmpList))

    totalWight = sum(freqDict.values())
    for dna_seq in freqDict:
        freqDict[dna_seq] = round(freqDict[dna_seq] / totalWight, 2)
    return freqDict


def reverse_complement(dna_seq, startPosition=0):

    reversed_sequence = dna_seq[::-1]
    complement_sequence = ""
    for i in range(startPosition, len(reversed_sequence)):
        nuc = reversed_sequence[i]
        if nuc == "A":
            complement_sequence += "T"
        elif nuc == "T":
            complement_sequence += "A"
        elif nuc == "G":
            complement_sequence += "C"
        elif nuc == "C":
            complement_sequence += "G"
    return translate_dna_seq(complement_sequence)

def generate_reading_frames(dna_seq):
    # Need to generate 6 reading frames, including the reverse complement.
    frames = []
    frames.append(translate_dna_seq(dna_seq,0))
    frames.append(translate_dna_seq(dna_seq,1))
    frames.append(translate_dna_seq(dna_seq,2))
    frames.append(reverse_complement(dna_seq,0))
    frames.append(reverse_complement(dna_seq,1))
    frames.append(reverse_complement(dna_seq,2))

    return frames
          
def allProteins(seq, startPos=0, endPos = 0, ordered = False ):
    if endPos > startPos:
        rfs = generate_reading_frames(seq[startPos:endPos])
    else:
        rfs = generate_reading_frames(seq)

    res = []
    for rf in rfs:
        proteins = proteins_from_rf(rf)
        for protein in proteins:
            res.append(protein)
    if ordered:
        return sorted(res,key=len, reverse=True)

    return res




# print(generate_reading_frames("AAACGTTTTGCAATCATAATATCGCACTGATGAGAGTGACTCCTAACTATAGGGCGAGTC"))
# print(translate_dna_seq("AAACGTTTTGCAATCATAATATCGCACTGATGAGAGTGACTCCTAACTATAGGGCGAGTC"))
# print(codon_usage("AAACGTTTTGCAATCATAATATCGCACTGATGAGAGTGACTCCTAACTATAGGGCGAGTC", "N"))