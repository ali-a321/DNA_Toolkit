# Initialize the variables as global
maxPercentage = 0
maxRosaId = ""

# Function to calculate GC content
def findGC(seq, rosalindId):
    totalLength = len(seq)
    count = 0
    for nuc in seq:
        if nuc == "G" or nuc == "C":
            count += 1
    gcPercentage = round((count / totalLength) * 100, 2)
    return gcPercentage, rosalindId

# Function to find the sequence with the maximum GC content
def maxGC(gcPercentage, rosalindId):
    global maxPercentage
    global maxRosaId
    
    if gcPercentage > maxPercentage:
        maxPercentage = gcPercentage
        maxRosaId = rosalindId

# Now you can use these functions while reading the file
# Open the file and read sequences, assuming your file format is in FASTA
file_path =  r"C:\Users\Ali\Documents\Code\Python\Lectures\bioinformatics\data_dna_sequences.txt"

id = ""
sequence = ""
dnaDic = {}
# Read the file and process sequences
with open(file_path, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            if id and sequence:
                gcPercentage, rosalindId = findGC(sequence, id)
                sequence_data = {"ID": rosalindId, "GC percentage": gcPercentage}
                dnaDic[id] = sequence_data
                maxGC(gcPercentage, rosalindId)
            id = line[1:]
            sequence = ""
        else:
            sequence += line

# Check the last sequence in the file
if id and sequence:
    gcPercentage, rosalindId = findGC(sequence, id)
    maxGC(gcPercentage, rosalindId)


output_file_path = r"C:\Users\Ali\Documents\Code\Python\Lectures\bioinformatics\data_dnaSequencesInfo.txt"
with open(output_file_path, 'w') as output_file:
    for id, data in dnaDic.items():
        output_file.write(f"{id}\n")
        output_file.write(f"GC percentage: {data['GC percentage']}%\n")
        output_file.write("\n") 

print(f"Data written to {output_file_path}")
print(f"The Rosalind ID with the highest GC content is {maxRosaId} with a GC content of {maxPercentage}%")