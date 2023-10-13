# # Specify the file path
# file_path = r"C:\Users\Ali\Documents\Code\Python\Lectures\bioinformatics\data_dna_sequences.txt"

# # Open the file in read mode and read the content
# with open(file_path, 'r') as file:
#     file_content = file.read()

# # Convert the content to uppercase
# uppercase_content = file_content.upper()

# # Open the file again in write mode and write the uppercase content back to the file
# with open(file_path, 'w') as file:
#     file.write(uppercase_content)

# print("File content has been converted to uppercase.")

# Specify the file path
file_path = r"C:\Users\Ali\Documents\Code\Python\Lectures\bioinformatics\data_dna_sequences.txt"

# Initialize a counter
counter = 1

# Read the file content
with open(file_path, 'r') as file:
    file_content = file.read()

# Define the search string pattern
search_pattern = ">random sequence {} consisting of 100000 bases."

# Replace occurrences and increment counter
while True:
    search_string = search_pattern.format(counter)
    if search_string in file_content:
        # Replace with ">Rosalind_#" without adding a period
        file_content = file_content.replace(search_string, f">Rosalind_{counter}\n", 1)
        counter += 1
    else:
        break

# Write the modified content back to the file
with open(file_path, 'w') as file:
    file.write(file_content)

print(f"Replaced {counter - 1} occurrences.")

