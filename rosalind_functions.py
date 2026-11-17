
# Rosalind DNA Counting Nucleotides of DNA
def dna_counts(sequence):
    A = sequence.count('A')
    C = sequence.count('C')
    G = sequence.count('G')
    T = sequence.count('T')
    return f"{A} {C} {G} {T}"

# Rosalind RNA Converting DNA to RNA
def dna_to_rna(sequence):
    return sequence.replace('T', 'U')


# Rosalind REVC The strands complement
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    result = ""

    for base in sequence[::-1]:
        result += complement[base]

    return result
