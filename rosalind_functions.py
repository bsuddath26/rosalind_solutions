
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

# Rosalind GC Calculating proportion of GC nucleotides
def gc_content(fasta_data):
    sequences = {}
    current_id = ""

    for line in fasta_data.strip().splitlines():
        line = line.strip()

        if line.startswith(">"):
            current_id = line[1:]
            sequences[current_id] = ""
        else:
            sequences[current_id] += line

    best_id = ""
    best_gc = -1

    for seq_id in sequences:
        sequence = sequences[seq_id]
        G = sequence.count('G')
        C = sequence.count('C')
        gc = (G + C) / len(sequence) * 100

        if gc > best_gc:
            best_gc = gc
            best_id = seq_id

    return f"{best_id}\n{best_gc:.6f}"


# Rosalind HAMM Counting Point Mutations
def mutation_counter(seq1, seq2):
    mutations = 0

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mutations += 1

    return mutations

