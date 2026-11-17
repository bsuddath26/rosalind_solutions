
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

# Rosalind PROT Translating RNA into Protein
def translate_rna(rna_sequence):
    code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '-', 'TAG': '-',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '-', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    dna_sequence = rna_sequence.replace('U', 'T')
    protein = ""
    i = 0

    while i < len(dna_sequence) - 2:
        codon = dna_sequence[i:i+3]

        if code[codon] == '-':
            return protein

        protein += code[codon]
        i += 3

    return protein


# Rosalind SUBS Finding a Motif in DNA
def find_motif(sequence, motif):
    positions = ""

    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i+len(motif)] == motif:
            positions += str(i + 1) + " "

    return positions.strip()

# Rosalind PRTM Calculating the Mass of Protein
def prot_mass(prot):
    mass = {
        'A': '71.03711',  'C': '103.00919', 'D': '115.02694',
        'E': '129.04259', 'F': '147.06841', 'G': '57.02146',
        'H': '137.05891', 'I': '113.08406', 'K': '128.09496',
        'L': '113.08406', 'M': '131.04049', 'N': '114.04293',
        'P': '97.05276',  'Q': '128.05858', 'R': '156.10111',
        'S': '87.03203',  'T': '101.04768', 'V': '99.06841',
        'W': '186.07931', 'Y': '163.06333'
    }

    total = 0

    for aa in prot:
        total = total + float(mass[aa])

    return total


# Rosalind REVP Locating Restriction Sites
def restriction_sites(fasta_data):
    sequence = ""

    for line in fasta_data.strip().split('\n'):
        line = line.strip()
        if not line.startswith(">"):
            sequence += line

    results = ""

    for i in range(len(sequence)):
        for length in [4, 6, 8, 10, 12]:
            part = sequence[i:i+length]

            if len(part) == length:
                if part == reverse_complement(part):
                    results += f"{i+1} {length}\n"

    return results.strip()
