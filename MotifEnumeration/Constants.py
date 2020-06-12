# define the possible substitutions of degenerate characters for each amino acid
map_deg_to_aa = {
    "X": ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
    "B": ["D", "N"],
    "Z": ["E", "Q"],
    "J": ["I", "L"],
    "l": ["I", "V", "L"],
    "@": ["Y", "H", "W", "F"],
    "h": ["W", "F", "Y", "M", "L", "I", "V", "A", "C", "T", "H"],
    "o": ["S", "T"],
    "p": ["D", "E", "H", "K", "N", "Q", "R", "S", "T"],
    "t": ["A", "G", "C", "S"],
    "s": ["A", "G", "C", "S", "V", "N", "D", "T", "P"],
    "b": ["E", "F", "I", "K", "L", "M", "Q", "R", "W", "Y"],
    "+": ["K", "R", "H"],
    "-": ["D", "E"],
    "c": ["D", "E", "K", "R", "H"]
}

# define the length of the sliding window
window_size = 8

# fasta file name
fasta_file_name = "Krogan_Protein_database_REDUCED.fasta"

# type of motifs
real_or_shuffled = "real"

# output file names
output_file_name = "Krogan_motif_annotations_" + real_or_shuffled + ".tsv"
