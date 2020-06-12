from Bio import SeqIO


def read_fasta_file(fasta_file_name):
    """Get the contents of the fasta file"""
    krogan_protein_list = list(SeqIO.parse(fasta_file_name, 'fasta'))
    return krogan_protein_list


def identify_motifs(krogan_protein_list, window_size):
    """Read the Krogan proteins with a certain window size and create list of unique motifs"""

    # Empty dictionary for motifs
    motif_dictionary = dict()

    # Add the motifs to the dictionary
    for seq_record in krogan_protein_list:
        window_begins = 0
        window_ends = window_size

        while len(seq_record.seq[window_begins:window_ends]) == window_size:
            # Add the motif to the dictionary if it's unique. ".update" already checks the uniqueness
            motif_dictionary.update({str(seq_record.seq[window_begins:window_ends]): []})
            # Move the search window by one character from the start and end positions
            window_begins += 1
            window_ends += 1

    return motif_dictionary
