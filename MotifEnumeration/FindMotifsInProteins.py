def search_protein_sequences_for_motifs(krogan_protein_list, window_size, motif_dictionary):
    """Finds motifs using a sliding window and inputs them in dictionary mapping to proteins in which they appear"""

    for seq_record in krogan_protein_list:
        window_begins = 0
        window_ends = window_size

        while len(seq_record.seq[window_begins:window_ends]) == window_size:
            # Making sure the proteins accession id are unique
            if seq_record.id[3:9] not in motif_dictionary[seq_record.seq[window_begins:window_ends]]:
                motif_dictionary[seq_record.seq[window_begins:window_ends]].append(seq_record.id[3:9])
            # Move the search window by one character from the start and end positions
            window_begins += 1
            window_ends += 1
    return motif_dictionary
