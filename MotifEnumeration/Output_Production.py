import csv


def output_annotations_csv(uniprot_mapping_file_name, motif_to_proteins_dictionary, output_name):
    """Outputs the motif to protein to a tsv file"""

    # A dictionary with the mapping between gene name and their IDs
    gene_name_id_dic=dict()
    # read the CSV file with the name to uniprot ID mapping
    # Uniprot accession are the keys and the names are the values
    with open(uniprot_mapping_file_name, mode='r') as csv_file:
        mapping_table = csv.reader(csv_file, delimiter="\t")
        next(mapping_table)  # Skips the header
        for row in mapping_table:
            gene_name_id_dic.update({row[0]: row[1]})

    # write the tsv file mapping motifs to gene IDs and gene symbols
    with open(output_name, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Motif', 'Placeholder', 'Placeholder', 'Placeholder', 'Placeholder',
                             'Placeholder', 'Placeholder', 'Gene IDs', 'Gene Symbols'])
        for elem in motif_to_proteins_dictionary:
            name_list = []
            # Search the mapping dictionary to identify the names fo the genes
            for elements in motif_to_proteins_dictionary[elem].split("|"):
                name_list.append(gene_name_id_dic[elements])

            # Format the Gene Names appropriately
            name_formatted = '|'.join([str(e) for e in name_list])
            motif = str(elem)
            # Gene IDs are already formatted appropriately
            gene_ids = motif_to_proteins_dictionary[elem]
            tsv_writer.writerow([motif, '', '', '', '', '', '', gene_ids, name_formatted])


def output_motifs_txt(all_motifs, output_name):
    # write all of the motifs to a file
    with open(output_name, 'wt') as out_file:
        for motif in all_motifs:
            out_file.write(motif + "\n")
