import csv
import math
import random
import regex
import Constants
import MotifDiscovery


def import_representative_motif_dictionary(file_name):
    """import a csv file with one or more representative motifs for each cluster of motifs
    Columns of the file:
        Cluster -- an integer corresponding to the cluster number
        Motifs -- the motif(s) with the lowest p-value, separated by "|"
    """

    # initialize an empty dictionary
    representative_motif_dictionary = dict()

    # opening the file
    with open(file_name) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            # for each row, retrieve the information from each column
            cluster = row["cluster_number"]
            motifs = row["motif(s)"]

            # update the dictionary; key is the motif and value is a list of the remaining columns
            representative_motif_dictionary.update({cluster: motifs.split("|")})

    return representative_motif_dictionary


def count_degenerate_substitutions(motif, map_deg_to_aa):
    """calculate the sum of the numbers of possible substitutions for each degenerate character in a motif sequence"""
    # initialize a sum
    sum_of_substitutions = 0

    # loop through each character and add the number of substitutions at that character to the sum
    for character in motif:
        if character in map_deg_to_aa:
            sum_of_substitutions += len(map_deg_to_aa[character])

    return sum_of_substitutions


def choose_representative_motif(motif_list, map_deg_to_aa):
    """from a list of motifs, return the one that is the least degenerate, based on the number of degenerate characters
     and, if this does not determine the least degenerate motif, the sum of the number of possible substitutions for
     each degenerate character; if two or more are identical in both measures, choose one randomly"""

    # keep track of the least degenerate motif
    least_degenerate_motif = []

    # first, try to get a least degenerate motif using the number of degenerate characters; keep track of that number
    lowest_number_of_degenerate_characters = math.inf

    for motif in motif_list:
        # count the number of degenerate characters
        number_of_degenerate_characters = 0
        for character in motif:
            if character in map_deg_to_aa:
                number_of_degenerate_characters += 1

        # if this is the motif with the lowest number of degenerate characters so far, update the variables
        if number_of_degenerate_characters < lowest_number_of_degenerate_characters:
            lowest_number_of_degenerate_characters = number_of_degenerate_characters
            least_degenerate_motif = [motif]
        # if it's tied with the current motif, add it to the list
        elif number_of_degenerate_characters == lowest_number_of_degenerate_characters:
            least_degenerate_motif.append(motif)

    # if there is a single motif with the lowest number of degenerate characters, return it as a string
    if len(least_degenerate_motif) == 1:
        print("1")
        return least_degenerate_motif[0]

    # if there is more than one motif, narrow it down using sum of substitutions for all degenerate characters
    # keep track of that sum
    lowest_number_of_degenerate_substitutions = math.inf

    least_degenerate_motif_candidates = least_degenerate_motif.copy()
    least_degenerate_motif = []

    for motif in least_degenerate_motif_candidates:
        # get the sum of substitutions for all degenerate characters
        sum_of_substitutions = count_degenerate_substitutions(motif, map_deg_to_aa)

        # if this is the motif with the lowest sum so far, update the variables
        if sum_of_substitutions < lowest_number_of_degenerate_substitutions:
            lowest_number_of_degenerate_substitutions = sum_of_substitutions
            least_degenerate_motif = [motif]
        # if it's tied with the current motif, add it to the list
        elif sum_of_substitutions == lowest_number_of_degenerate_substitutions:
            least_degenerate_motif.append(motif)

    # if there is a single motif with the lowest number of degenerate substitutions, return it as a string
    if len(least_degenerate_motif) == 1:
        print("2")
        return least_degenerate_motif[0]

    # if there is no single motif identified, return a random one as a string
    print("3")
    return random.sample(least_degenerate_motif, 1)[0]


def motif_to_regex(deg_motif):
    regex = ""
    for char in deg_motif:
        if char in Constants.map_deg_to_aa:
            regex = regex + "[" + "".join([i for i in Constants.map_deg_to_aa[char]]) + "]"
        else:
            regex = regex + char
    return str(regex)


def find_corresponding_original_motifs(degenerate_motif_regex, krogan_protein_list):
    """find all original motifs corresponding to a degenerate motifs, returned as a list"""

    # initialize a list of all matches to the degenerate motif
    corresponding_original_motifs = []

    # loop through each protein and find all occurrences of the degenerate motif
    for seq_record in krogan_protein_list:
        corresponding_original_motifs += regex.findall(degenerate_motif_regex, str(seq_record.seq))

    # return the list of corresponding original motifs as a "|"-separated list
    return corresponding_original_motifs


def output_corresponding_original_motifs(file_name, representative_motif_dictionary,
                                         motif_to_corresponding_original_motifs_dictionary):
    with open(file_name, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['cluster_number', 'representative_motif', 'original_motifs'])
        for cluster_id in representative_motif_dictionary:
            # get the representative and original motifs corresponding to the cluster
            representative_motif = representative_motif_dictionary[cluster_id]
            original_motifs = motif_to_corresponding_original_motifs_dictionary[cluster_id]

            # write the data to the file, with original motifs as a "|"-separated string
            tsv_writer.writerow([cluster_id, representative_motif, "|".join(original_motifs)])


# import the dictionary of representative motifs per cluster
representative_motif_dictionary = import_representative_motif_dictionary("clusters_sig_motifs_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted_9clusters_completelinkage.tsv")

# make a new dictionary of corresponding original motifs for each cluster
motif_to_corresponding_original_motifs_dictionary = dict()

# import a fasta file of protein sequences
krogan_protein_list = MotifDiscovery.read_fasta_file(Constants.fasta_file_name)

# loop through the clusters in the dictionary
for cluster_id in representative_motif_dictionary:
    # get the representative motif(s) of the cluster
    motifs = representative_motif_dictionary[cluster_id]

    # get a single motif to represent the cluster
    if len(motifs) == 1:
        # if there is only one motif, it is the representative motif
        motif = motifs[0]
    else:
        # if there is more than one motif, a representative motif is chosen
        motif = choose_representative_motif(motifs, Constants.map_deg_to_aa)

    # replace the motif list in the dictionary with the string representing the unique representative motif
    representative_motif_dictionary[cluster_id] = motif

    # turn the representative motif into a regular expression
    motif_regex = motif_to_regex(motif)

    # find the original motifs corresponding to the representative degenerate motif
    corresponding_original_motifs = find_corresponding_original_motifs(motif_regex, krogan_protein_list)

    # add the list of corresponding original motifs to the dictionary as the value for the cluster_id key
    motif_to_corresponding_original_motifs_dictionary.update({cluster_id: corresponding_original_motifs})

# output the dictionary to a file
output_corresponding_original_motifs("clusters_corresponding_original_motifs_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted_9clusters_completelinkage.tsv", representative_motif_dictionary,
                                     motif_to_corresponding_original_motifs_dictionary)
