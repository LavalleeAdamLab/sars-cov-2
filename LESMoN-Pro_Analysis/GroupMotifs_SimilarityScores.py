import csv
import time
import resource
import pandas
import numpy
from math import comb


def import_significant_motifs_dictionary(file_name):
    """imports a csv file as a dictionary
    Columns of the dictionary:
        "Motif" -- 8-letter string
        "P-value" -- a double
        "NumOfProteins" -- an integer
        "ProtAccessions" -- a "|"-separated string
        "ProteinNames" -- a "|"-separated string
    """

    # initialize an empty dictionary
    significant_motifs_dictionary = dict()

    # opening the file
    with open(file_name) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            # for each row, retrieve the information from each column
            motif = row["Motif"]
            tpd = row["TPD"]
            p_value = row["P-value"]
            num_proteins = row["NumOfProteins"]
            protein_accessions = row["ProtAccessions"]
            protein_names = row["ProteinNames"]

            # update the dictionary; key is the motif and value is a list of the remaining columns
            significant_motifs_dictionary.update({motif: [tpd, p_value, num_proteins, protein_accessions, protein_names]})

    return significant_motifs_dictionary


def parse_identity_matrix(file_name, number_of_proteins):
    """
    parse a raw identity matrix and returns a panda dataframe with its values
    :param name_of_matrix:
    :return:
    """
    identity_matrix = pandas.DataFrame(columns=range(0, number_of_proteins))
    accession_name = []
    with open(file_name, "r") as file:
        reader = csv.reader(file, delimiter="\n")

        # pass the first 6 lines
        for _ in range(6):
            next(reader)

        # A counter for inserting each line to the DataFrame
        counter = 0

        # for each row
        for row in reader:
            a_row = row[0].split()  # Split the row by whitespace
            accession = a_row[1][3:9]  # Get the accession ID of the protein
            identity_matrix.loc[counter] = a_row[2:]  # Add the data to the identity matrix
            accession_name.append(accession)  # add the accession name eto the list

            counter += 1

    identity_matrix.columns = accession_name
    identity_matrix.index = accession_name

    identity_matrix = identity_matrix.fillna(0)

    return identity_matrix


def determine_fraction_of_homologous_proteins(percent_identity_matrix, protein_list, identity_threshold):
    """determines the number of homologous protein pairs in a list of proteins
    Arguments:
        percent_identity_matrix -- a DataFrame of percent identity values of pairwise sequence alignments
        protein_list -- a list of proteins which are being tested for fraction of homologous pairs
        threshold -- the fraction identity required to consider a pair of proteins homologous
    """

    # initialize a counter of homologous pairs
    homology_count = 0

    # for each protein in the list
    for index, protein_1 in enumerate(protein_list):
        # get a list of proteins to compare to the outer protein, to ensure each pairwise comparison is only done once
        remaining_proteins_list = protein_list[(index + 1):]
        for protein_2 in remaining_proteins_list:
            # if the percentage identity of the pairwise comparison is greater than the threshold, increment the count
            if float(percent_identity_matrix.loc[protein_1, protein_2]) >= identity_threshold:
                homology_count += 1

    # return the fraction of protein pairs that are homologous
    return homology_count / comb(len(protein_list), 2)


def produce_augmented_deg_to_aa_mapping(map_deg_to_aa):
    """for each degenerate character, add to the mapping any other degenerate characters it fully contains"""

    # copy the dictionary
    augmented_map_deg_to_aa = map_deg_to_aa.copy()

    # for each degenerate character
    for deg_char in map_deg_to_aa:
        # nested selection of another degenerate character
        for other_deg_char in map_deg_to_aa:
            # ensure not to add a character to its own entry
            if deg_char != other_deg_char:
                # if the inner degenerate character is a subset of the outer degenerate character, add to its mapping
                if set(map_deg_to_aa[other_deg_char]).issubset(map_deg_to_aa[deg_char]):
                    augmented_map_deg_to_aa[deg_char].append(other_deg_char)

    return augmented_map_deg_to_aa


def calculate_similarity_score_proteins(significant_motifs_dictionary, motif_1, motif_2):
    """calculates a similarity score for the proteins of two motifs"""

    # get the list of associated proteins for both motifs
    proteins_list_1 = significant_motifs_dictionary[motif_1][3].split("|")
    proteins_list_2 = significant_motifs_dictionary[motif_2][3].split("|")

    # calculate the number of proteins that are shared between both lists
    intersection = len(set(proteins_list_1) & set(proteins_list_2))

    # retrieve the length of the smaller of the two protein lists
    minimum_proteins_length = min(len(proteins_list_1), len(proteins_list_2))

    # calculate and return the similarity score
    # a similarity score of 1 means that the protein lists are identical
    return intersection / minimum_proteins_length


def calculate_distance_metric(similarity_score):
    """calculates a distance metric for the proteins of two motifs, based on a similarity score"""
    return (1 / similarity_score) - 1


def fill_lower_triangular_dataframe(df):
    """fills the upper triangular part of a lower triangular DataFrame"""
    df_nan_to_0 = df.fillna(0)  # fill NaN with 0
    df_transposed = df_nan_to_0.transpose()  # transpose
    df_filled = df_nan_to_0 + df_transposed  # add to the transpose
    numpy.fill_diagonal(df_filled.values, df.iloc[0, 0])  # restore original diagonal values

    return df_filled


# get start time
start_time = time.time()

# import dictionary
significant_motifs_dictionary = import_significant_motifs_dictionary(
    "krogan_significantMotifs_details_FDR0.05_proteinsMapped.txt")
print("Imported the dictionary.")

# parse a percentage identity matrix
percent_identity_matrix = parse_identity_matrix("clustalo-E20200609-163132-0117-39629537-p1m.pim", 195)

# define parameters for filtering based on homology
percent_identity_is_homologous = 0.25
accepted_fraction_homologous = 0.75


# filter the dictionary to proteins that are not homologous
filtered_significant_motifs_dictionary = dict()  # a new dictionary for motifs passing the homology test
for motif in significant_motifs_dictionary:
    if not determine_fraction_of_homologous_proteins(percent_identity_matrix,
                                                     significant_motifs_dictionary[motif][3].split("|"),
                                                     percent_identity_is_homologous) > accepted_fraction_homologous:
        filtered_significant_motifs_dictionary.update({motif: significant_motifs_dictionary[motif]})

print("Length before filtering: " + str(len(significant_motifs_dictionary)))
print("Length after filtering: " + str(len(filtered_significant_motifs_dictionary)))

# export the filtered dictionary to a file
with open("krogan_significantMotifs_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted.txt", 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Motif', 'TPD', 'P-value', 'NumOfProteins', 'ProtAccessions', 'ProteinNames'])
    for motif in filtered_significant_motifs_dictionary:
        # get the data associated to the motif from the full dictionary previously read from a file
        line_in_full_dictionary = significant_motifs_dictionary[motif]

        # extract each element of the list
        tpd = line_in_full_dictionary[0]
        pvalue = line_in_full_dictionary[1]
        num_proteins = line_in_full_dictionary[2]
        accessions = line_in_full_dictionary[3]
        names = line_in_full_dictionary[4]

        # write the data to the file
        tsv_writer.writerow([motif, tpd, pvalue, num_proteins, accessions, names])

# create the pandas DataFrame with similarity score
significant_motifs_list = list(set(filtered_significant_motifs_dictionary))  # list of motifs for row and column names
lower_triangular_similarity_matrix_list = []  # initialize the list representation of the DataFrame
for index_row, motif_row in enumerate(significant_motifs_list):  # loop through the rows
    similarity_scores = []  # initialize the list representation of one row of the DataFrame
    significant_motifs_for_columns = significant_motifs_list[:(index_row + 1)]  # row length for lower triangular matrix
    for motif_column in significant_motifs_for_columns:  # loop through the columns
        # calculate similarity score using the protein method and append it to the row list
        similarity_scores.append(
            calculate_similarity_score_proteins(filtered_significant_motifs_dictionary, motif_row, motif_column))
    lower_triangular_similarity_matrix_list.append(similarity_scores)  # append a row to the DF list representation
# create a DataFrame from its list representation
lower_triangular_similarity_df = pandas.DataFrame(lower_triangular_similarity_matrix_list,
                                                  significant_motifs_list, significant_motifs_list)
# fill the full matrix from the lower triangular portion
similarity_df = fill_lower_triangular_dataframe(lower_triangular_similarity_df)

# create the pandas DataFrame with distance metric
distance_df = calculate_distance_metric(similarity_df)

# output the distance matrix as a csv file
distance_df.to_csv("SignificantMotifs_DistanceMatrix_0.05FDR_0.25identity_0.75accepted.csv")

# print running time in seconds and peak memory usage
print(time.time() - start_time)
print(str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024 ** 2)))
