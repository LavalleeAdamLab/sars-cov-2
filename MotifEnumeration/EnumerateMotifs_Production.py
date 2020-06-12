import MotifDiscovery
import Constants
import CreateDegenerateMotifs
import FindMotifsInProteins
from multiprocessing import Pool
import Output_Production
import time


def update_motif_to_protein_dictionary(result):
    # get the original motif and the set of degenerate motifs from the return tuple
    original_motif = result[0]
    degenerate_motifs = result[1]

    # get the proteins associated to the original motif
    associated_proteins = motif_to_proteins_dictionary[original_motif]
    associated_proteins_list = associated_proteins.split("|")

    # for each degenerate motif
    for motif in degenerate_motifs:
        # if the motif is already in the dictionary
        if motif in motif_to_proteins_dictionary:
            # add each protein to the list of proteins in that entry, without duplicates
            for protein in associated_proteins_list:
                if protein not in motif_to_proteins_dictionary[motif]:
                    motif_to_proteins_dictionary[motif] = motif_to_proteins_dictionary[motif] + "|" + protein
        # if the motif is not yet in the dictionary
        else:
            # add the motif to the dictionary with all of its associated proteins
            motif_to_proteins_dictionary[motif] = associated_proteins

    # increment the counter
    global counter
    counter += 1

    # print a descriptive dictionary update statement
    print("Just updated dictionary for degenerate motifs of " + original_motif + ", which is original motif #" + str(counter))


# PART 1: GETTING THE ORIGINAL MOTIFS AND MAPPINGS
# read fasta file of proteins
krogan_protein_list = MotifDiscovery.read_fasta_file(Constants.fasta_file_name)

# find original motifs from protein sequence sliding windows
window_size = Constants.window_size
motif_dictionary = MotifDiscovery.identify_motifs(krogan_protein_list, window_size)

# count motifs in proteins
motif_to_proteins_dictionary = FindMotifsInProteins.search_protein_sequences_for_motifs(krogan_protein_list,
                                                                                        window_size,
                                                                                        motif_dictionary)

# in the dictionary, change lists to strings separated by "|"
for key in motif_to_proteins_dictionary:
    motif_to_proteins_dictionary[key] = "|".join(motif_to_proteins_dictionary[key])

# get the list of motifs as a set in preparation for finding all degenerate motifs
motifs = set(motif_dictionary)


# PART 2: GETTING THE DEGENERATE MOTIFS AND MAPPINGS
# get start time for finding the degenerate motifs
start_time = time.time()

# define the possible substitutions of degenerate characters for each amino acid
map_aa_to_deg = CreateDegenerateMotifs.invert_dic_with_lists(Constants.map_deg_to_aa)

# set a counter to keep track of the original motifs
counter = 0

# find all of the possible degenerate motifs, parallel implementation using multiprocessing pool
with Pool(6) as p:
    [p.apply_async(CreateDegenerateMotifs.find_degenerate_motifs_single,
                   args=(motif, map_aa_to_deg),
                   callback=update_motif_to_protein_dictionary)
     for motif in motifs]
    p.close()
    p.join()

# get the final time for part 2
degenerate_motif_time = time.time() - start_time
print("It took " + str(degenerate_motif_time) + " seconds (" + str(degenerate_motif_time / 60) + " minutes) " +
      "to find the degenerate motifs and their associated proteins.")


# PART 3: OUTPUTTING THE ANNOTATION TABLE
# get start time for outputting the annotation table
start_time = time.time()

# output the annotation table
Output_Production.output_annotations_csv("name_to_uniprot_mapping.txt", motif_to_proteins_dictionary, Constants.output_file_name)

# get the final time for part 3
output_time = time.time() - start_time
print("It took " + str(output_time) + " seconds (" + str(output_time / 60) + " minutes) " + 
      "to output the annotation table.")
