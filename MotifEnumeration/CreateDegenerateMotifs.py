def degenerate_substitutions(motif, mapping):
    """Compute the possible substitutions of degenerate characters at each position of a motif
    Arguments:
        motif -- a string of one-letter amino acid codes
        mapping -- a dictionary (char:list) of all the degenerate character substitutions for each amino acid
    Returns: a list of possible substitutions, each entry being a position followed by degenerate character (sep ":")
    """
    substitutions = []
    for i in range(len(motif)):
        for x in mapping[motif[i]]:
            substitutions.append(str(i) + ":" + x)
    return substitutions


def invert_dic_with_lists(original_dic):
    """Invert a dictionary mapping a value to a list"""
    inverse = {}
    for key in original_dic:
        for x in original_dic[key]:
            if x in inverse:
                inverse[x].append(key)
            else:
                inverse[x] = [key]
    return inverse


def perform_substitution(motif, substitution):
    """Perform a single substitution of a single amino acid for a degenerate character"""
    position = int(substitution.split(":")[0])
    char = substitution.split(":")[1]
    motif_as_list = list(motif)
    substituted_motif_as_list = motif_as_list
    substituted_motif_as_list[position] = char
    substituted_motif = "".join(substituted_motif_as_list)
    return substituted_motif


def determine_remaining_substitutions(substitutions, performed):
    """Identify substitutions on residues that have not yet been converted to a degenerate character"""
    new_substitutions = []
    position = performed.split(":")[0]
    for x in substitutions:
        if x.split(":")[0] != position:
            new_substitutions.append(x)
    return new_substitutions


def find_degenerate_motifs_single(motif, map_aa_to_deg):
    """Find all of the possible degenerate motifs with up to four degenerate characters for a given motif
    Arguments:
        motif -- a string of one-letter amino acid codes
        map_aa_to_deg -- a dictionary (char:list) of all the degenerate character substitutions for each amino acid
    Returns: a list of degenerate motifs, not including the original motif
    """

    # create a list to store all original and degenerate motifs
    degenerate_motifs = {motif}
    # compute all possible degenerate motifs for the given motif
    # compute a nested list of all possible substitutions for each position in the motif
    first_substitutions = degenerate_substitutions(motif, map_aa_to_deg)
    # begin code for first substitution
    for i in range(len(first_substitutions)):
        # perform substitution of the motif
        single_substituted_motif = perform_substitution(motif, first_substitutions[i])
        # add the motif to the list of degenerate motifs, if it is not a duplicate
        if single_substituted_motif not in degenerate_motifs:
            degenerate_motifs.add(single_substituted_motif)
        # determine the possible second substitutions
        second_substitutions = determine_remaining_substitutions(first_substitutions, first_substitutions[i])
        # begin code for second level substitution
        for j in range(len(second_substitutions)):
            # perform substitution of the motif
            double_substituted_motif = perform_substitution(single_substituted_motif, second_substitutions[j])
            # add the motif to the list of degenerate motifs, if it is not a duplicate
            if double_substituted_motif not in degenerate_motifs:
                degenerate_motifs.add(double_substituted_motif)
            # determine the possible third substitutions
            third_substitutions = determine_remaining_substitutions(second_substitutions, second_substitutions[j])
            # begin code for third level substitution
            for k in range(len(third_substitutions)):
                # perform substitution of the motif
                triple_substituted_motif = perform_substitution(double_substituted_motif, third_substitutions[k])
                # add the motif to the list of degenerate motifs, if it is not a duplicate
                if triple_substituted_motif not in degenerate_motifs:
                    degenerate_motifs.add(triple_substituted_motif)
                # determine the possible fourth substitutions
                fourth_substitutions = determine_remaining_substitutions(third_substitutions,
                                                                         third_substitutions[k])
                # begin code for fourth level substitution
                for m in range(len(fourth_substitutions)):
                    # perform substitution of the motif
                    quad_substituted_motif = perform_substitution(triple_substituted_motif, fourth_substitutions[m])
                    # add the motif to the list of degenerate motifs, if it is not a duplicate
                    if quad_substituted_motif not in degenerate_motifs:
                        degenerate_motifs.add(quad_substituted_motif)
    # remove the original motif from the degenerate motifs set
    degenerate_motifs.remove(motif)

    # return a tuple of the original motif and the set of degenerate motifs
    return motif, degenerate_motifs
