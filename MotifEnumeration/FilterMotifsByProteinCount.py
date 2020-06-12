import csv
import time


ORIGINAL_FILE_NAME = "Krogan_motif_annotations_real.tsv"
FILTERED_FILE_NAME = "Krogan_motif_annotations_real_FILTERED.tsv"
NUMBER_PROTEIN_THRESHOLD = 3


def main():
    read_filter_write(ORIGINAL_FILE_NAME,FILTERED_FILE_NAME, NUMBER_PROTEIN_THRESHOLD)


def read_filter_write(input_file_name, output_file_name, threshold):
    """
    Function that reads the original motif file line by line and
    only keeps the motifs with more than NUMBER_PROTEIN_THRESHOLD proteins
    the motifs that pass the threshold will be written into a new file
    :param input_file_name: Input file name
    :param output_file_name: Output file name
    :param threshold: the threshold for filtering the motifs
    """
    with open(output_file_name, mode='w',newline='', buffering=131072) as output_file:  # 128 KiB buffer
        with open(input_file_name, mode='r') as input_file:

            reader = csv.reader(input_file, delimiter="\t")  # The reader object
            writer = csv.writer(output_file, delimiter="\t")  # The writer object

            writer.writerow(next(reader))  # Write the header into the output file

            # For each row, if there are equal or more proteins than the threshold present
            # write them into the output file
            for row in reader:
                if len(row[7].split("|")) >= threshold:
                    writer.writerow(row)


if __name__ == '__main__':
    start_time = time.time()
    main()
    print("Took " + str(time.time() - start_time) + " seconds to filter the file.")
