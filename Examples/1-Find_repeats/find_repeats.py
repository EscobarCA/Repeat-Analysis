# -*- coding: utf-8 -*-
"""Module containing functions that find repeat sequences using regular
expressions.

by Cristian Escobar
cristian.escobar.b@gmail.com

Last update 3/10/2022
"""

import time
import re


def open_sequence(file_name):
    """Open Nucleotide sequence file.

    Method that opens genomic files as either simple sequence or with multiple
    fasta sequences. It will return a dictionary with each fasta sequence
    in the file.

    Parameters
    ----------
    file_name : String
        Sequence file name.

    Returns
    -------
    data : Dict
        Dictionary containing each sequence in the genomic file..
    """

    print('Opening file...')

    # open file
    seq_file = open(file_name, 'r')

    sequence = seq_file.read()

    seq_file.close()

    # Dictionary containing each sequence
    data = {}

    if sequence[0] == '>':

        # divide the file by fasta header
        sequences = sequence.split('>')

        for seq in sequences:

            if seq != '':

                enter_index = seq.find('\n')

                first_line = seq[:enter_index]

                # get accesion code
                accesion = first_line.split()[0]

                # get nucleotide sequence after first \n
                fasta_seq = seq[enter_index:].replace('\n', '')

                fasta_seq = fasta_seq.upper()

                data[accesion] = fasta_seq

    else:

        sequence = sequence.replace('\n', '')
        sequence.upper()

        data[file_name] = sequence

    print('Sequence file opened')

    return data


def find_repeats(sequence, pattern):
    """Find repeat pattern in sequence.

    Function that finds repeats sequences using regular expressions and
    and the input pattern.

    Parameters
    ----------
    sequence : String
        DNA sequence string.
    pattern : String
        String used to search the sequence.

    Returns
    -------
    repeat_list : List
        List containing matches for the pattern sequence, each match
        corresponds to a tuple containing the sequence start, end and
        number of repeats.

    """

    # List containing repeats
    repeat_list = []

    # Find repeats ussing regular expressions and pattern
    hit_list = re.finditer(pattern, sequence)

    for repeat in hit_list:

        rep_start, rep_end = repeat.span()

        rep_num = (rep_end - rep_start)/2

        repeat_list.append((rep_start+1, rep_end, rep_num))

    return repeat_list


def print_repeats(sequence, repeat_list, filename):
    """Print repeats found to a file.


    Parameters
    ----------
    sequence : String
        DNA sequence string.
    repeat_list : List
        List with repeat data.
    filename : String
        File name for output file.

    Returns
    -------
    None.

    """

    file = open(filename, 'w')

    line = '{:<10} {:<10} {:<5}\n'.format('Start', 'End', 'n')

    file.write(line)

    for repeat in repeat_list:

        rep_start, rep_end, rep_num = repeat

        rep_seq = sequence[rep_start-1: rep_end]

        line = '{:<10} {:<10} {:<5}  {}\n'.format(rep_start, rep_end, rep_num,
                                                  rep_seq)

        file.write(line)

    file.close()


def get_histogram(*args):
    """Count repeats by length.

    Parameters
    ----------
    *args : List
        List containing repeat lists.

    Returns
    -------
    hist_dict : Dictionary
        Dictionary containing number of repeats with certain length.

    """

    hist_dict = {}

    for repeat_list in args:

        for repeat in repeat_list:

            rep_num = repeat[2]

            if rep_num not in hist_dict:

                hist_dict[rep_num] = 1

            else:

                hist_dict[rep_num] += 1

    return hist_dict


def print_totals(hist_dict, filename):
    """


    Parameters
    ----------
    hist_dict : Dictionary
        Histogram dictionary containing number of repeats with certain length.
    filename : String
        Output file name.

    Returns
    -------
    total : Int
        Total number of repeats found.

    """

    keys = list(hist_dict.keys())  # list of repeat lengths

    keys.sort()

    file = open(filename, 'w')

    file.write('rep_num  total\n')

    total = 0  # count repeats

    for rep_num in keys:

        repeats = hist_dict[rep_num]

        line = '{:7}  {}\n'.format(rep_num, repeats)

        file.write(line)

        total += repeats

    # print total number of repeats
    file.write('\nTotal: {}'.format(total))

    file.close()

    return total


def genomic_sequence_analysis(seq_file_list, pattern_dict, out_filename,
                              count_only=False):
    """Analize genome sequence for repeats.


    Parameters
    ----------
    seq_file_list : List
        List of file names corresponding to sequence files.
    pattern_dict : Dictionary
        Dictionary containing search patterns for search.
    out_filename : String
        File name for histogram analysis and totals.
    count_only : Bool, optional
        If true, the method won't print output_file containing repeats found.
        If False, the file won't be printed. This is useful to when the total
        number of repeats is needed. The default is False.

    Returns
    -------
    total : Int
        Total number of repeats.

    """
    # Dictionary containing sequence data
    data_dict = {}

    # List containing repeat list from all files
    repeat_results = []

    # Open files and collect sequences to data_dict
    for file_name in seq_file_list:

        time1 = time.time()

        data = open_sequence(file_name)

        data_dict.update(data)

        time2 = time.time()

        print(file_name, 'file opened in: ',  round(time2-time1, 2), 's')

    # Analysis of sequences to search for repeats
    for data in data_dict:

        sequence = data_dict[data]

        # Analyze each sequece with each pattern
        for key in pattern_dict.keys():

            repeat_list = find_repeats(sequence, pattern_dict[key])

            repeat_results.append(repeat_list)

            if not count_only:    # print only  if count_only is False

                repeats_file = 'Repeats_{}_{}.txt'.format(data, key)

                print(repeats_file)

                print_repeats(sequence, repeat_list, repeats_file)

    time2 = time.time()

    print('\nAnalysis in: ', round((time2-time1), 4), 's\n')

    # Calcualte number of repeats at all possible length found
    histogram = get_histogram(*repeat_results)

    total = print_totals(histogram, out_filename)

    return total


# Run main program

if __name__ == '__main__':

    # Sequence patterns used by regular expressions
    GT_pattern = 'G(TG){11,}T?'  # For GT 11.5 and larger
    AC_pattern = 'A?(CA){11,}C'  # For AC 11.5 and larger

    # Mapping of strand type to sequence pattern
    pattern_dict = {'+': GT_pattern, '-': AC_pattern}

    # List of sequence files to analyze
    chromosomes_list = ['Homo_sapiens.GRCh38.dna.chromosome.Y.fa']

    total = genomic_sequence_analysis(chromosomes_list, pattern_dict,
                                      'Histogram_totals.txt',
                                       count_only=False)
