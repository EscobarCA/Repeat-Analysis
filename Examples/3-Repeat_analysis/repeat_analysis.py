# -*- coding: utf-8 -*-
"""
Module with functions that analyze repeat sequences in genes.

by Cristian Escobar
cristian.escobar.b@gmail.com

Last update 3/12/2022
"""

import sequence as sq
import json


def open_genes(file_name):
    """"Open genes annotation file in json format.

    Parameters
    ----------
    file_name : String
        JSON file name.

    Returns
    -------
    data_dict : Dictionary
        Dictionary with Gene instances {geneID: Gene_instance, ...}.
    """

    # Open annotation file in json format
    file = open(file_name, 'r')

    data = json.load(file)

    file.close()

    # Gene instances disctionary
    data_dict = {}

    chromosome = file_name[:-11]  # chromosome name from file name

    # read annotation
    for geneID in data:

        gene_data = data[geneID]

        name = gene_data['Name']
        seq_type = gene_data['seq_type']
        seq_range = gene_data['seq_range']
        orientation = gene_data['orientation']
        transcripts = gene_data['transcripts']

        # Create gene instance
        gene = sq.Gene(seq_type, seq_range, orientation)

        gene.set_geneID(geneID)
        gene.set_name(name)
        gene.set_chromosome(chromosome)

        for ID in transcripts:

            transcript_data = transcripts[ID]

            transcript_type = transcript_data['seq_type']
            transcript_range = transcript_data['seq_range']
            exons = transcript_data['exons']

            if transcript_type == 'mRNA':

                UTR5 = transcript_data['5UTR']
                CDSs = transcript_data['CDS']
                UTR3 = transcript_data['3UTR']

                # create mRNA instance
                transcript = sq.mRNA(
                    transcript_type, transcript_range, orientation)

                transcript.set_ID(ID)
                transcript.set_geneID(geneID)

                for cds_range in CDSs:

                    # Create CDS sequence instance
                    CDS = sq.Sequence('CDS', cds_range, orientation)

                    CDS.set_parent(ID)
                    CDS.set_geneID(geneID)

                    transcript.add_CDS(CDS)

                for UTR_range in UTR5:

                    # Create UTR sequence instance
                    UTR = sq.Sequence('5UTR', UTR_range, orientation)

                    UTR.set_parent(ID)
                    UTR.set_geneID(geneID)

                    transcript.add_5UTR(UTR)

                for UTR_range in UTR3:

                    # Create UTR sequence instance
                    UTR = sq.Sequence('3UTR', UTR_range, orientation)

                    UTR.set_parent(ID)
                    UTR.set_geneID(geneID)

                    transcript.add_3UTR(UTR)

            else:

                # Create transcript instance
                transcript = sq.Transcript(transcript_type, transcript_range,
                                           orientation)

                transcript.set_geneID(geneID)
                transcript.set_ID(ID)

            # Add exons to transcript
            for exon_range in exons:

                # Create exon sequence instance
                exon = sq.Sequence('exon', exon_range, orientation)

                exon.set_parent(ID)
                exon.set_geneID(geneID)

                transcript.add_exon(exon)

            # add transcript to gene
            gene.add_transcript(ID, transcript)

        # Add gene instance to dictionary
        data_dict[geneID] = gene

    return data_dict


def open_repeat_analysis(file_name, strand):
    """Open repeat files.


    Parameters
    ----------
    file_name :  String
        File containing information about the reepat position and number
        of repeats in tabular format.
    strand :  String
        Use '+' for sense strand and '-' for antisense.

    Returns
    -------
    data : List
        List of repeat Sequence instances.

    """

    file = open(file_name, 'r')

    lines = file.readlines()

    file.close()

    data = []

    for index in range(1, len(lines)):

        line = lines[index].split()

        position = int(line[0])

        end = int(line[1])

        repeat_num = line[2]

        seq_range = (position, end)

        # Create sequence instance
        repeat = sq.Sequence('repeat', seq_range, strand)

        repeat.set_ID(repeat_num)

        data.append(repeat)

    return data


def repeat_analysis(genes_dictionary, repeats_list):
    """Find repeat-gene matches.

    Method that finds repeats to gene matches. If repeat falls into
    a gene, the position of the repeat is evaluated to see whether matches
    a exon or intron. Method returns a list of genes with repeats.

    Parameters
    ----------
    genes_dictionary : TYPE
        Dictionary containing genes instances.
    repeats_list : List
        List of repeats [(repeat_position, repeat_number), ...].

    Returns
    -------
    genes_with_repeats : List
        List containing genes with repeats.
    """

    # List of gene instances that matched a repeat
    genes_with_repeats = []

    for geneID in genes_dictionary.keys():

        gene = genes_dictionary[geneID]

        gene_end = gene.get_seq_range()[1]

        for repeat in repeats_list:

            rep_start = repeat.get_seq_range()[0]

            # check if repeat is inside the gene
            match = gene.is_subsequence(repeat)

            if match:

                # analyze where the repeat falls
                if gene.get_transcripts() != {}:

                    gene.position_analysis(repeat)

                    gene.position_analysis_mRNA(repeat)

                    # Add geneID to result list

                    if gene not in genes_with_repeats:

                        genes_with_repeats.append(gene)

            elif not match and gene_end < rep_start:

                break


    return genes_with_repeats


def print_info_exon_intron(genes_with_repeats, filename):
    """Print detailed data for repeat-gene match.

    Parameters
    ----------
    genes_with_repeats : List
        List of genes with repeat match.
    filename : String
        Filename for output file.

    Returns
    -------
    None.
    """
    # counters
    total_exon_rep = 0
    total_intron_rep = 0

    # output file names
    exon_filename = '{}_{}.txt'.format(filename, 'exon')
    intron_filename = '{}_{}.txt'.format(filename, 'intron')

    # files
    exon_outfile = open(exon_filename, 'w')
    intron_outfile = open(intron_filename, 'w')

    # First line in file
    line = '{:10} {:10} {:7} {:<8} {:<10} {:<14} {:8}  {:10} {:10}   {:12}'.format('rep_start', 'rep_end',
                                                                                   'rep_num', 'distance',
                                                                                   'geneID', 'name',
                                                                                   'seq_type', 'gene_start',
                                                                                   'gene_end', 'chromosome')

    exon_outfile.write(line + '\n')
    intron_outfile.write(line + '\n')

    # add info for each gen in list
    for gene in genes_with_repeats:

        results_intron = gene.get_analysis_result('intron')

        total_intron_rep += len(results_intron)

        for result in results_intron:

            intron_outfile.write(result+'\n')

        results_exon = gene.get_analysis_result('exon')

        total_exon_rep += len(results_exon)

        for result in results_exon:

            exon_outfile.write(result+'\n')

    # print totals
    exon_outfile.write(
        'Total number of exons with repeats: {}'.format(total_exon_rep))
    intron_outfile.write(
        'Total number of introns with repeats: {}'.format(total_intron_rep))

    exon_outfile.close()
    intron_outfile.close()

    print('Printed file: {}'.format(exon_filename))
    print('Printed file: {}\n'.format(intron_filename))


def print_info_mRNA(genes_with_repeats, filename):
    """Print repeat analysis result for mRNA.


    Parameters
    ----------
    genes_with_repeats : List
        List of genes with repeat match.
    filename : String
        Filename for output file.

    Returns
    -------
    None.
    """

    # counters
    total_UTR5_rep = 0
    total_UTR3_rep = 0
    total_CDS_rep = 0

    # output file
    outfile = open(filename, 'w')

    # First line in file
    line = '{:10} {:10} {:7} {:<8} {:<10} {:<14} {:8}  {:10} {:10}   {:12}'.format('rep_start', 'rep_end',
                                                                                   'rep_num', 'distance',
                                                                                   'geneID', 'name',
                                                                                   'seq_type', 'gene_start',
                                                                                   'gene_end', 'chromosome')

    outfile.write(line + '\n')

    for gene in genes_with_repeats:

        if gene.has_mRNA_repeats():

            UTR5_results = gene.get_analysis_result_mRNA('5UTR')

            UTR3_results = gene.get_analysis_result_mRNA('3UTR')

            CDS_results = gene.get_analysis_result_mRNA('CDS')

            total_UTR5_rep += len(UTR5_results)
            total_UTR3_rep += len(UTR3_results)
            total_CDS_rep += len(CDS_results)

            for result in UTR5_results:

                outfile.write(result+'\n')

            for result in UTR3_results:

                outfile.write(result+'\n')

            for result in CDS_results:

                outfile.write(result+'\n')

    # print totals
    outfile.write(
        'Total number of 5pUTR with repeats: {}\n'.format(total_UTR5_rep))
    outfile.write(
        'Total number of 3pUTR with repeats: {}\n'.format(total_UTR3_rep))
    outfile.write('Total number of CDS with repeats: {}'.format(total_CDS_rep))

    outfile.close()


def distance_analysis(genes_with_repeats, seq_type, filename):
    """Performs distance analysis for repeat positions at introns.

    Method determines the distance from a repeat start to either end of an
    sequence. A positive distance indicates the repeat is closer to the sequence
    start, while a negative value indicates it is closer to the sequence 3' end
    All the distances are the collected as a histogram, so it couns how many
    repeats start at the same position. This data is then printed to a file.

    Parameters
    ----------
    genes_with_repeats : List
        List of gene instances with repeats
    seq_type : String
        Sequence type can be 'intron' or 'exon'.
    filename : String
        output file name.

    Returns
    -------
    None.
    """

    # open files
    info_file = open(filename, 'w')

    # first line
    info_file.write('{:<10} {:<4} {}\n'.format('Distance', 'Hits',
                                               'geneID-repeat position-rep num-Chromosome'))

    histogram_counter = {}

    for gene in genes_with_repeats:

        counter_dict = gene.get_distance_result(seq_type)

        for distance in counter_dict:

            if distance not in histogram_counter:

                histogram_counter[distance] = counter_dict[distance]

            else:

                histogram_counter[distance] += counter_dict[distance]

    # sort distances in the histogram data
    ordered_distances = list(histogram_counter.keys())

    ordered_distances.sort()

    # Print histogram data
    for distance in ordered_distances:

        hits = len(histogram_counter[distance])

        # line with info
        info_line = '{:<10} {:<4} '.format(distance, hits)

        for info in histogram_counter[distance]:

            info_line += info + ' '

        info_file.write(info_line+'\n')

    info_file.close()

    print('Distance analysis printed for ' + seq_type)


def main(gene_files, pos_files, neg_files):
    """Analize repeat data for genes.

    Main method that runs the repeat analysis.

    Parameters
    ----------
    gene_files : List
        List of file names containing gene information for a given chromosome
        in JSON format.
    pos_files : List
        List of repeats found in the sense strand (+).
    neg_files : List
        List of repeats found in the antisense strand (-).

    Returns
    -------
    None.
    """

    genes_with_repeats = []

    for gene_file, pos_file, neg_file in zip(gene_files, pos_files, neg_files):

        # open gene file
        genes_dictionary = open_genes(gene_file)

        print('Gene file {} opened \n'.format(gene_file))

        # open repeats files

        pos_repeats = open_repeat_analysis(pos_file, '+')
        neg_repeats = open_repeat_analysis(neg_file, '-')

        print('Repeat files opened: {} {}\n'.format(pos_file, neg_file))

        # Analysis of repeats on sense strand
        pos_gene_list= repeat_analysis(genes_dictionary, pos_repeats)
        # update counters
        genes_with_repeats += pos_gene_list


        # Analysis of repeats on anti-sense strand
        neg_gene_list = repeat_analysis(genes_dictionary, neg_repeats)

        # update counters
        genes_with_repeats += neg_gene_list


    # Print detailed results for repeats in intron and exons
    print_info_exon_intron(genes_with_repeats, '_Intron-Exon_results')

    # print analysis for mRNA
    print_info_mRNA(genes_with_repeats, '_mRNA_result.txt')

    # print distance analysis for intron
    distance_analysis(genes_with_repeats, 'intron',
                      '_Distance_analysis_intron.txt')


# Anotated files in JSON format
gene_files = ['NC_000024.10_genes.json']

# List of files with repeats found in sense strand
pos_files = ['Repeats_Y_+.txt']

# List of files with repeats found in antisense strand
neg_files = ['Repeats_Y_-.txt']


if __name__ == '__main__':

    main(gene_files, pos_files, neg_files)
