# -*- coding: utf-8 -*-
"""
Module for getting annotation for genes and psudogenes from gff annotation files

by Cristian Escobar
cristian.escobar.b@gmail.com

Last update 3/10/2022
"""

import sequence as sq
import json


def get_line_info(line):
    """Get gene information from line in data file.

    Parameters
    ----------
    line : String
        Line from gene anotation file.
    offset : Integer, optional
        Number used to offset the position of data in cases the columns
        do not match. The default is 0.

    Returns
    -------
    Instance of a gene, transcript or exon.

    None
        Returns None in case the sequence type is in the excluded list.

    """
    excluded = ( 'centromere', 'match', 'cDNA_match', 'recombination_feature',
                'sequence_feature', 'region', 'enhancer', 'D_loop')


    if 'Curated Genomic' in line:  # this sequence reference changes the column data order

        offset=1

    else:

        offset=0

    line = line.split()

    line_dict = {}



    # get seqeunce type
    line_dict['seq_type'] = line[2 + offset]

    if line_dict['seq_type'] in excluded:  # Return if  sequence type on the excluded list

        return None

    # get chromosome ID
    line_dict['chromosome'] = line[0]

    # Get sequence start and end
    line_dict['seq_range'] = (int(line[3+offset]), int(line[4+offset]))

    # Get sequence orientation (+: sense, -:antisense)
    line_dict['orientation'] = line[6+offset]

    # get general seqeunce description
    description = line[8+offset].split(';')

    line_dict['ID'] = description[0][3:]

    if 'Parent' in description[1]:

        line_dict['parent'] = description[1][7:]
        geneID = description[2]

        index = geneID.find('GeneID:')

        line_dict['geneID'] = geneID[index+7:].split(',')[0]
        line_dict['name'] = ''

    else:

        line_dict['name'] = description[2][5:]
        line_dict['parent'] = ''
        geneID = description[1]
        index = geneID.find('GeneID:')
        line_dict['geneID'] = geneID[index+7:].split(',')[0]

    return line_dict


def get_genes(file_name):

    # Open annotation file
    file = open(file_name, 'r')

    data_lines = file.readlines()

    file.close()

    print('Data file opened....')

    coding = ['mRNA', 'V_gene_segment', 'C_gene_segment', 'J_gene_segment',
              'D_gene_segment']

    data_dict = {}

    for line in data_lines:

        if line[0] == '#':

            continue

        line_dict = get_line_info(line)

        if line_dict == None:
            continue

        ID = line_dict['ID']
        geneID = line_dict['geneID']
        seq_type = line_dict['seq_type']
        seq_range = line_dict['seq_range']
        orientation = line_dict['orientation']
        chromosome = line_dict['chromosome']
        parent = line_dict['parent']

        if chromosome not in data_dict:

            data_dict[chromosome] = {}

            print(chromosome, ' started...')


        # select only for genes
        if 'ID=gene-' in line:

            gene = sq.Gene(seq_type, seq_range, orientation)

            gene.set_ID(ID)
            gene.set_geneID(geneID)
            gene.set_parent(parent)
            gene.set_name(line_dict['name'])
            gene.set_chromosome(chromosome)

            data_dict[chromosome][geneID] = gene

        # Select for general transcripts except mRNA
        elif geneID in data_dict[chromosome] and seq_type not in ['exon', 'CDS'] + coding:

            transcript = sq.Transcript(seq_type, seq_range, orientation)

            transcript.set_ID(ID)
            transcript.set_geneID(geneID)
            transcript.set_parent(parent)
            transcript.set_chromosome(chromosome)

            gene = data_dict[chromosome][geneID]
            gene.add_transcript(ID, transcript)


        elif geneID in data_dict[chromosome] and seq_type in coding:

            mrna = sq.mRNA(seq_type, seq_range, orientation)

            mrna.set_ID(ID)
            mrna.set_geneID(geneID)
            mrna.set_parent(parent)
            mrna.set_chromosome(chromosome)

            gene = data_dict[chromosome][geneID]
            gene.add_transcript(ID, mrna)


        elif geneID in data_dict[chromosome] and seq_type in ['exon', 'CDS']:

            sequence = sq.Sequence(seq_type, seq_range, orientation)
            sequence.set_ID(ID)
            sequence.set_geneID(geneID)
            sequence.set_parent(parent)
            sequence.set_chromosome(chromosome)

            gene = data_dict[chromosome][geneID]

            if seq_type == 'exon':

                gene.add_exon_to_transcript(parent, sequence)

            elif seq_type == 'CDS':

            # try:
                gene.add_CDS_to_mRNA(parent, sequence)

            # except:

                # print('Other coding region')
                # print(line)
                # print('\n')

    return data_dict


def print_data(data_dict):
    """Print gene data to json files per chromosome.

    Parameters
    ----------
    data_dict : Dict
        Dictionary containing gene data.

    Returns
    -------
    None.
    """

    for chromosome in data_dict:

        print(chromosome)
        print('Total genes:', len(data_dict[chromosome]))

        out_file = chromosome+'_genes.json'

        file = open(out_file, 'w')

        chromosome_dict = {}

        # Get gene info for each in dictionary format
        for geneID in data_dict[chromosome]:

            gene = data_dict[chromosome][geneID]

            if gene.get_seq_type() == 'gene':  # only collect genes

                gene_dict = gene.get_data_dict()

                chromosome_dict[geneID] = gene_dict

        # Dump chromosoem dictionary to Json file
        json.dump(chromosome_dict, file, indent=4)

        file.close()

        print('\n' + out_file + ' file printed')


def main(anotation_file):
    """Run main program.

    Parameters
    ----------
    anotation_file : String
        Name of file containing annotation data.

    Returns
    -------
    None.

    """

    gene_dictionary = get_genes(anotation_file)



    print_data(gene_dictionary)


if __name__ == '__main__':

    anotation_file = 'Annotation_fragment.txt'

    main(anotation_file)
