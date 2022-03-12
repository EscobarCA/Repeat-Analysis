# -*- coding: utf-8 -*-
"""
Module with class definitions for Sequence, Gene, mRNA and Transcript

by Cristian Escobar
cristian.escobar.b@gmail.com

Last update 3/10/2022
"""


class Sequence():
    """DNA sequence Class"""

    def __init__(self, seq_type, seq_range, orientation):
        """Sequence class constructor.

        Parameters
        ----------
        seq_type : String
            Description of sequence type (gene, mrna, etc).
        seq_range : Tuple
            Tuple containing sequence start and end.
        orientation : String
            Sequence orientation can be '+' for sense strand and '-' for
            antisense strand.

        Returns
        -------
        None.
        """

        self.ID = ''
        self.seq_type = seq_type
        self.seq_range = seq_range
        self.parent = ''
        self.geneID = ''
        self.chromosome = ''
        self.orientation = orientation
        self.repeats = []

    def __eq__(self, other_sequence):
        """Compare if two sequences are equal.

        Method that evaluates if two sequences are equal based on the
        sequence start and end.

        Parameters
        ----------
        other_sequence : Sequence
            DESCRIPTION.

        Returns
        -------
        Bool
            Returns True if the tow sequences have the same sequence range.
        """

        other_seq_range = other_sequence.get_seq_range()

        return self.seq_range == other_seq_range

    def __str__(self):
        """Generate string representation of Sequence.

        Returns
        -------
        String
            Returns a string corresponding to the sequence range.
        """

        return str(self.seq_range)

    def add_repeat(self, repeat_seq):
        """Append a repeat sequence to repeat list.

        Parameters
        ----------
        repeat_seq : sequence
            Repeat sequence.

        Returns
        -------
        None.
        """

        self.repeats.append(repeat_seq)

    def get_chromosome(self):
        """Get chromosome ID.

        Returns
        -------
        String
            Chromosome identification.
        """

        return self.chromosome

    def get_distance(self, another_sequence):
        """Get another_sequence distance from to 5' or 3'.

        Method to get the sistance of another_sequence start to either 5' or 3'
        end of this sequence instance.

        Parameters
        ----------
        another_sequence : Sequence
            Sequence instance.

        Returns
        -------
        distance : Int
            Distance from start or end. A positive value indicates that
            position is closer to the 5' start of this sequence, while
            a negative value indicates that the other sequence is closer to
            the 3' end of this sequence.
        """
        # get other sequence range
        other_start, other_end = another_sequence.get_seq_range()

        # start of sequence depends on orientation
        if self.orientation == '-':

            end, start = self.seq_range
            position = other_end

        else:
            start, end = self.seq_range
            position = other_start

        # calculate distance to start and end
        dist_to_start = abs(position - start)+1

        dist_to_end = abs(end - position)+1

        # check which distance is smaller
        if dist_to_start <= dist_to_end:

            distance = dist_to_start

        else:

            distance = -dist_to_end

        return distance

    def get_geneID(self):
        """Get geneID.

        Returns
        -------
          String
              GeneID.
        """

        return self.geneID

    def get_ID(self):
        """Get sequence ID.

        Returns
        -------
        String
            Sequence ID.
        """

        return self.ID

    def get_orientation(self):
        """Get sequence orientation in the genome.

        Returns
        -------
        String
            Sequence orientation can be '+' for sense strand and '-' for
            antisense strand.
        """

        return self.orientation

    def get_parent(self):
        """Get parent sequence ID.

        Returns
        -------
        String
            parent sequence ID.
        """

        return self.parent

    def get_repeats(self):
        """Get repeat list.

        Method to get the list containing the repeats found in the
        sequence.

        Returns
        -------
        List
            List of repeat sequence instances.
        """

        return self.repeats

    def get_seq_range(self):
        """Get sequence range.

        Returns
        -------
        Tuple
            Tuple containing sequence start and end.
        """

        return self.seq_range

    def get_seq_type(self):
        """Get sequence type.

        Returns
        -------
        string
            Description of sequence type (mRNA, gene, etc).
        """

        return self.seq_type

    def in_sequence(self, position):
        """Return True if sequence position inside sequence range.

        Parameters
        ----------
        position : Int
            Position number in sequence.

        Returns
        -------
        Bool
            Returns True if positions is whithin the sequence start and
            end, False otherwise.

        """
        start, end = self.seq_range

        return start <= position and end >= position

    def is_subsequence(self, other_sequence):
        """Check if other sequence is subsequence of this sequence.

        To evaluate if the other sequence is subsequence of this sequence
        instance both sequences must have the same orientation and other
        sequence range must be within this sequence start and end.

        Parameters
        ----------
        other_sequence : Sequence
            Another Sequence instance.

        Returns
        -------
        Bool
            Returns True if other_sequence is subsequence of this instance.
        """

        other_start, other_end = other_sequence.get_seq_range()

        start, end = self.seq_range

        same_orientation = self.orientation == other_sequence.get_orientation()

        return start <= other_start and end >= other_end and same_orientation

    def set_orientation(self, orientation):
        """Set sequence orientation.

        Parameters
        ----------
        orientation : String
            Sequence orientation can be '+' for sense strand and '-' for
            antisense strand.

        Returns
        -------
        None.
        """

        self.orientation = orientation

    def set_chromosome(self, chromosome):
        """Set chromosome ID where sequence is found.

        Parameters
        ----------
        chromosome : String
            Chromosome ID

        Returns
        -------
        None.
        """
        self.chromosome = chromosome

    def set_ID(self, ID):
        """Set sequence ID.

        Parameters
        ----------
        ID : String
            Sequence identification.

        Returns
        -------
        None.
        """

        self.ID = ID

    def set_geneID(self, geneID):
        """Set geneID to sequence instance.

        Parameters
        ----------
        geneID : String
            Sequence geneID.

        Returns
        -------
        None.
        """

        self.geneID = geneID

    def set_parent(self, parent):
        """Set parent sequence ID.

        Parameters
        ----------
        parent : String
            parent sequence ID.

        Returns
        -------
        None.
        """

        self.parent = parent


class Transcript(Sequence):

    def __init__(self, seq_type, seq_range, orientation):
        """Transcript class constructor.

        Parameters
        ----------
        seq_type : String
            Description of sequence type (gene, mrna, etc).
        seq_range : Tuple
            Tuple containing sequence start and end.
        orientation : String
            Sequence orientation can be '+' for sense strand and '-' for
            antisense strand.

        Returns
        -------
        None.
        """

        super().__init__(seq_type, seq_range, orientation)

        self.exons = []

    def __str__(self):
        """String represenatation of Transcript.

        Returns
        -------
        line : String
            String represenatation of Transcript including sequence type, ID,
            Start and end.
        """
        # get seqeunce range
        start, end = self.seq_range

        line = '{:19} {:19} {:10} {:10} exons: '.format(self.seq_type,
                                                        self.ID, start, end)
        # Add exon range to line
        for exon in self.exons:

            start, end = exon.get_seq_range()

            line += '{}-{} '.format(start, end)

        return line

    def add_exon(self, exon):
        """Append exon sequence to list of exons

        Parameters
        ----------
        exon : Sequence
            Exon sequence instance.

        Returns
        -------
        None.
        """

        self.exons.append(exon)

    def get_exons(self):
        """Get exon list.

        Returns
        -------
        List
            List of exon sequences in the transcript.
        """

        return self.exons


class Gene(Sequence):

    def __init__(self, seq_type, seq_range, orientation):
        """Gene class constructor.

        Parameters
        ----------
        seq_type : String
            Description of sequence type (gene, mrna, etc).
        seq_range : Tuple
            Tuple containing sequence start and end.
        orientation : String
            Sequence orientation can be '+' for sense strand and '-' for
            antisense strand.

        Returns
        -------
        None.
        """

        super().__init__(seq_type, seq_range, orientation)

        self.orientation = orientation
        self.name = ''
        self.exons = []
        self.unspliced = []
        self.transcripts = {}
        self.repeats = []
        self.mRNA_repeats = []
        self.mRNAs = []

    def __str__(self):
        """String representation of gene instance.

        Returns
        -------
        line : String
            String representation of gene regarding geneID, name, sequence type,
            range, orientation and chromosome.
        """

        start, end = self.seq_range

        line = '{:10} {:12} {:8}  {:<10} {:<10} {} {:<12} '.format(self.geneID,
                                                                   self.name,
                                                                   self.seq_type,
                                                                   start, end,
                                                                   self.orientation,
                                                                   self.chromosome)

        return line

    def get_data_dict(self):
        """Get data in the Gene instance as dictionary.

        Method that collects all data in the Gene instance. This includes each
        Transcript and mRNA and their derived exons, UTRs and CDS. This data is
        used to generate the Json output file.

        Returns
        -------
        gene_dict : Dictionary
            Dictionary with all gene data.
        """

        gene_dict = {'geneID': self.geneID, 'Name': self.name,
                     'seq_type': self.seq_type, 'seq_range': self.seq_range,
                     'orientation': self.orientation,
                     'chromosome': self.chromosome, 'transcripts': {}}

        transcripts = self.get_transcripts()

        # fetch each transcript data
        for ID in transcripts:

            transcript = transcripts[ID]

            transcript_dict = {}

            transcript_dict['seq_type'] = transcript.get_seq_type()
            transcript_dict['ID'] = transcript.get_ID()
            transcript_dict['seq_range'] = transcript.get_seq_range()

            # get exons
            transcript_dict['exons'] = []

            for exon in transcript.get_exons():

                transcript_dict['exons'].append(exon.get_seq_range())

            # get info for mRNA
            if transcript.get_seq_type() == 'mRNA':

                if transcript.get_CDSs() != []:

                    # Calculate UTRs
                    transcript.get_5pUTR()
                    transcript.get_3pUTR()

                # get all sequences as list
                mrna_subseq = transcript.get_mRNA_subsequences()

                for seq in mrna_subseq:

                    seq_type = seq.get_seq_type()

                    seq_range = seq.get_seq_range()

                    if seq_type not in transcript_dict:

                        transcript_dict[seq_type] = [seq_range]

                    else:
                        transcript_dict[seq_type].append(seq_range)

            gene_dict['transcripts'][ID] = transcript_dict

        return gene_dict

    def get_exons(self):
        """Get exon list.

        Returns
        -------
        List
            List of exon sequences in the transcript.
        """

        return self.exons

    def get_unspliced(self):
        """Get exon-intron list.

        Returns
        -------
        List
            List of unique exon and intron sequences in the gene.
        """

        return self.unspliced

    def get_name(self):
        """Gete gene name.

        Returns
        -------
        String
            Gene name.
        """

        return self.name

    def set_name(self, name):
        """Set gene name.

        Parameters
        ----------
        name : String
            Gene name.

        Returns
        -------
        None.
        """

        self.name = name

    def get_transcripts(self):
        """Get transcript dictionary'

        Returns
        -------
        Dictionary
            Dictionary containing all transcript instances by ID.
        """

        return self.transcripts

    def get_mRNA_repats(self):
        """Get list of repeats found at mRNAs

        Returns
        -------
        List
            List of repeats found in mRNAs.
        """

        return self.mRNA_repeats

    def add_transcript(self, ID, transcript):
        """Add transcript instance to the gene.

        Parameters
        ----------
        ID : String
            Transcript ID.
        transcript : Transcript
            Transcript Instance to be added to the transcripts dictioanry and/or
            mRNAs list.

        Returns
        -------
        None.
        """

        self.transcripts[ID] = transcript

        seq_type = transcript.get_seq_type()

        if seq_type == 'mRNA':

            self.mRNAs.append(transcript.get_ID())

    def add_exon_to_transcript(self, transcript_ID, exon):
        """Add exon sequence to a transcript based on ID.

        Parameters
        ----------
        transcript_ID : String
            Transcript ID where exon will be added
        exon : Sequence
            Exon sequence instance.

        Returns
        -------
        None.
        """

        if transcript_ID in self.transcripts:

            self.transcripts[transcript_ID].add_exon(exon)

    def add_CDS_to_mRNA(self, mRNA_ID, CDS):
        """Add CDS sequence to a mRNAt based on ID.

        Parameters
        ----------
        mRNA_ID : String
            mRNA ID where exon will be added.
        CDS : Sequence
            CDS sequence instance.

        Returns
        -------
        None.
        """

        if mRNA_ID in self.transcripts:

            self.transcripts[mRNA_ID].add_CDS(CDS)

    def get_unique_exons(self, exon_list):
        """Create a list of unique_exons in a gene.

        This method creates a list of exons where exons that are subseuqence
        of another are removed. Usually these are exons smaller that the
        original. this method uses a recursive aproach.

        Parameters
        ----------
        exon_list : List
            List of possible exons in the gene from different transcripts.

        Returns
        -------
        List
            List of exons where most subsequences are removed.
        """

        if len(exon_list) == 0:

            return []

        if len(exon_list) == 1:

            return exon_list

        # compare two adjacent exons
        exon1 = exon_list[0]

        exon2 = exon_list[1]

        if exon1.is_subsequence(exon2):

            return [exon1] + self.get_unique_exons(exon_list[2:])

        elif exon2.is_subsequence(exon1):

            return [exon2] + self.get_unique_exons(exon_list[2:])

        else:

            return [exon1] + self.get_unique_exons(exon_list[1:])

    def get_sorted_exons(self):
        """Get list unique exons sorted by start position.

        Method creates a list of unique exons taken from all transcripts in
        the gene. Then it sorts this list based on exon start position.

        Returns
        -------
        exons : List
            List of unique Exons in the gene sorted by position.
        """

        exons = []

        for ID in self.transcripts:

            exons_in_transcript = self.transcripts[ID].get_exons()

            for exon in exons_in_transcript:

                if exon not in exons:

                    exons.append(exon)

        # sort by exon start
        exons.sort(key=lambda exon: exon.get_seq_range()[0])

        self.exons = exons

        return exons

    def has_mRNA_repeats(self):
        """Check if gene has repeats at mRNA


        Returns
        -------
        Bool
            Returns True if there repeats associated with mRNA.
        """

        return self.mRNA_repeats != []

    def has_mRNA(self):
        """Check whether gene has mRNA transcripts.

        Returns
        -------
        bool
            Returns True if there are mRNA transcripts associated to the gene.
        """

        if len(self.mRNAs) == 0:

            return False

        else:

            return True

    def transcrips_to_exons(self):
        """Create a list of unique exons from transcripts.

        This method creates a list of unique exons, where exon variants are
        removed by get_unique_exons method. Since this last method can't
        remove all subsequences, this method repeats the process until all
        subsequences are removed.

        Returns
        -------
        None.
        """

        # initial list of unique exons sorted by position
        exons = self.get_sorted_exons()

        difference = 1

        # run until there is no diference between exons and new exons
        while difference != 0:

            new_exons = self.get_unique_exons(exons)

            difference = len(exons) - len(new_exons)

            exons = new_exons

        self.exons = exons

    def get_exon_intron_list(self, exons):
        """Get a ordered list of exon and intron sequences.

        Method takes a list of unique exons and calculates intron sequences
        between exons. it returns a list or ordered exons and introns by
        possition. This method uses recursion to perform calculations


        Parameters
        ----------
        exons : List
            List of unique exons.

        Returns
        -------
        List
            Ordered list of exon and introns {exon1, intron1, exon2...].
        """

        if len(exons) == 1:
            self.unspliced += exons
            return exons

        # Get two adjacent exons and calculate intron range
        exon1_range = exons[0].get_seq_range()
        exon2_range = exons[1].get_seq_range()

        intron_range = (exon1_range[1] + 1, exon2_range[0] - 1)

        # Create an intron instance
        intron = Sequence('intron', intron_range, self.orientation)

        intron.set_ID(self.ID)
        intron.set_geneID(self.geneID)
        intron.set_chromosome(self.chromosome)
        intron.set_parent(self.ID)

        self.unspliced += [exons[0], intron]

        return [exons[0], intron] + self.get_exon_intron_list(exons[1:])

    def position_analysis(self, repeat):
        """Check where a repeat is located in the gene.

        Parameters
        ----------
        repeat : Sequence
            Repeat sequence instance.

        Returns
        -------
         bool
             Returns True if thereis a match
        """

        # define unique exons based on all different transcripts
        self.transcrips_to_exons()

        # define introns from exons

        if self.unspliced == []:

            subsequences = self.get_exon_intron_list(self.exons)

        else:

            subsequences = self.unspliced

        # check if repeat matches any subsequence
        for sequence in subsequences:

            if sequence.is_subsequence(repeat):

                sequence.add_repeat(repeat)

                self.repeats.append((repeat, sequence))

                return True

    def position_analysis_mRNA(self, repeat):
        """Check where a repeat is located in mRNA.

        Parameters
        ----------
        repeat : TYPE
            DESCRIPTION.

        Returns
        -------
        bool
            Returns True if thereis a match.
        """

        if self.has_mRNA():

            for mRNA_ID in self.mRNAs:

                mRNA = self.transcripts[mRNA_ID]

                # get 5' UTR + DCs, 3' UTRs
                subsequences = mRNA.get_mRNA_subsequences()

                for sequence in subsequences:

                    if sequence.is_subsequence(repeat):

                        sequence.add_repeat(repeat)

                        self.mRNA_repeats.append((repeat, sequence))

                        return True

    def get_analysis_result(self, seq_type):
        """Get result from repeat analysis for introns and exons.

        Method that returns a list of strings with information regarding
        repeats found in the gene.


        Parameters
        ----------
        seq_type : String
            Sequence type to filter result output, it can be 'exon' or 'intron'.

        Returns
        -------
        results : List
            List of results strings.
        """

        results = []

        for repeat, sequence in self.repeats:

            if seq_type == sequence.get_seq_type():

                # info about repeat
                repeat_start, repeat_end = repeat.get_seq_range()

                num_repeats = repeat.get_ID()

                distance = sequence.get_distance(repeat)

                line = '{:<10} {:<10} {:<7} {:<8} {} {:6} {}'.format(repeat_start,
                                                                     repeat_end,
                                                                     num_repeats,
                                                                     distance,
                                                                     self.__str__(),
                                                                     seq_type,
                                                                     sequence)

                results.append(line)

        return results

    def get_distance_result(self, seq_type):
        """Get distance result for all repeats in the specified sequence type.

        Parameters
        ----------
        seq_type : String
            Sequence type can be 'intron' or 'exon'.

        Returns
        -------
        counter : Dictionary
            Dictionary containing the distance to repeat info as key-value
            pairs.
        """

        counter = {}

        for repeat, sequence in self.repeats:

            if seq_type == sequence.get_seq_type():

                repeat_start, repeat_end = repeat.get_seq_range()

                num_repeats = repeat.get_ID()

                distance = sequence.get_distance(repeat)

                info = '{}-{}-{}-{}'.format(self.geneID, repeat_start,
                                            num_repeats, self.chromosome)

                if distance not in counter:

                    counter[distance] = [info]

                else:

                    counter[distance].append(info)

        return counter

    def get_analysis_result_mRNA(self, seq_type):
        """Get result from repeat analysis for mRNA sequence.

        Method that returns a list of strings with information regarding
        repeats found in mRNA.


        Parameters
        ----------
        seq_type : String
            Sequence type to filter result output, it can be 'CDS', '5UTR' or
            '3UTR.

        Returns
        -------
        results : List
            List of results strings
        """

        results = []

        for repeat, sequence in self.mRNA_repeats:

            if seq_type == sequence.get_seq_type():

                # info about repeat
                repeat_start, repeat_end = repeat.get_seq_range()

                num_repeats = repeat.get_ID()

                distance = sequence.get_distance(repeat)

                line = '{:<10} {:<10} {:<7} {:<8} {} {:6} {}'.format(repeat_start,
                                                                     repeat_end,
                                                                     num_repeats,
                                                                     distance,
                                                                     self.__str__(),
                                                                     seq_type,
                                                                     sequence)

                results.append(line)

        return results


class mRNA(Sequence):

    def __init__(self, seq_type, seq_range, orientation):

        super().__init__(seq_type, seq_range, orientation)

        self.sequence = ''

        self.CDSs = []
        self.exons = []
        self.UTR5 = []
        self.UTR3 = []

    def __str__(self):
        """String representation of mRNA sequence.

        Returns
        -------
        line : String
            String representation of mRNA sequence containing information
            about sequence type, ID, range, exon, CDS and UTRs.
        """

        start, end = self.seq_range

        line = '{:19} {:19} {:10} {:10} exons: '.format(self.seq_type,
                                                        self.ID, start, end)
        # Add exon info
        for exon in self.exons:

            start, end = exon.get_seq_range()

            line += '{}-{} '.format(start, end)

        cds_line = 'CDSs: '

        # Add  CDS info
        for CDS in self.CDSs:

            start, end = CDS.get_seq_range()

            cds_line += '{}-{} '.format(start, end)

        line += cds_line

        # Add  UTRs info
        UTR_line = '5pUTR: '

        for UTR in self.UTR5:

            start, end = UTR.get_seq_range()

            UTR_line += '{}-{} '.format(start, end)

        line += UTR_line

        UTR_line = '3pUTR: '

        for UTR in self.UTR3:

            start, end = UTR.get_seq_range()

            UTR_line += '{}-{} '.format(start, end)

        line += UTR_line

        return line

    def get_mRNA_subsequences(self):
        """Get lit of UTRs and CDS associated with mRNA.

        Returns
        -------
        List
            List of UTRS and CDS sequences in the mRNA.

        """

        return self.UTR5 + self.CDSs + self.UTR3

    def get_5pUTR(self):
        """Get 5' UTR sequences.

        Method to define 5' UTRs in a mRNA based on exon and CDS sequences.
        It takes into consideration gene sequence orientation.

        Returns
        -------
        None
        """

        for exon in self.exons:

            seq_range = ()

            exon_start, exon_end = exon.get_seq_range()

            if self.orientation == '+':

                CDS_start = self.CDSs[0].get_seq_range()[0]

                in_seq = exon.in_sequence(CDS_start)

                if CDS_start > exon_end and not in_seq:

                    seq_range = (exon_start, exon_end)

                elif in_seq:

                    if exon_start == CDS_start:

                        seq_range = (0, 0)

                    else:

                        seq_range = (exon_start, CDS_start-1)

            elif self.orientation == '-':

                CDS_start = self.CDSs[0].get_seq_range()[1]

                in_seq = exon.in_sequence(CDS_start)

                if CDS_start < exon_start and not in_seq:

                    seq_range = (exon_start, exon_end)

                elif in_seq:

                    if CDS_start == exon_end:

                        seq_range = (0, 0)
                    else:
                        seq_range = (CDS_start+1, exon_end)

            if seq_range != ():

                UTR = Sequence('5UTR', seq_range,  self.orientation)

                UTR.set_ID(exon.get_ID()+'-5UTR')

                UTR.set_geneID(exon.get_geneID())

                self.UTR5.append(UTR)

    def get_3pUTR(self):
        """Get 3' UTR sequences.

        Method to define 3' UTRs in a mRNA based on exon and CDS sequences.
        It takes into consideration gene sequence orientation.

        Returns
        -------
        None
        """

        for exon in self.exons:

            exon_start, exon_end = exon.get_seq_range()

            seq_range = ()

            if self.orientation == '+':

                CDS_end = self.CDSs[-1].get_seq_range()[1]

                in_seq = exon.in_sequence(CDS_end)

                if CDS_end < exon_start and not in_seq:

                    seq_range = (exon_start, exon_end)

                elif in_seq:

                    if CDS_end == exon_end:

                        seq_range = (0, 0)

                    else:
                        seq_range = (CDS_end+1, exon_end)

            elif self.orientation == '-':

                CDS_end = self.CDSs[-1].get_seq_range()[0]

                in_seq = exon.in_sequence(CDS_end)

                if CDS_end > exon_end and not in_seq:

                    seq_range = (exon_start, exon_end)

                elif in_seq:

                    if exon_start == CDS_end:

                        seq_range = (0, 0)

                    else:
                        seq_range = (exon_start, CDS_end-1)

            if seq_range != ():

                UTR = Sequence('3UTR', seq_range,  self.orientation)

                UTR.set_ID(exon.get_ID()+'-3UTR')

                UTR.set_geneID(exon.get_geneID())

                self.UTR3.append(UTR)

    def set_sequence(self, sequence):
        """Set mRNA DNA sequence.

        Parameters
        ----------
        sequence : String
            mRNA DNA sequence.

        Returns
        -------
        None.
        """

        self.sequence = sequence

    def get_sequence(self):

        return self.sequence

    def add_CDS(self, CDS):
        """Append CDS sequence to CDS list.

        Parameters
        ----------
        CDS : Sequence
            CDS sequence instance.

        Returns
        -------
        None.
        """

        self.CDSs.append(CDS)

    def add_5UTR(self, UTR):
        """Append 5' UTR sequence to 5' UTR list.

        Parameters
        ----------
        UTR : Sequence
            5' UTR sequence instance.

        Returns
        -------
        None.
        """

        self.UTR5.append(UTR)

    def add_3UTR(self, UTR):
        """Append 3' UTR sequence to 3' UTR list.

        Parameters
        ----------
        UTR : Sequence
            3' UTR sequence instance.

        Returns
        -------
        None.
        """

        self.UTR3.append(UTR)

    def get_CDSs(self):
        """Get list of CDS instances in the mRNA.

        Returns
        -------
        List
            list of CDS instances in the mRNA.
        """

        return self.CDSs

    def add_exon(self, exon):
        """Append exon sequence to list of exons.

        Parameters
        ----------
        exon : Sequence
            Exon sequence instance.

        Returns
        -------
        None.
        """

        self.exons.append(exon)

    def get_exons(self):
        """Get exon list.

        Returns
        -------
        List
            List of exon sequences in the mRNA.
        """

        return self.exons
