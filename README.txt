# Repeat-Analysis

Overview
This project contains code used to:
1-	Find dinucleotide repeat sequences in genome sequences, for example in the human genome. 
2-	Analyze if those repeats are in genes and other sequences such as exons, introns, coding sequences (CDS), 
    and untranslated regions (UTR).

Using This Code:

Finding repeats (find_repeats.py).

This code requires a list of file names corresponding to the DNA sequence where the search will be performed. 
Files can be in either Fasta format or simple sequence files (no header). If a file contains more than one 
Fasta sequence, the program will produce an output file for each of those sequences independently. Since the 
code uses regular expression, this code needs two regular expression patterns, one for the sense strand and 
another for the antisense strand (see example). The code will output a file with repeats found at each sequence 
describing where the repeat sequence starts and ends, number of repeats, and the actual nucleotide sequence.

Analysis for Human genome was performed only for chromosomes 1 to 22, X and Y of assembly GRCh38.p14. Files can be downloaded at https://www.ncbi.nlm.nih.gov/genome/51?genome_assembly_id=1820449

Parsing genome annotation (get_genes.py and sequence.py)

Repeat analysis mainly focuses on DNA sequences that produce transcripts, therefore the Human genome annotation
was filtered to only consider sequences classified as ‘gene’. The get_genes.py code requires an annotation file
in gff format and the sequence.py file, which contains object definitions. This code will output the annotation
for each genomic sequence (chromosomes and other) in JSON format for easy access. The full repeat analysis only
considered the annotation for chromosomes 1 to 22, X and Y 
(see https://www.ncbi.nlm.nih.gov/genome/51?genome_assembly_id=1820449 for each chromosome RefSeq code)

Annotation of human genome can be downloaded here:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz

Repeat analysis (repeat_analysis.py, sequence.py)

This code will perform the actual analysis of repeat position in genes. As inputs, this code needs a list of files
containing repeats sequences for a given chromosome (for sense and antisense sequence), and the respective annotation. 
The output is divided into 3 files, one for repeats found in exons, introns and mRNA. In the case of mRNAs, the 
analysis returns information as to whether the repeat was in UTRs or CDSs. These files have information regarding 
the repeat, the gene where it was found, the exact sequence where it was found (intron, exon, UTRs or CDS) and the 
distance to either end of the sequence. For example, a positive signal indicates the repeat is closer to the 5’ end 
while a negative values correspond to repeats closer to the 3’ end.  An additional analysis is performed for introns
where distances for all introns and their respective repeat sequences are accumulated as a histogram, where repeats 
with the same distance  are put together (see examples).
