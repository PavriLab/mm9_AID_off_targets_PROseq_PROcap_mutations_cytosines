#!/usr/bin/env python3

'''
Description:
    the python script takes as input a list of genes in bed format and the name of the output file.
    It will then iterate over the genes and print out to a new file the sites that fall into one of the
    two motifs (RGYW and WRCY), along with the name of the gene, the strand, the base, and the motif itself.
    For the script to work, you need the fasta files of the mouse mm9 genome placed in a specified folder named
    "mm9_genome" two directories above the directory of the script (see 'folder_structure') and the files of the
    chromosomes should be called "Mus_musculus.NCBIM37.67.dna.chromosome.*.fa", where the asterisk (*) stands for the chromosome.
    The gene list should preferably be ordered by chromosomes.

Bio.Seq and Bio.Alphabet are required to run the script.

Usage: ./motif_based_sampling.py input output

'''
import re
import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

if len(sys.argv)!=3:
    print('Invalid number of arguments.\nUsage: ./motif_based_sampling.py input output')
    sys.exit()

input_file=sys.argv[1]
output_file=sys.argv[2]

motifG=re.compile('[RAG]G[YCT][WAT]') #the RGYW motif
motifC=re.compile('[WAT][RAG]C[YCT]') #the WRCY motif

previous_chr=''
with open(input_file, 'r') as genes, \
open(output_file,'w') as out:
    for i in genes:
        i=i.strip()
        i=re.split('\t',i)
        chr=i[0].replace('chr','')
        start_coord=int(i[1])
        end_coord=int(i[2])

        if previous_chr==chr:
            pass
        if previous_chr!=chr:
            print(f'Chromosome {chr}')
            with open(f'../../mm9_genome/Mus_musculus.NCBIM37.67.dna.chromosome.{chr}.fa', 'r') as chr_fasta:
                chr_list=[]
                for line in chr_fasta:
                    line=line.strip()
                    if not line.startswith('>'):
                        chr_list.append(line)
            chr_list=''.join(chr_list)

        if (i[5]=='+'):
            if (end_coord - start_coord >=500): #if the length of the gene is more than or equal to 500 bp, extract the Cs and Gs in the WRCY and RGYW from the first 500 bp
                piece=chr_list[start_coord-1:start_coord+499].upper()
            else: #if the length of the gene is less than 500 bp, then use the start and end coordinates of the gene to extract the Cs and Gs from
                piece=chr_list[start_coord-1:end_coord].upper()#-1 because indexing starts with 0 in python. So the nth bp will be the nth-1 element in the python object
            for m in motifG.finditer(piece):
                print(chr,start_coord+m.start()+1,'1',m.group()[1], 0, m.group(), i[3], i[7], file=out, sep='\t')
            for m in motifC.finditer(piece):
                print(chr,start_coord+m.start()+2,'1',m.group()[2], 0, m.group(), i[3], i[7], file=out, sep='\t')
        elif (i[5]=='-'): #if we have a minus (-) gene, then the START and END COORDINATES indicate actual END and START OF THE GENE respectively
            if (end_coord - start_coord >=500):
                minus_gene=Seq(chr_list[end_coord-501:end_coord].upper(), IUPAC.ambiguous_dna)
            else:
                minus_gene=Seq(chr_list[start_coord-1:end_coord].upper(), IUPAC.ambiguous_dna)
            piece=str(minus_gene.reverse_complement())
            for m in motifG.finditer(piece):
                print(chr,end_coord-m.start()-1,'-1',m.group()[1], 0, m.group(), i[3], i[7], file=out, sep='\t')
            for m in motifC.finditer(piece):
                print(chr,end_coord-m.start()-2,'-1',m.group()[2], 0, m.group(), i[3], i[7], file=out, sep='\t')


        previous_chr=chr
