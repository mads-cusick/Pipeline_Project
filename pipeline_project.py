# Loyola University Chicago
# Comp 383 - Pipeline Project
# Mads Cusick

# import necessary packages
import os
import Bio
from Bio import Entrez
from Bio import SeqIO
import csv

# make a directory to run everything in
os.system('mkdir PipelineProject_Mads_Cusick')
# enter into the direstory
os.chdir('PipelineProject_Mads_Cusick')

# open the log file to write to
log = open('PipelineProject.log', 'w')

# make a list of sample SRR names and IDs from two patient donors
sample_srr = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
sample_id = ['Donor 1 (2dpi)', 'Donor 1 (6dpi)', 'Donor 3 (2dpi)', 'Donor 3 (6dpi)']

# retrieve the transcriptomes from the two patient donors from SRA using wget
for sample in sample_srr:
    os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/'+sample+'/'+sample)

# convert to paired-end fastq files with fastq-dump
for sample in sample_srr:
    os.system('fastq-dump -I --split-files '+sample)



# make a small subset (10000 reads) from the paired-end fastq files
for sample in sample_srr:
    for i in range(1,3):
        os.system('head -10000 '+sample+'_'+i+'.fastq > test_'+sample+'_'+i+'.fastq')



# obtain HCMV genome fasta from NCBI and save as fasta file
Entrez.email = "mcusick2@luc.edu"
# use efetch to retrieve the HCMV record in the fasta format
handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype="fasta")
# parse the handle into a FastaIterator object
seq_record_iterator = SeqIO.parse(handle, "fasta")
# write the sequence from the FastaIterator object out to a fasta file
SeqIO.write(seq_record_iterator, "NC_006273.2_genome.fasta", "fasta")

# create index for HCMV with Bowtie2
os.system('bowtie2-build NC_006273.2_genome.fasta HCMV')

# map to index with Bowtie2 for each of the paired-end reads
for sample in sample_srr:
    os.system('bowtie2 -x HCMV -1 '+sample+'_1.fastq -2 '+sample+'_2.fastq -S HCMVmap_'+sample+'_.sam --al-conc-gz '+sample+'_mapped_%.fq.gz')

# write to log file the number of reads in each transcriptome before and after the Bowtie2 mapping
for i in range(len(sample_srr)):
    # count the lines of the before and after fastq files and divide by 4 to get the number of reads
    before = os.system('wc '+sample_srr[i]+'_1.fastq') / 4
    after = os.system('wc '+sample_srr[i]+'_mapped_1.fq.gz') / 4
    # write to log file
    log.write(sample_id[i]+' had '+str(before)+' read pairs before Bowtie2 filtering and '+str(after)+' read pairs after.\n')
log.write('\n')



# assemble all 4 transcriptomes together to produce 1 assembly via SPAdes
os.system('spades.py -k 77,99,127 -t 2 --only-assembler '+
          '--pe-1 1 SRR5660030_mapped_1.fq.gz --pe-2 1 SRR5660030_mapped_2.fq.gz '+
          '--pe-1 2 SRR5660033_mapped_1.fq.gz --pe-2 2 SRR5660033_mapped_2.fq.gz '+
          '--pe-1 3 SRR5660044_mapped_1.fq.gz --pe-2 3 SRR5660044_mapped_2.fq.gz '+
          '--pe-1 4 SRR5660045_mapped_1.fq.gz --pe-2 4 SRR5660045_mapped_2.fq.gz '+
          '-o SRR56600_assembly/')

# write SPAdes command to log file
log.write('spades.py -k 77,99,127 -t 2 --only-assembler '+
          '--pe-1 1 SRR5660030_mapped_1.fq.gz --pe-2 1 SRR5660030_mapped_2.fq.gz '+
          '--pe-1 2 SRR5660033_mapped_1.fq.gz --pe-2 2 SRR5660033_mapped_2.fq.gz '+
          '--pe-1 3 SRR5660044_mapped_1.fq.gz --pe-2 3 SRR5660044_mapped_2.fq.gz '+
          '--pe-1 4 SRR5660045_mapped_1.fq.gz --pe-2 4 SRR5660045_mapped_2.fq.gz '+
          '-o SRR56600_assembly/\n')
log.write('\n')



# copy the scaffolds.fasta file into the current working directory
os.system('cp SRR56600_assembly/scaffolds.fasta scaffolds.fasta')

# read in the scaffolds.fasta file
handle = open('scaffolds.fasta')
# parse the scaffolds.fasta file with SeqIO and create a list of SeqRecord objects
seqIOlist = list(SeqIO.parse(handle, 'fasta'))

# extract the id and sequence of each contig from the SeqRecord objects and store in a dictionary
seq_dict = {str(contig.id):str(contig.seq) for contig in seqIOlist}

# store the id and sequences of only the contigs with lengths > 1000 in a dictionary (from the dictionary of all the contigs)
long_contigs_only = {contig_id:contig_seq for contig_id, contig_seq in seq_dict.items() if len(contig_seq) > 1000}

# write to log file the number of contigs > 1000 bp in the assembly (from the length of the long_contigs_only dictionary)
log.write('There are ' + str(len(long_contigs_only)) + ' contigs > 1000 bp in the assembly\n')

# calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length)
assembly_length = sum(len(contig_seq) for contig_seq in long_contigs_only.values())

# write to log file the total length of the assembly
log.write('There are ' + str(assembly_length) + ' bp in the assembly\n')
log.write('\n')



# retrieve the longest contig from the SPAdes assembly
# loop through all of the contigs in the long_contigs_only dictionary and only store the contig if it has the longest length of all the contigs
longest_contig = {contig_id:contig_seq for contig_id, contig_seq in long_contigs_only.items() if contig_seq == max(long_contigs_only.values())}

# write the longest contig to a fasta file
with open('longest_contig.fasta', 'w') as f:
    # write the first line with > and the id
    f.write('>' + list(longest_contig.keys())[0])
    # write the sequence on the next line
    f.write('\n' + list(longest_contig.values())[0])



# obtain Betaherpesvirinae subfamily sequences from NCBI and save as fasta file
Entrez.email = "mcusick2@luc.edu"
# use esearch to retrieve all the genbank IDs of the betaherpesvirinae subfamily
handle = Entrez.esearch(db="nucleotide", term='Betaherpesvirinae[Organism]', retmax=20000)
# use Entrez to read the handle and extract the IDs and save them to a list
id_list = Entrez.read(handle)['IdList']
# use efetch to retrieve the sequences from the list of IDs of the betaherpesvirinae subfamily
handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta")
# parse the handle into a FastaIterator object
seq_record_iterator = SeqIO.parse(handle, "fasta")
# write the sequences from the FastaIterator object out to a fasta file
SeqIO.write(seq_record_iterator, "Betaherpesvirinae_sequences.fasta", "fasta")



# use donwloaded Betaherpesvirinae sequences to make a local database
os.system('makeblastdb -in Betaherpesvirinae_sequences.fasta -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl')

# Use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily
# only keep the best alignment for any single query-subject pair of sequences
# write all the hits to a tab-delimited file with the specified info
os.system('blastn -query longest_contig.fasta -db betaherpesvirinae -out betaherpesvirinae_blastn_results.csv -max_hsps 1 '+
          '-outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"')

# write to log file the top 10 hits
# make a list of names for the header row and write them to the log file separated by tabs
header_row = ['sacc','pident','length','qstart','qend','sstart','send','bitscore','evalue','stitle']
log.write('\t'.join(header_row))
log.write('\n')

# iterate through the blastn results file and write to the log file
with open('betaherpesvirinae_blastn_results.csv', 'r') as file:
    # create a list of the first 10 lines
    head = [next(file) for i in range(10)]
    # write only these first 10 to the log file
    for line in head:
        log.write(line)

# close the log file
log.close()







