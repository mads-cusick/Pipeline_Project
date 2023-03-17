# Pipeline_Project
Python wrapper to automate the execution of bioinformatics software tools

## Description
The goal of this project is to create a python wrapper that will automate the execution of various bioinformatics tools used for genome assembly. The transcriptomes of Human herpesvirus 5 (HCMV) from two patient donors 2- and 6-days post-infection were analyzed. The strains that are most similar to these patient samples were determined by using Bowtie2 to create an index for HCMV and saving the reads that mapped to this index for use in assembly. Then all four transcriptomes were mapped together to produce the single assembly using SPAdes. This assembly was aligned with other virus strains using NCBI’s BLAST algorithm via blast+, using the longest contig from the assembly to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily.

## Dependencies
* macOS or similar Unix-based OS
* Bowtie2 installed
* SPAdes installed
* blast+ installed
* Biopython installed

## Installing
Download the pipeline_project.py file and save to your working directory of your terminal

## Executing program
Run the program from your command line by calling it with python pipeline_project.py
