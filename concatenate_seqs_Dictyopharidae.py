#! /usr/bin/env python3

### This script takes a list of multifasta files containing alignments, where individual sequences are labelled SAMPLE_NAME + "_" + GENE,
### and produces a multifasta where all gene sequences for a given sample are concatenated.
### It assumes that sequence names follow one of the below formats:
### >OG0000450.GCF_002291445.1_ASM229144v1.locus_tag=AWC36_RS17565
### >DICMUL_00057



### Usage: run in the directory containing files to concatenate

import glob, sys, os, re

if len(sys.argv) != 2:
    sys.exit('This script concatenates all fasta sequences within a directory\n'
	         'Sequence names as in 20210701 alignments! The script outputs concatenation,\n'
	         ' as well as partition files, in "./Concatenated" directory in the working directory.\n'
	         'Usage: ./concatenate_seqs.py <directory> \n')

Script, work_dir = sys.argv


def ImportFastaAsList(fasta_file):
    FASTA = open(fasta_file, 'r')
    Seq_list = []
    Sequence = ''
    Seq_heading = ''
    for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
        line = line.strip()
        if line.startswith('>'):
            if Sequence != '':   # Saves the existing Seq_heading and Sequence to a list before overwriting them
                Seq_list.append([Seq_heading, Sequence])
                Sequence = ''
                
            if line.startswith(">OG0"):
                Seq_heading = re.sub(r'.*\.(GC.*)\.locus.*', r'\1', line)
            elif re.search(r'>\w{6}_\d+', line):
                Seq_heading = line[1:7]
            else:
                sys.exit(f"ERROR!"
                         f"Unexpected sequence name - {line.strip()} in file {fasta_file}"
                         "Exiting.....")

        else:
            Sequence += line.strip().upper()
            
    Seq_list.append([Seq_heading, Sequence]) # Saves the final sequence (Seq_heading and Sequence) to a list
    
    ### Check gene lengths... should all be the same!
    ali_len = 0
    for seq in Seq_list:
        if ali_len != 0 and ali_len != len(seq[1]):
            sys.exit(f"ERROR!\n"
                     f"Not all genes in alignment {fasta_file} have the same length!\n"
                     f"Alignment length = {ali_len}, but the length of {seq[0]} = {len(seq[1])}")
        ali_len = len(seq[1])    
    
    print("Gene successfully added from file %s: %d sequences of %d nucleotides" % (fasta_file, len(Seq_list), ali_len))

    FASTA.close()
    return(Seq_list)



# 0) Print general info

print("This script takes a list of multifasta files containing alignments and produces a multifasta where all gene sequences for a given sample are concatenated.\n")
print("Working directory: ", end = "")
os.system("pwd")


# 1) List genes to be concatenated
print("Step 1) Listing fasta files containing alignments that should be concatenated: ")
gene_list = []
for file in glob.glob("%s/*.fasta" % work_dir):
    gene_list.append(file.split("/")[-1].split(".")[0])
    print(" - ", gene_list[-1])
print("\n")
if gene_list == []:
    sys.exit("No fasta files found. Are you pointing at the directory where files are located? Aborting....")


# 1a) make table where samples are represented by rows and genes by columns
# First, making header row, where gene names are listed
all_genes_dict = {'Sample':[]}
for gene in gene_list:
    all_genes_dict['Sample'].append(gene)

# 2) read fasta files, fill gene_table:
# [['Sample', 'gene1', 'gene2'],
#  ['RANSCY', 'ACTG', 'AGTATG']]

print("\n##########\n###### Now, constructing gene table - and adding contents of individual fasta files\n##########\n")

gene_lens = []

for gene in gene_list:
    seq_list = ImportFastaAsList(gene + ".fasta")
    gene_lens.append(len(seq_list[0][1]))

    
    #### Populating the first column of gene table: [['Sample', ...], ['TETUND1', ...], ['TETULN', ....]]
    if len(all_genes_dict) == 1:
        for seq in seq_list:
            all_genes_dict[seq[0]] = [] ### for each strain, creating a separate entry in the dictionary
    
    #### Adding sequences of all genes to all_genes_dict
    print("   Genomes where the gene was found: ", end = "")
    genomes_where_gene_found = []
    
    for seq in seq_list:                     ###  for each sequence (two-element list: ['Seq1', 'ACTGTA'])
        seq_name, sequence = seq
        if seq_name in all_genes_dict.keys():          ### if strain name in dictionary:
            all_genes_dict[seq_name].append(sequence)  ###    just add the sequence to list for that strain
            genomes_where_gene_found.append(seq_name)  ###    and add genome name to the list
            #print(f"-{seq_name} in all_genes_dict.keys()")
            
        else:                                ### otherwise:
            all_genes_dict[seq_name] = []    ###     create new dictionary item for that strain
            for gene_len in gene_lens[:-1]:  ###     and fill up previous gene info with dashes!
                all_genes_dict[seq_name].append("-"*gene_len)
            all_genes_dict[seq_name].append(sequence)  ###    and add latest sequence to list
            genomes_where_gene_found.append(seq_name)  ###    and add genome name to the list
            #print(f"---- {seq_name} NOT in all_genes_dict.keys(), added to all_genes_dict")
        
        print(seq_name, end = ", ")
    print("\n")

    print("   Genomes where the gene not found: ", end = "")    
    for genome_name in all_genes_dict:
        if genome_name != 'Sample' and genome_name not in genomes_where_gene_found:
            all_genes_dict[genome_name].append("-"*gene_lens[-1])
            print(genome_name, end = ", ")
    print("\n")
    #print("   All genomes in the dict", all_genes_dict.keys(), "\n")
    
    print(f"--------------------------------------------------------\n"
          f"---- end of gene {gene} ----\n\n")


##### Printing some more summaries....   
print("Genomes in gene_table: ")
for genome_name in all_genes_dict.keys():
    if not genome_name == "Sample":
        print(f"{genome_name}", end = ", ")
print(f"\n ----- Total genomes: {len(all_genes_dict.keys())-1}  ----- \n")


print("Genes in gene_table: ")
for gene_name in all_genes_dict["Sample"]:
    print(f"{gene_name}", end = ", ")
print(f"\n ----- Total genes: {len(all_genes_dict['Sample'])}  ----- \n")



# 3) export sequence concatenation
print("\n\n##########\n###### Now, exporting sequence concatenation....\n##########\n")

print("Now exporting gene concatenation........", end = "")
if not glob.glob("%s/Concatenated/" % work_dir):
    os.system("mkdir %s/Concatenated" % work_dir)

CONC = open("%s/Concatenated/Concatenated.fasta" % work_dir, 'w')
for genome_name in all_genes_dict.keys():
    if not genome_name == "Sample":
        print(f">{genome_name}", file = CONC)
        for gene in all_genes_dict[genome_name]:
            print(gene, end = "", file = CONC)
        print("", file = CONC)

CONC.close()

print("DONE!\nConcatenation successfully exported as './Concatenated/Concatenated.fasta'\n")




# 4) export DNA partitions file

print("\n\n##########\n###### Now exporting partitions file as './Concatenated/Partitions.txt'... Contents as below.\n##########\n")

PARTITIONS_FILE = open("%s/Concatenated/Partitions.txt" % work_dir, "w")
pos_ct = 1

for i in range(len(gene_list)):
    print("DNA, %s=%d-%d" % (gene_list[i], pos_ct, pos_ct+gene_lens[i]-1), file=PARTITIONS_FILE)
    print("DNA, %s=%d-%d" % (gene_list[i], pos_ct, pos_ct+gene_lens[i]-1))
    pos_ct += gene_lens[i]

PARTITIONS_FILE.close()
    
    
# 5) export PartitionFinder partitions file
print("\nNow exporting PartitionFinder partitions file as './Concatenated/PartitionFinder.txt'... Contents as below.")

PARTITIONS_FILE = open("%s/Concatenated/PartitionFinder.txt" % work_dir, "w")
pos_ct = 1

for i in range(len(gene_list)):
    print("%s_1=%d-%d\\3;" % (gene_list[i], pos_ct, pos_ct+gene_lens[i]-1), file=PARTITIONS_FILE)
    print("%s_2=%d-%d\\3;" % (gene_list[i], pos_ct+1, pos_ct+gene_lens[i]-1), file=PARTITIONS_FILE)
    print("%s_3=%d-%d\\3;" % (gene_list[i], pos_ct+2, pos_ct+gene_lens[i]-1), file=PARTITIONS_FILE)
    print("%s_1=%d-%d\\3;" % (gene_list[i], pos_ct, pos_ct+gene_lens[i]-1))
    print("%s_2=%d-%d\\3;" % (gene_list[i], pos_ct+1, pos_ct+gene_lens[i]-1))
    print("%s_3=%d-%d\\3;" % (gene_list[i], pos_ct+2, pos_ct+gene_lens[i]-1))
    pos_ct += gene_lens[i]

PARTITIONS_FILE.close()

print("\nDONE! Script executed successfully! :D\n")

    
    
    
    
#### Now, all I need to do is run RAxML - using one of the following commands
# raxmlHPC -m GTRGAMMA -p 12345 -s Concatenated.fasta -n Conc
# raxmlHPC -m GTRGAMMA -p 12345 -x 12345 -# 100 -s ileS.fasta -n ileS_boot.tre
# raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s ileS.fasta -n ileS3.tre
