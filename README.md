## phylogenomics_Gamma_symbs

2021-07-04
This is just a draft list of scripts and commands that I have used for my quick-and-dirty alignment of Gammaproteobacterial sequences.

To begin, we have used nucleotide and amino acid alignments of 129 protein-coding genes from symbionts listed in https://www.sciencedirect.com/science/article/pii/S0960982219303306
(see supplementary material for details of strains).  
  
We also used lists of nucleotide and amino acid sequences of protein-coding genes extracted from CALKRU/DICMUL/RANSCY Gamma-symbiont genomes.  

```
for gene in *.hmm; do
    SampleName=`basename $gene .hmm`
    echo ------ Searching for gene $SampleName ------
    cat "$SampleName".nucl.al.fa > ../new_alignments/$SampleName.fasta
    for genome in ranscy dicmul calkru; do
        echo searching for $gene in $genome;
        hmmsearch $gene ../prokka_out_$genome/$genome.faa | grep -m 1 -A 2 "E-value" | tail -1 | awk '{print $9}' > ../temp.txt
        samtools faidx ../prokka_out_$genome/$genome.ffn -r ../temp.txt >> ../new_alignments/$SampleName.fasta
    done;
   
    align_nucl_by_codon2.py ../new_alignments/$SampleName.fasta ../new_alignments/"$SampleName"_aligned.fasta
   
done
```

I had to change `align_nucl_by_codon2.py`, because the original alignments included ambiguous IUPAC bases, which were skipped by the original script. Current script version attached!  
  
Then, I substantially rewrote the original `concatenate_seqs.py`. Script `concatenate_seqs_Dictyopharidae.py` attached!  

Working directory: ~/symbio/phylogenomics/new_alignments/

### Tree-building - raxml:
```
raxmlHPC-PTHREADS-SSE3 -f a -T 60 -p 35 -x 345 -#100 -m GTRGAMMA -c 6 -q Parts_codon.txt -s Concatenated.fasta -n Conc_PTHREADS_codons.nwk
```
... using the following Partitions file:
```
DNA, OG0000450_aligned_nucl_1 = 1-147174\3
DNA, OG0000450_aligned_nucl_2 = 2-147174\3
DNA, OG0000450_aligned_nucl_3 = 3-147174\3
```

I then moved on to remove third codon pos from the alignment:
```
for ((i=3;i<=147174;i+=3)); do     echo $i-$i; done > Exclude.txt

raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -s Concatenated.fasta -E Exclude.txt -n Reduce
*You have excluded 49058 out of 147174 columns*
*An alignment file with excluded columns is printed to file Concatenated.fasta.Exclude.txt*
```

Before rerunning the analyses with renamed files:

```
raxmlHPC-PTHREADS-SSE3 -f a -T 60  -p 35 -x 345 -#100  -m GTRGAMMA -c 6 -q Parts2Codon.txt -s Conc2Codon.txt -n Conc2Codon.nwk
```

