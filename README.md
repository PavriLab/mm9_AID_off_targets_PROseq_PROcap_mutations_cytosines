# PROcap, PROseq, mutation frequencies of the mm9 AID off-target genes (Pavri lab)
Steps to produce the PROcap/PROseq/mutations plots for the mm9 AID off-target genes (Pavri Lab, IMP, Vienna, Austria)

### Software used

* Python (v.3.7.4)
  * Biopython’s Bio.SeqIO


* R (v3.6.2) with RStudio (v1.1.463)

* R packages:
  * dplyr
  * ggplot2
  * readODS
  * vioplot
  * ggpubr
  * stringr
  * openxlsx
  * gridExtra
  * ggforce
  * OneR

* Awk (v4.1.4)

* Bash (v4.4.20)

* Ubuntu (v18.04.3 LTS)

* SAMtools (v1.7)
  * merge
  * view
  * reheader
  * index
  * mpileup

* BEDtools (v2.26.0)
  * bamtobed
  * slopBed
  * genomecov

* featureCounts (v2.0.0)

* deepTools2 (v3.3.1)
  * plotCorrelation
  * multiBamSummary
  * bamCoverage

* USCS browser (online tool)
* IGV (v2.5.3)

### Directory structure

```
mkdir 0.External_input_data 1.Replicates_bam 2.Replicates_bedgraph 3.Results 4.All_replicates_to_be_merged 5.Merged_replicates 6.Results_from_merged_enhancers 7.Graphs 8.Results_from_merged_enhancers 9.Graphs_enhancers 10.Graphs_mutations mm9_genome Pre_analysis

cd 1.Replicates_bam
mkdir pro_cap_data pro_seq_data corrected_bg
cd corrected_bg
mkdir up_bg

cd ../../2.Replicates_bedgraph
mkdir Consistency_test

cd ../5.Merged_replicates
mkdir proseq_GCB procap_GCB

cd ../TSS_scripts
```
Then, clone the repository here (in the TSS_scripts directory) and name it TSS_analysis

### BAM files used (from the Pavri lab (GCB mm9 PRO-Seq and PRO-Cap)

The mm9 BAM files of the PRO-Cap and PRO-Seq replicates should be in the directory called 1.Replicates_bam

Samples used:
* **GCB PRO-Seq: 50486, 50487** (named 50486.bam and 50487.bam respectively), in the 1.Replicates_bam/pro_seq subdirectory

* **GCB PRO-cap: 135692** (named gcb_procap.bam), in the 1.Replicates_bam/pro_cap_data subdirectory

### External (mm9) data used
These data are found in the directory called 0.External_input_data. See below how to obtain each one of them.

1. We need the list containing the **mouse (mm9) genes**. The file is called refseq.mm9.2014_0110.imp.withIgTcr.merged.bed. To change the format from bed to gft (necessary step), do:
```
cd path_to_project
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,"REFSEQ\texon",$2+1,$3,".",$6,".\tgene_id \""$4"\";"}' refseq.mm9.2014_0110.imp.withIgTcr.merged.bed > refseq.mm9_updated.merged.gtf
mv refseq.mm9_updated.merged.gtf 0.External_input_data
mv refseq.mm9.2014_0110.imp.withIgTcr.merged.bed 0.External_input_data
```
Now we have our list of mm9 genes in gtf format.

_Note: a similar file called Bcell_pro_seq_expression_levels_cutoff_0_only_tested_percentiles25_75.tpm is used for the analysis. It contains the genes and their coordinates but also a column with gene expression categorization; this column doesn't play a role in the analysis and it is just used to split the final pdf output to multiple pdf files (for easier navigation) according to (arbitrary) groups of genes defined by these categories._

2. We also need a list with **AID off-target genes**. We will use the AID off-target gene list by Álvarez-Prado et al 2018. Using this list we will differentiate between the AID off-target genes and the rest. The file is called "Ramiro_2018 AID targets_JEM_20171738_TableS2.xlsx" and is the original one.


3. For the mutation analysis:
  * the file "PerBase_CG_targets_dKO_PavriLab.ods" was provided to us by Álvarez-Prado and contains positional information about mutation sites and frequencies of mutation per site for the tested genes.

  * we also need the fasta files containing the mm9 genomic sequences. They were downloaded from: ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/ and they are the unmasked versions. Do:
  ```
  wget -m ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/ -P mm9_genome
  ```
  * we need a list of tested genes. We manually edited the "JEM_20171738_TableS1.txt" file changing the names of specific genes with their synomyms (see the "synonyms.txt" file in this folder).

  * we created a list containing the Cs and Gs in WRCY/RGYW hotspots in the first 500bp of each AID off-target gene using a custom python script. For this, do:
```
cd TSS_scripts/TSS_analysis
./motif_based_sampling.py ../../0.External_input_data/Bcell_pro_seq_expression_levels_cutoff_0_only_tested_percentiles25_75.tpm ../../0.External_input_data/wrcy_sites_newpriB_tested_tpm0.txt
cd ../../
```

### Merging the two GCB PROseq replicates
To create the merged replicates, do:
```
samtools merge 1.Replicates_bam/pro_seq/50486_50487_merged.bam 1.Replicates_bam/pro_seq/50486.bam 1.Replicates_bam/pro_seq/50487.bam
rm 1.Replicates_bam/pro_seq/50486.bam 1.Replicates_bam/pro_seq/50487.bam
```

### From PRO-Cap/Seq BAM replicates to coverage (BEDGRAPH) files

Go to the 1.Replicates_bam directory and get the chromosome sizes of the mm9 assembly (the file should be named mm9.genome)
```
cd 1.Replicates_bam
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
"select chrom, size from mm9.chromInfo" > mm9.genome
```
Also, place a copy of the split_and_correct_bg.sh script in this directory:
```
cp ../TSS_scripts/TSS_analysis/split_and_correct_bg.sh .
```

Now we are ready to proceed with the generation of the BEDGRAPH files. The split_and_correct_bg.sh script
splits the strands in minus and plus and calculates the coverage using only 3' ends for PRO-Seq or
only 5' ends for PRO-Cap using bedtool's genomecov, and corrects the "multisignals" (PROCap and PROSeq offer single nucleotide resolution, however, consecutive bases with exactly the same coverage get merged and we need to correct for this) using the
custom sign_correction.py script. We will then normalize by library size using awk.

* For the mm9 GCB PRO-Cap
```
mv pro_cap_data/gcb_procap.bam .
./split_and_correct_bg.sh procap mm9.genome
#correcting some incompatibility issues and removing non-mm9 reads
cd corrected_bg
for file in *.bedgraph; do sudo cat $file | sed 's/^mm9_/chr/g' | grep -v dm5 | grep -v NT_ > upg_"$file"; done
mv upg* up_bg
#normalizing
cd up_bg
for file in *minus.bedgraph;
do sum_m=$(awk '{sum+=$4;} END{print sum;}' "$file");
v=$(echo $file | sed 's/minus/plus/');
sum_p=$(awk '{sum+=$4;} END{print sum;}' "$v");
sum_total=$((sum_m+sum_p));
awk '{print ($1, $2, $3, -1000000*$4/'$sum_total')}' 'OFS=\t' "$file" > normalized_"$file";
sed -i 's/chr//g' normalized_"$file";
awk '{print ($1, $2, $3, 1000000*$4/'$sum_total')}' 'OFS=\t' "$v" > normalized_"$v";
sed -i 's/chr//g' normalized_"$v";
sum_p=0;
sum_m=0;
sum_total=0;
done
mv normalized* ../../../5.Merged_replicates/procap_GCB
cd ../../
mv gcb_procap.bam pro_cap_data
```

  and "clean" the corrected_bg, up_bg, and 1.Replicates_bam directories from any by-product files created.
  ```
  rm corrected_bg/* up_bg/*
  ```

* Now for the mm9 GCB PRO-Seq
```
mv pro_seq_data/50486_50487_merged.bam .
./split_and_correct_bg.sh proseq mm9.genome
#correcting some incompatibility issues and removing non-mm9 reads
cd corrected_bg
for file in *.bedgraph; do sudo cat $file | sed 's/^mm9_/chr/g' | grep -v dm5 | grep -v NT_ > upg_"$file"; done
mv upg* up_bg
#normalizing
cd up_bg
for file in *minus.bedgraph;
do sum_m=$(awk '{sum+=$4;} END{print sum;}' "$file");
v=$(echo $file | sed 's/minus/plus/');
sum_p=$(awk '{sum+=$4;} END{print sum;}' "$v");
sum_total=$((sum_m+sum_p));
awk '{print ($1, $2, $3, -1000000*$4/'$sum_total')}' 'OFS=\t' "$file" > normalized_"$file";
sed -i 's/chr//g' normalized_"$file";
awk '{print ($1, $2, $3, 1000000*$4/'$sum_total')}' 'OFS=\t' "$v" > normalized_"$v";
sed -i 's/chr//g' normalized_"$v";
sum_p=0;
sum_m=0;
sum_total=0;
done
mv normalized* ../../../5.Merged_replicates/
cd ../../
mv 50486_50487_merged.bam pro_seq_data
```

  and "clean" the corrected_bg, up_bg, and 1.Replicates_bam folder from any by-product files created.
  ```
  rm corrected_bg/* up_bg/*
  cd ..
  ```


Now the folder 5.Merged_replicates is populated with the normalized BEDGRAPH files GCB PROCap and PROSeq.

### mm9 PROCap, PROSeq, and mutation tracks for the AID off-targets

To generate the region-specific plots of the mutated genes, do:
```
R
source("TSS_scripts/TSS_analysis/mutations_distance_procap_proseq.R")
```

The output will be placed in the 10.Graphs_mutations directory.

### mm9 PROSeq (GCB) C:Gs binned by PROseq RPM and mutation frequency
```
R
source("TSS_scripts/TSS_analysis/mutations_analysis_all_categories.R")
```
