# The script will take as input the per base C and G mm9 mutations, the list of tested genes
# with categorization of the genes into low, moderate and highly mutated, and a list with
# the positions of wrcy motifs (was produced in an earlier step). It will
# also ask for the PROCap & PROSeq replicates (plus and minus strands
# in BEDGRAPH format).

# It will then create a series of plots for the [-100,+550] region of each gene combining the signals
# from the aforementioned replicates (and the mutations for the [0,+500] region) and positioning them accordingly to the sense
# and antisense strand and the annotated start of the gene. On the tracks produced, the WRCY motifs will be highlighted.

# The script will also check whether the max distance between two mutation sites in all genes
# is less than around 600bp and remove from the analysis genes with all of their mutations more than 500 bp away from the annotated gene start.

# The script uses the custom-made Internal.transcription.sites() function.

if (!require(readODS)) install.packages('readODS', repos = 'http://cran.us.r-project.org')
library('readODS')

if (!require(ggplot2)) install.packages('ggplot2', repos = 'http://cran.us.r-project.org')
library('ggplot2')

if (!require(dplyr)) install.packages('dplyr', repos = 'http://cran.us.r-project.org')
library('dplyr')

if (!require(vioplot)) install.packages('vioplot', repos = 'http://cran.us.r-project.org')
library('vioplot')

#args = commandArgs(trailingOnly=TRUE)

categorizing=1 #as.integer(readline('Do you want the genes split into different categories? If yes, enter 1: ')), used to give the user a choice, but I removed it since we always use categorizing now

#importing the sheets
#distance_mutations<-read.ods(file.choose()) # the file with the per-base mm9 mutations, PerBase_CG_targets_dKO_PavriLab.ods
distance_mutations<-read.ods("0.External_input_data/PerBase_CG_targets_dKO_PavriLab.ods")
#distance_mutations<-read_ods(args[1])

#formatting the data frames
for (i in 1:length(distance_mutations)){
  distance_mutations[[i]]<-distance_mutations[[i]][c(-5,-6,-7,-8,-12,-13)]
  distance_mutations[[i]]<-distance_mutations[[i]][-1,]
  colnames(distance_mutations[[i]])<-c('chr','end','strand','base','motif','gene','freqA','freqT')
}

#combining the data frames
conc_mutations<-distance_mutations[[1]]
for (i in 2:length(distance_mutations)){
  conc_mutations<-rbind(conc_mutations,distance_mutations[[i]])
}
conc_mutations_C<-conc_mutations[which(conc_mutations$base=='C' & conc_mutations$freqT>0),-7]
conc_mutations_G<-conc_mutations[which(conc_mutations$base=='G' & conc_mutations$freqA>0),-8]
colnames(conc_mutations_C)<-c('chr','end','strand','base','motif','gene','freq')
colnames(conc_mutations_G)<-colnames(conc_mutations_C)

conc_mutations_C$freq<-as.numeric(conc_mutations_C$freq)
conc_mutations_G$freq<-as.numeric(conc_mutations_G$freq)
conc_mutations_G$freq<-log10(conc_mutations_G$freq+1)*(-1)
conc_mutations_C$freq<-log10(conc_mutations_C$freq+1)

conc_mutations<-rbind(conc_mutations_C,conc_mutations_G)

conc_mutations<-conc_mutations[,c('chr','end','strand','base','freq','motif','gene')]

conc_mutations$gene[]<-lapply(conc_mutations$gene, gsub, pattern='_[A-Za-z0-9]+',replacement='')
conc_mutations$gene[]<-lapply(conc_mutations$gene, gsub, pattern='-[0-9]+',replacement='')
conc_mutations$strand[which(conc_mutations$strand=='1')]<-'+'
conc_mutations$strand[which(conc_mutations$strand=='-1')]<-'-'

distance_mutations<-conc_mutations[c(1,2,3,5,7)]

distance_mutations$gene[which(distance_mutations$gene=='Rad9')]<-list('Rad9a')
distance_mutations$gene[which(distance_mutations$gene=='H47')]<-list('Vimp')

sorted_df<-dplyr::arrange(distance_mutations,chr,end)
sorted_df$end<-as.integer(sorted_df$end)
genes<-unique(sorted_df$gene)

#importing the refseq file
refseq<-read.csv("0.External_input_data/Bcell_pro_seq_expression_levels_cutoff_0_only_tested_percentiles25_75.tpm", sep='\t', header=FALSE) #IMPORTANT: if you chose to categorize, be sure to provide a refseq file with categorization of the genes
#refseq<-read.csv(args[2], sep='\t', header=FALSE) #IMPORTANT: if you chose to categorize, be sure to provide a refseq file with categorization of the genes

colnames(refseq)<-c('chr','Start','End','Geneid','zero','strand','tpm','exp_lvl')
refseq$Geneid<-as.character(refseq$Geneid)
refseq$Geneid[]<-lapply(refseq$Geneid, gsub, pattern='_[A-Za-z0-9]+',replacement='')
refseq$Geneid[]<-lapply(refseq$Geneid, gsub, pattern='-[0-9]+',replacement='')


#the results' data frame
results_df<-as.data.frame(matrix(NA,nrow=length(genes),ncol=13))
colnames(results_df)<-c('gene','chr','first_position','last_position','mean_distance_between_mutations','difference_between_first_and_last','distance_from_start','all_positions','all_positions_from_start','gene_start','gene_end','strand','frequencies')
results_df$all_positions<-as.list(results_df$all_positions)
results_df$all_positions_from_start<-as.list(results_df$all_positions_from_start)
results_df$frequencies<-as.list(results_df$frequencies)
j=1

for (i in genes){
  results_df$gene[j]<-i
  results_df$chr[j]<-sorted_df[which(sorted_df$gene==i),]$chr[1]
  results_df$last_position[j]<-tail(sorted_df[which(sorted_df$gene==i),][[2]], n=1)
  results_df$first_position[j]<-sorted_df[which(sorted_df$gene==i),][[2]][1]
  results_df$mean_distance_between_mutations[j]<-mean(diff(sorted_df$end[which(sorted_df$gene==i)]))
  results_df$difference_between_first_and_last[j]<-tail(sorted_df[which(sorted_df$gene==i),][[2]], n=1) - sorted_df[which(sorted_df$gene==i),][[2]][1]
  results_df$all_positions[[j]]<-sorted_df$end[which(sorted_df$gene==i)]
  results_df$frequencies[[j]]<-sorted_df$freq[which(sorted_df$gene==i)]
  if (i %in% refseq$Geneid){
    results_df$strand[j]<-sorted_df$strand[which(sorted_df$gene==i)][1]
    results_df$gene_start[j]<-refseq$Start[which(refseq$Geneid==i)]
    results_df$gene_end[j]<-refseq$End[which(refseq$Geneid==i)]
    if (sorted_df$strand[which(sorted_df$gene==i)][1]=='-'){
      results_df$distance_from_start[j]<-refseq$End[which(refseq$Geneid==i)] - sorted_df[which(sorted_df$gene==i),][[2]][1]
      results_df$all_positions_from_start[[j]]<-refseq$End[which(refseq$Geneid==i)] - results_df$all_positions[[j]]
    } else {
      results_df$distance_from_start[j]<-tail(sorted_df[which(sorted_df$gene==i),][[2]], n=1) - refseq$Start[which(refseq$Geneid==i)]
      results_df$all_positions_from_start[[j]]<- results_df$all_positions[[j]] - refseq$Start[which(refseq$Geneid==i)]
    }
  }
  j=j+1
}

results_df<-na.omit(results_df)
results_df_500bp<-results_df[which(results_df$distance_from_start<501 & results_df$distance_from_start>(-100)),]
row.names(results_df_500bp)<-NULL

#extracting the TSSs
test<-read.delim("5.Merged_replicates/procap_GCB/normalized_upg_corrected_gcb_procap_minus.bedgraph", header=FALSE) #load the PROCap minus strand
#test<-read.delim(args[3], header=FALSE) #load the PROCap minus strand

results_minus<-matrix(NA, nrow=length(results_df_500bp$gene),ncol=6)
results_minus<-as.data.frame(results_minus)
intermediate<-as.data.frame(matrix(0, nrow=1, ncol=5))
source("TSS_scripts/TSS_analysis/Internal.transcription.sites.R") #sourcing the Internal.transcription.sites.R function
for (i in 1:length(results_df_500bp$gene)){
  if (results_df_500bp$strand[i]=="+"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_start[i])-100),(as.integer(results_df_500bp$gene_start[i])+550)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_minus$V1[i]<-intermediate[1]
    results_minus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_minus$V3[i]<-intermediate[3]
    results_minus$V4[i]<-intermediate[4]
    results_minus$V5[i]<-intermediate[5]
    results_minus$V6[i]<-results_df_500bp$gene[i]
    results_minus$V7[i]<-as.integer(results_df_500bp$gene_start[i])
    results_minus$V8[i]<-results_df_500bp$strand[i]
  } else if (results_df_500bp$strand[i]=="-"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_end[i])-550),(as.integer(results_df_500bp$gene_end[i])+100)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_minus$V1[i]<-intermediate[1]
    results_minus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_minus$V3[i]<-intermediate[3]
    results_minus$V4[i]<-intermediate[4]
    results_minus$V5[i]<-intermediate[5]
    results_minus$V6[i]<-results_df_500bp$gene[i]
    results_minus$V7[i]<-as.integer(results_df_500bp$gene_end[i])
    results_minus$V8[i]<-results_df_500bp$strand[i]
  }
}
row.names(results_minus)<-NULL

test<-read.delim("5.Merged_replicates/procap_GCB/normalized_upg_corrected_gcb_procap_plus.bedgraph", header=FALSE) #load the PROCap plus strand
#test<-read.delim(args[4], header=FALSE) #load the PROCap plus strand

results_plus<-matrix(NA, nrow=length(results_df_500bp$gene),ncol=8)
results_plus<-as.data.frame(results_plus)
intermediate<-as.data.frame(matrix(0, nrow=1, ncol=5))
for (i in 1:length(results_df_500bp$gene)){
  if (results_df_500bp$strand[i]=="+"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_start[i])-100),(as.integer(results_df_500bp$gene_start[i])+550)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_plus$V1[i]<-intermediate[1]
    results_plus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_plus$V3[i]<-intermediate[3]
    results_plus$V4[i]<-intermediate[4]
    results_plus$V5[i]<-intermediate[5]
    results_plus$V6[i]<-results_df_500bp$gene[i]
    results_plus$V7[i]<-as.integer(results_df_500bp$gene_start[i])
    results_plus$V8[i]<-results_df_500bp$strand[i]
  } else if (results_df_500bp$strand[i]=="-"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_end[i])-550),(as.integer(results_df_500bp$gene_end[i])+100)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_plus$V1[i]<-intermediate[1]
    results_plus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_plus$V3[i]<-intermediate[3]
    results_plus$V4[i]<-intermediate[4]
    results_plus$V5[i]<-intermediate[5]
    results_plus$V6[i]<-results_df_500bp$gene[i]
    results_plus$V7[i]<-as.integer(results_df_500bp$gene_end[i])
    results_plus$V8[i]<-results_df_500bp$strand[i]
  }
}
row.names(results_plus)<-NULL

colnames(results_minus)<-c('hits','chr','starting_position','ending_position','signal_strength','gene','annotated_start','strand')
colnames(results_plus)<-c('hits','chr','starting_position','ending_position','signal_strength','gene','annotated_start','strand')

results_plus$distance_from_start<-list(NA)
for (i in 1:length(results_plus$hits)){
  if (results_plus$strand[i]=='+'){
    results_plus$distance_from_start[[i]]<-(as.numeric(results_plus$ending_position[[i]])-as.numeric(results_plus$annotated_start[i]))
  } else {
    results_plus$distance_from_start[[i]]<-as.numeric(results_plus$annotated_start[i])-(as.numeric(results_plus$ending_position[[i]]))
  }
}
results_minus$distance_from_start<-list(NA)
for (i in 1:length(results_minus$hits)){
  if (results_minus$strand[i]=='+'){
    results_minus$distance_from_start[[i]]<-(as.numeric(results_minus$ending_position[[i]])-as.numeric(results_minus$annotated_start[i]))
  } else {
    results_minus$distance_from_start[[i]]<-as.numeric(results_minus$annotated_start[i])-(as.numeric(results_minus$ending_position[[i]]))
  }
}

results_sense_minus<-results_minus[which(results_minus$strand=="-"),]
results_antisense_for_plus_genes<-results_minus[which(results_minus$strand=="+"),]

results_sense_plus<-results_plus[which(results_plus$strand=="+"),]
results_antisense_for_minus_genes<-results_plus[which(results_plus$strand=="-"),]

results_sense<-rbind(results_sense_minus,results_sense_plus)
results_sense$index<-as.numeric(rownames(results_sense)) #getting the index
results_sense<-results_sense[order(results_sense$index),] #ordering according to the index

results_antisense<-rbind(results_antisense_for_minus_genes,results_antisense_for_plus_genes)
results_antisense$index<-as.numeric(rownames(results_antisense)) #getting the index
results_antisense<-results_antisense[order(results_antisense$index),] #ordering according to the index

results_antisense$signal_strength<-lapply(results_antisense$signal_strength,as.numeric)
results_sense$signal_strength<-lapply(results_sense$signal_strength,as.numeric)

results_antisense$distance_from_start[which(is.na(results_antisense$distance_from_start))]<-0
results_sense$distance_from_start[which(is.na(results_sense$distance_from_start))]<-0

results_df_500bp$frequencies<-lapply(results_df_500bp$frequencies,as.numeric)
results_df_500bp$sum_of_freq<-lapply(results_df_500bp$frequencies,sum)
results_df_500bp$sum_of_freq<-as.numeric(results_df_500bp$sum_of_freq)

results_df_500bp$sum_strength_sense<-lapply(results_sense$signal_strength,sum)
results_df_500bp$sum_strength_sense<-lapply(results_df_500bp$sum_strength_sense, abs)

results_df_500bp$sum_strength_antisense<-lapply(results_antisense$signal_strength,sum)
results_df_500bp$sum_strength_antisense<-lapply(results_df_500bp$sum_strength_antisense, abs)


results_df_500bp$freq_hits<-lapply(results_df_500bp$all_positions,length)

results_df_500bp$sense_hits<-results_sense$hits
results_df_500bp$antisense_hits<-results_antisense$hits

results_df_500bp$freq_hits<-as.numeric(results_df_500bp$freq_hits)
results_df_500bp$sense_hits<-as.numeric(results_df_500bp$sense_hits)
results_df_500bp$antisense_hits<-as.numeric(results_df_500bp$antisense_hits)

#extracting the pauses
test<-read.delim("5.Merged_replicates/proseq_GCB/normalized_upg_corrected_50486_50487_merged_minus.bedgraph", header=FALSE) #load the PROSeq minus strand
#test<-read.delim(args[5], header=FALSE) #load the PROSeq minus strand

results_proseq_minus<-matrix(NA, nrow=length(results_df_500bp$gene),ncol=6)
results_proseq_minus<-as.data.frame(results_proseq_minus)
intermediate<-as.data.frame(matrix(0, nrow=1, ncol=5))
for (i in 1:length(results_df_500bp$gene)){
  if (results_df_500bp$strand[i]=="+"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_start[i])-100),(as.integer(results_df_500bp$gene_start[i])+550)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_proseq_minus$V1[i]<-intermediate[1]
    results_proseq_minus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_proseq_minus$V3[i]<-intermediate[3]
    results_proseq_minus$V4[i]<-intermediate[4]
    results_proseq_minus$V5[i]<-intermediate[5]
    results_proseq_minus$V6[i]<-results_df_500bp$gene[i]
    results_proseq_minus$V7[i]<-as.integer(results_df_500bp$gene_start[i])
    results_proseq_minus$V8[i]<-results_df_500bp$strand[i]
  } else if (results_df_500bp$strand[i]=="-"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_end[i])-550),(as.integer(results_df_500bp$gene_end[i])+100)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_proseq_minus$V1[i]<-intermediate[1]
    results_proseq_minus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_proseq_minus$V3[i]<-intermediate[3]
    results_proseq_minus$V4[i]<-intermediate[4]
    results_proseq_minus$V5[i]<-intermediate[5]
    results_proseq_minus$V6[i]<-results_df_500bp$gene[i]
    results_proseq_minus$V7[i]<-as.integer(results_df_500bp$gene_end[i])
    results_proseq_minus$V8[i]<-results_df_500bp$strand[i]
  }
}
row.names(results_proseq_minus)<-NULL

test<-read.delim("5.Merged_replicates/proseq_GCB/normalized_upg_corrected_50486_50487_merged_plus.bedgraph", header=FALSE) #load the PROSeq plus strand
#test<-read.delim(args[6], header=FALSE) #load the PROSeq plus strand

results_proseq_plus<-matrix(NA, nrow=length(results_df_500bp$gene),ncol=8)
results_proseq_plus<-as.data.frame(results_proseq_plus)
intermediate<-as.data.frame(matrix(0, nrow=1, ncol=5))
for (i in 1:length(results_df_500bp$gene)){
  if (results_df_500bp$strand[i]=="+"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_start[i])-100),(as.integer(results_df_500bp$gene_start[i])+550)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_proseq_plus$V1[i]<-intermediate[1]
    results_proseq_plus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_proseq_plus$V3[i]<-intermediate[3]
    results_proseq_plus$V4[i]<-intermediate[4]
    results_proseq_plus$V5[i]<-intermediate[5]
    results_proseq_plus$V6[i]<-results_df_500bp$gene[i]
    results_proseq_plus$V7[i]<-as.integer(results_df_500bp$gene_start[i])
    results_proseq_plus$V8[i]<-results_df_500bp$strand[i]
  } else if (results_df_500bp$strand[i]=="-"){
    intermediate<-Internal.trascription.sites(1,0,0,0,20,results_df_500bp$chr[i],(as.integer(results_df_500bp$gene_end[i])-550),(as.integer(results_df_500bp$gene_end[i])+100)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
    results_proseq_plus$V1[i]<-intermediate[1]
    results_proseq_plus$V2[i]<-results_df_500bp$chr[i] ###CHECK THIS###
    results_proseq_plus$V3[i]<-intermediate[3]
    results_proseq_plus$V4[i]<-intermediate[4]
    results_proseq_plus$V5[i]<-intermediate[5]
    results_proseq_plus$V6[i]<-results_df_500bp$gene[i]
    results_proseq_plus$V7[i]<-as.integer(results_df_500bp$gene_end[i])
    results_proseq_plus$V8[i]<-results_df_500bp$strand[i]
  }
}
row.names(results_proseq_plus)<-NULL

colnames(results_proseq_minus)<-c('hits','chr','starting_position','ending_position','signal_strength','gene','annotated_start','strand')
colnames(results_proseq_plus)<-c('hits','chr','starting_position','ending_position','signal_strength','gene','annotated_start','strand')

results_proseq_plus$distance_from_start<-list(NA)
for (i in 1:length(results_proseq_plus$hits)){
  if (results_proseq_plus$strand[i]=='+'){
    results_proseq_plus$distance_from_start[[i]]<-(as.numeric(results_proseq_plus$ending_position[[i]])-as.numeric(results_proseq_plus$annotated_start[i]))
  } else {
    results_proseq_plus$distance_from_start[[i]]<-as.numeric(results_proseq_plus$annotated_start[i])-(as.numeric(results_proseq_plus$ending_position[[i]]))
  }
}

results_proseq_minus$distance_from_start<-list(NA)
for (i in 1:length(results_proseq_minus$hits)){
  if (results_proseq_minus$strand[i]=='+'){
    results_proseq_minus$distance_from_start[[i]]<-(as.numeric(results_proseq_minus$ending_position[[i]])-as.numeric(results_proseq_minus$annotated_start[i]))
  } else {
    results_proseq_minus$distance_from_start[[i]]<-as.numeric(results_proseq_minus$annotated_start[i])-(as.numeric(results_proseq_minus$ending_position[[i]]))
  }
}

results_proseq_sense_minus<-results_proseq_minus[which(results_proseq_minus$strand=="-"),]
results_proseq_antisense_for_plus_genes<-results_proseq_minus[which(results_proseq_minus$strand=="+"),]

results_proseq_sense_plus<-results_proseq_plus[which(results_proseq_plus$strand=="+"),]
results_proseq_antisense_for_minus_genes<-results_proseq_plus[which(results_proseq_plus$strand=="-"),]

results_proseq_sense<-rbind(results_proseq_sense_minus,results_proseq_sense_plus)
results_proseq_sense$index<-as.numeric(rownames(results_proseq_sense)) #getting the index
results_proseq_sense<-results_proseq_sense[order(results_proseq_sense$index),] #ordering according to the index

results_proseq_antisense<-rbind(results_proseq_antisense_for_minus_genes,results_proseq_antisense_for_plus_genes)
results_proseq_antisense$index<-as.numeric(rownames(results_proseq_antisense)) #getting the index
results_proseq_antisense<-results_proseq_antisense[order(results_proseq_antisense$index),] #ordering according to the index

results_proseq_antisense$signal_strength<-lapply(results_proseq_antisense$signal_strength,as.numeric)
results_proseq_sense$signal_strength<-lapply(results_proseq_sense$signal_strength,as.numeric)

results_proseq_antisense$distance_from_start[which(is.na(results_proseq_antisense$distance_from_start))]<-0
results_proseq_sense$distance_from_start[which(is.na(results_proseq_sense$distance_from_start))]<-0
### --- ###


if (categorizing==1){
  for (i in 1:nrow(results_df_500bp)){
    if (results_df_500bp$gene[i] %in% refseq$Geneid){
      results_df_500bp$category[i]<-as.character(refseq$exp_lvl[which(refseq$Geneid==results_df_500bp$gene[i])])
    }
  }
  quant_mut_str<-summary(results_df_500bp$sum_of_freq)
  attr(quant_mut_str, "names") <- NULL

  # categorizing according to the sum of mutation frequency and gene expression levels
  for (i in 1:length(results_df_500bp$gene)){

    if (results_df_500bp$sum_of_freq[i]>=quant_mut_str[2] & results_df_500bp$sum_of_freq[i]<=quant_mut_str[5]){
      results_df_500bp$category[i]<-paste0(results_df_500bp$category[i],'m')
    } else if (results_df_500bp$sum_of_freq[i]<quant_mut_str[2]){
      results_df_500bp$category[i]<-paste0(results_df_500bp$category[i],'l')
    } else if (results_df_500bp$sum_of_freq[i]>quant_mut_str[5]){
      results_df_500bp$category[i]<-paste0(results_df_500bp$category[i],'h')
    }

    results_sense$category<-results_df_500bp$category
    results_antisense$category<-results_df_500bp$category
    results_proseq_sense$category<-results_df_500bp$category
    results_proseq_antisense$category<-results_df_500bp$category
  }
  #creating a list of data.frames deriving from the different categories
  categories<-sort(unlist(unique(results_df_500bp$category)))
  df_list<-list(NA,length(categories))
  j=1
  for (category in categories){
    df_list[[j]]<-results_df_500bp[which(results_df_500bp$category==category),]
    j=j+1
  }

  df_sense_list<-list(NA,length(categories))
  j=1
  for (category in categories){
    df_sense_list[[j]]<-results_sense[which(results_sense$category==category),]
    j=j+1
  }

  df_antisense_list<-list(NA,length(categories))
  j=1
  for (category in categories){
    df_antisense_list[[j]]<-results_antisense[which(results_antisense$category==category),]
    j=j+1
  }

  df_proseq_sense_list<-list(NA,length(categories))
  j=1
  for (category in categories){
    df_proseq_sense_list[[j]]<-results_proseq_sense[which(results_proseq_sense$category==category),]
    j=j+1
  }

  df_proseq_antisense_list<-list(NA,length(categories))
  j=1
  for (category in categories){
    df_proseq_antisense_list[[j]]<-results_proseq_antisense[which(results_proseq_antisense$category==category),]
    j=j+1
  }

}

par(mar=c(5.1,4.1,4.1,2.1))


# the positions of WRCY motifs
motifs<-read.csv("0.External_input_data/wrcy_sites_newpriB_tested_tpm0.txt", sep='\t', header=FALSE)
#motifs<-read.csv(args[7], sep='\t', header=FALSE)

colnames(motifs)<-c('chr','end','strand','base','zero','motif','gene','exp_lvl')
motifs$strand[which(motifs$strand=='1')]<-'+'
motifs$strand[which(motifs$strand=='-1')]<-'-'

#motif_colors=c(col=rgb(0.6,0.2,0, alpha=0.8),col=rgb(0,0.4,0.1, alpha = 0.7),rgb(0.6,0,0.4, alpha=0.8),rgb(0,0.4,0.6, alpha = 0.7))
motif_colors=c(rgb(0,0.4,0.6, alpha = 0.7),rgb(0,0.4,0.6, alpha = 0.7),rgb(0,0.4,0.6, alpha = 0.7),rgb(0,0.4,0.6, alpha = 0.7))

# plotting time

a3 = c(0.1,10,6.7,7.44)/10
b3 = c(0.1,10,6.06,6.7)/10
all3 = c(0.1,10,6.06,7.44)/10
a4 = c(0.1,10,5.42,6.16)/10
b4 = c(0.1,10,4.78,5.42)/10
all4 = c(0.1,10,4.78,6.16)/10
a5 = c(0.1,10,4.14,4.88)/10
b5 = c(0.1,10,3.5,4.14)/10
all5 = c(0.1,10,3.5,4.88)/10
a6 = c(0.1,10,2.86,3.6)/10
b6 = c(0.1,10,2.22,2.86)/10
all6 = c(0.1,10,2.22,3.6)/10

#pdf(paste0('../../10.Graphs_mutations/mut_freq_TSSs_pauses_AID_off.pdf'))
if (categorizing==1){
  for (i in (1:length(df_list))){
    pdf(paste0('10.Graphs_mutations/',df_list[[i]]$category[[1]],'_mut_freq_TSSs_pauses_AID_off_updated.pdf'), pointsize = 12)
    #par(mfrow=c(6,1), mai = c(0, 0, 0.15, 0.05), oma=c(2,2,0.3,0))
    if(length(df_list[[i]]$gene)>0){
      for (k in seq(1,length(df_list[[i]]$gene),by=1)){
        layout(matrix(c(1,2,3,4,5), nrow = 5, ncol = 1, byrow = T))
        par(mar=c(0,6,1.5,0.5))
        par(fig=c(0.1,10,8.62,10)/10)


        plot(log10(abs(df_sense_list[[i]]$signal_strength[[k]])+1)~(df_sense_list[[i]]$distance_from_start[[k]]), pch=20, col = rgb(0.3, 0.3, 0.8, alpha = 0), lwd=1.5,
             ylim=c(-max(log10(abs(df_antisense_list[[i]]$signal_strength[[k]])+1)),max(log10(abs(df_sense_list[[i]]$signal_strength[[k]])+1))),
             xlab='',ylab='', las=1,
             xlim=c(-100,550), cex=1, xaxt='n', cex.axis=1)
        axis(side=1, labels = FALSE)
        title(ylab='GCB PRO-cap\n(log10 RPM)', line=4, cex.lab=1)
        segments((df_sense_list[[i]]$distance_from_start[[k]]), 0, (df_sense_list[[i]]$distance_from_start[[k]]), log10(abs(df_sense_list[[i]]$signal_strength[[k]])+1))
        #points(-log10(abs(df_antisense_list[[i]]$signal_strength[[k]])+1)~(df_antisense_list[[i]]$distance_from_start[[k]]), pch=20, col='lightblue4',lwd=1.5, cex=1)
        segments((df_antisense_list[[i]]$distance_from_start[[k]]), 0, (df_antisense_list[[i]]$distance_from_start[[k]]), -log10(abs(df_antisense_list[[i]]$signal_strength[[k]])+1))
        mtext(df_sense_list[[i]]$gene[k], side=3, cex = 0.7, font =3)

        par(fig=c(0.1,10,7.34,8.72)/10)
        par(new=T)

        plot(log10(abs(df_proseq_sense_list[[i]]$signal_strength[[k]])+1)~(df_proseq_sense_list[[i]]$distance_from_start[[k]]), pch=20, col = rgb(0.3, 0.3, 0.8, alpha = 0), lwd=1.5,
             ylim=c(-max(log10(abs(df_proseq_antisense_list[[i]]$signal_strength[[k]])+1)),max(log10(abs(df_proseq_sense_list[[i]]$signal_strength[[k]])+1))),
             xlab=NA,ylab='', las=1,
             xlim=c(-100,550), cex=1, xaxt='n', cex.axis=1)
        axis(side=1, labels = FALSE)
        title(ylab='GCB PRO-seq\n(log10 RPM)', line=4, cex.lab=1)
        segments((df_proseq_sense_list[[i]]$distance_from_start[[k]]), 0, (df_proseq_sense_list[[i]]$distance_from_start[[k]]), log10(abs(df_proseq_sense_list[[i]]$signal_strength[[k]])+1))
        #points(-log10(abs(df_proseq_antisense_list[[i]]$signal_strength[[k]])+1)~(df_proseq_antisense_list[[i]]$distance_from_start[[k]]), pch=20, col='lightblue4',lwd=1.5, cex=1)
        segments((df_proseq_antisense_list[[i]]$distance_from_start[[k]]), 0, (df_proseq_antisense_list[[i]]$distance_from_start[[k]]), -log10(abs(df_proseq_antisense_list[[i]]$signal_strength[[k]])+1))


        par(mar=c(0,6,1.5,0.5))
        #par(mar=c(0,6,0,0.5))

        par(fig=c(0.010,1,0.606,0.744))
        par(new=T)

        #getting the coordinates from the start of the gene of the WRCY motifs
        if (df_sense_list[[i]]$strand[k]=='+'){
          motif_distance_from_start_sense<-motifs$end[which(motifs$gene==df_sense_list[[i]]$gene[k] & motifs$base=='C')]-refseq$Start[which(refseq$Geneid==df_sense_list[[i]]$gene[k])]
        }else if (df_sense_list[[i]]$strand[k]=='-'){
          motif_distance_from_start_sense<-refseq$End[which(refseq$Geneid==df_sense_list[[i]]$gene[k])]-motifs$end[which(motifs$gene==df_sense_list[[i]]$gene[k] & motifs$base=='C')]
        }
        if (df_sense_list[[i]]$strand[k]=='+'){
          motif_distance_from_start_antisense<-motifs$end[which(motifs$gene==df_sense_list[[i]]$gene[k] & motifs$base=='G')]-refseq$Start[which(refseq$Geneid==df_sense_list[[i]]$gene[k])]
        }else if (df_sense_list[[i]]$strand[k]=='-'){
          motif_distance_from_start_antisense<-refseq$End[which(refseq$Geneid==df_sense_list[[i]]$gene[k])]-motifs$end[which(motifs$gene==df_sense_list[[i]]$gene[k] & motifs$base=='G')]
        }


        plot((abs(df_list[[i]]$frequencies[[k]]))~df_list[[i]]$all_positions_from_start[[k]], pch=20, col = rgb(0.3, 0.3, 0.8, alpha = 0),
             xlab='Position in relation to the annotated start of the gene (bp)',ylab='', yaxt='n',
             xlim=c(-100,550), cex=1, cex.axis=1, xaxt='n',
             ylim=c(0,1.01*max(abs(df_list[[i]]$frequencies[[k]]))))
        axis(side=1, cex.axis=1, cex=1, labels=F)
        TSS_zero<-axTicks(1, axp = NULL, usr = NULL, log = NULL, nintLog = NULL)
        axis(side=1, cex.axis=1, cex=1, at=TSS_zero[which(TSS_zero!=0)])
        mtext('TSS',side=1, cex=0.6, at=0, padj = 1.9)
        mtext('(bp)',side=1, cex=0.6, at=550, padj = 1.7)

        segments((df_list[[i]]$all_positions_from_start[[k]]), 0, (df_list[[i]]$all_positions_from_start[[k]]), (abs(df_list[[i]]$frequencies[[k]])))
        axis(side=2, labels =T, cex.axis=1, las=1)
        title(ylab='Mutation freq\n(log10)', line=4, cex.lab=1)
        points(motif_distance_from_start_sense,rep(0,length(motif_distance_from_start_sense)), col=rgb(1,0,0, alpha=0.9), pch=20, cex=0.8)
        points(motif_distance_from_start_antisense,rep(0,length(motif_distance_from_start_antisense)), col=rgb(1,0,0, alpha = 0.9), pch=20, cex=0.8)

        #points(cold_motif_distance_from_start_sense,rep(0,length(cold_motif_distance_from_start_sense)), col=rgb(0.6,0,0.4, alpha=0.8), pch=20, cex=1)
        #points(cold_motif_distance_from_start_antisense,rep(0,length(cold_motif_distance_from_start_antisense)), col=rgb(0,0.4,0.6, alpha = 0.7), pch=20, cex=1)

        legend('topleft', legend = c('WRCY/RGWY'), border=c(0,0,0,0),  box.col = rgb(1,1,1,0),bg =  rgb(1,1,1,0), col=c(rgb(1,0,0, alpha = 0.9)), cex=0.9, pch=19)
      }
    }
    dev.off()
  }
}
#dev.off()
par(mfrow=c(1,1))
