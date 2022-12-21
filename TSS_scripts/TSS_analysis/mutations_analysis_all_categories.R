if (!require(readODS)) install.packages('readODS')
library('readODS')

if (!require(ggplot2)) install.packages('ggplot2')
library('ggplot2')

if (!require(vioplot)) install.packages('vioplot')
library('vioplot')

if (!require(ggpubr)) install.packages('ggpubr')
library('ggpubr')

if (!require(tidyverse)) install.packages('tidyverse')
library('dplyr')

if (!require(stringr)) install.packages('stringr')
library('stringr')

if (!require(openxlsx)) install.packages('openxlsx')
library('openxlsx')

if (!require(gridExtra)) install.packages('gridExtra')
library("gridExtra")

if (!require(ggforce)) install.packages('ggforce')
library('ggforce')

if (!require(OneR)) install.packages('OneR')
library('OneR')

# replicate_name<-readline(prompt = 'Choose a name for the replicate (procap, proseq_newpriB, etc): ')
replicate_name<-"proseq_GCB"

#importing the replicates
# test1<-read.delim(file.choose(), header=FALSE) #load the bedgraph file of one of the strands
# test2<-read.delim(file.choose(), header=FALSE) #load the bedgraph file of the other strand
test1<-read.delim("5.Merged_replicates/proseq_GCB/normalized_upg_corrected_50486_50487_merged_minus.bedgraph", header = F)
test2<-read.delim("5.Merged_replicates/proseq_GCB/normalized_upg_corrected_50486_50487_merged_minus.bedgraph", header = F)

# use_motifs<-readline(prompt = 'Are you using the WRCY/RGYW motifs? Enter 1 if yes, or any number for no: ')
use_motifs<-0
# non_mutated<-read.csv(file.choose(), header=FALSE, sep='\t') #choose the file with the C/Gs in the tested genes
non_mutated<-read.csv("0.External_input_data/c_sites_GCB_tested_tpm20.txt", header = F, sep = "\t")
# distance_cutoff<-as.integer(readline(prompt = 'Enter the maximum distance upstream and downstream from the C/G site: '))
distance_cutoff<-10

#importing the sheets
# mutations<-read.ods(file.choose()) #the mutations data are found in our case in the file PerBase_CG_targets_dKO_PavriLab.ods
mutations<-read.ods("0.External_input_data/PerBase_CG_targets_dKO_PavriLab.ods")#the mutations data are found in our case in the file PerBase_CG_targets_dKO_PavriLab.ods
mutations_backup <- mutations

#formatting the data frames
for (i in 1:length(mutations)){
  mutations[[i]]<-mutations[[i]][c(-5,-6,-7,-8,-12,-13)]
  mutations[[i]]<-mutations[[i]][-1,]
  colnames(mutations[[i]])<-c('chr','end','strand','base','motif','gene','freqA','freqT')
}


#combining the data frames
conc_mutations<-mutations[[1]]
for (i in 2:length(mutations)){
  conc_mutations<-rbind(conc_mutations,mutations[[i]])
}
conc_mutations_C<-conc_mutations[which(conc_mutations$base=='C' & conc_mutations$freqT>0),-7]
conc_mutations_G<-conc_mutations[which(conc_mutations$base=='G' & conc_mutations$freqA>0),-8]
colnames(conc_mutations_C)<-c('chr','end','strand','base','motif','gene','freq')
colnames(conc_mutations_G)<-colnames(conc_mutations_C)

conc_mutations<-rbind(conc_mutations_C,conc_mutations_G)

conc_mutations<-conc_mutations[,c('chr','end','strand','base','freq','motif','gene')]

if (use_motifs=='1'){
  all_motifs<-unique(conc_mutations$motif)
  rgyw<-str_subset(all_motifs,'[AG]G[CT][AT]')
  wrcy<-str_subset(all_motifs,'[AT][AG]C[CT]')
  two_motifs<-c(rgyw,wrcy)
  conc_mutations<-conc_mutations[which(conc_mutations$motif %in% two_motifs),]
}

colnames(non_mutated)<-c('chr','end','strand','base','freq','motif','gene','exp_lvl')

#formatting the genes
conc_mutations$gene[]<-lapply(conc_mutations$gene, gsub, pattern='_[A-Za-z0-9]+',replacement='')
conc_mutations$gene[]<-lapply(conc_mutations$gene, gsub, pattern='-[0-9]+',replacement='')

non_mutated$category<-0
for (genes in non_mutated$gene){
  if (genes %in% conc_mutations$gene){
    non_mutated$category[which(non_mutated$gene==genes)]<-1
  }
}

for (i in 1:nrow(conc_mutations)){
  print(i)
  non_mutated$category[which(non_mutated$chr==conc_mutations$chr[i] & non_mutated$end==conc_mutations$end[i])]<-2
  non_mutated$freq[which(non_mutated$chr==conc_mutations$chr[i] & non_mutated$end==conc_mutations$end[i])]<-
    conc_mutations$freq[i]
}

#we will add a starting position for the mutations (to be in accordance with the format we had be using so far)
non_mutated$end<-as.numeric(non_mutated$end)
non_mutated$start<-non_mutated$end-1
#we will also change the way we address the strand (we will use + and -)
non_mutated$strand[which(non_mutated$strand=='1')]<-'+'
non_mutated$strand[which(non_mutated$strand=='-1')]<-'-'

#reordering the columns
non_mutated<-non_mutated[c(1,10,2,3,4,5,6,7,8,9)]

#empty data frame that we will use to store the data from the expanding region around the mutations
exp_bubble<-as.data.frame(matrix(data=NA, nrow=distance_cutoff+1, ncol=3))
colnames(exp_bubble)<-c('region_size','tss_found','tss_not_found')
a<-distance_cutoff
#importing one strand of the replicate
test<-test1
# rm(test1)
intermediate<-as.data.frame(matrix(0, nrow=1, ncol=5))
source("TSS_scripts/TSS_analysis/Internal.transcription.sites.R") #sourcing the Internal.transcription.sites.R function
non_mutated$chr<-as.character(non_mutated$chr)

for (i in 1:length(non_mutated[,1])){
  intermediate<-Internal.trascription.sites(1,0,0,0,20,non_mutated[i,1],(as.integer(non_mutated[i,2])-distance_cutoff),(as.integer(non_mutated[i,3]))+distance_cutoff) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
  non_mutated$hits[i]<-intermediate[1]
  non_mutated$starting_position[i]<-intermediate[3]
  non_mutated$ending_position[i]<-intermediate[4]
  non_mutated$signal_strength[i]<-intermediate[5]
  cat(100*i/length(non_mutated[,1]),'% of the regions processed\n')
}
#importing the other strand of the replicate
test<-test2
# rm(test2)

intermediate<-as.data.frame(matrix(0, nrow=1, ncol=5))

for (i in 1:length(non_mutated[,1])){
  intermediate<-Internal.trascription.sites(1,0,0,0,20,non_mutated[i,1],(as.integer(non_mutated[i,2])-distance_cutoff),(as.integer(non_mutated[i,3])+distance_cutoff)) #the first input value should be 1, the second is the cutoff percentage, the third can be 1 or 0 if you want to see some plots (1) or not (0), the fourth can be 1 (and will create "super-regions", which means it will connect peaks the distance between which can should be less than the fifth parameter) or 0 (you don't want to assemble "super-regions"), the sixth is the chromosome, the seventh and eighth are the coordinates
  non_mutated$hits2[i]<-intermediate[1]
  non_mutated$starting_position2[i]<-intermediate[3]
  non_mutated$ending_position2[i]<-intermediate[4]
  non_mutated$signal_strength2[i]<-intermediate[5]
  cat(100*i/length(non_mutated[,1]),'% of the regions processed\n')
}
  
# write the results to a file
backup <- non_mutated
save.image(file="pavri_Cs_and_Gs.RData")

#we are treating each mutation point as a unique category with its own characteristics (that is position and frequency of mutation)

#adding the hits from the + and - strands
non_mutated$total_hits<-as.numeric(non_mutated$hits)+as.numeric(non_mutated$hits2)


exp_bubble$region_size[distance_cutoff+1]<-distance_cutoff
exp_bubble$tss_found[distance_cutoff+1]<-100*length(non_mutated$start[which(non_mutated$total_hits>0)])/length(non_mutated[,1]) #how many mutations had at least one TSS nearby
exp_bubble$tss_not_found[distance_cutoff+1]<-100*length(non_mutated$start[which(non_mutated$total_hits==0)])/length(non_mutated[,1]) #how many mutations didn't have any TSS nearby
length(non_mutated[,1]) #total mutations

non_mutated$signal_strength<-lapply(non_mutated$signal_strength, as.numeric)
non_mutated$signal_strength2<-lapply(non_mutated$signal_strength2, as.numeric)
sum_of_strength1<-lapply(non_mutated$signal_strength, sum)
sum_of_strength1<-lapply(sum_of_strength1, abs)
sum_of_strength2<-lapply(non_mutated$signal_strength2, sum)
sum_of_strength2<-lapply(sum_of_strength2, abs)
sum_of_strength1<-as.numeric(sum_of_strength1)
sum_of_strength2<-as.numeric(sum_of_strength2)
non_mutated$sum_of_strength<-sum_of_strength1+sum_of_strength2

non_mutated$total_hits<-as.integer(non_mutated$total_hits)
non_mutated$freq<-as.numeric(non_mutated$freq)

non_mutated$category<-as.factor(non_mutated$category)

non_mutated$exp_lvl <- factor(non_mutated$exp_lvl, levels = c('low','moderate','high'))

my_comparison<-list(c('0','1'),c('1','2'),c('0','2'))

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'.pdf'))
pdf(paste0(replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'.pdf'))

bp <- ggplot(non_mutated, aes(x=category, y=log2(total_hits+1), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
  #stat_compare_means(label.y = max(non_mutated$total_hits)+3.5)
bp$layers[[2]]$aes_params$textsize <- 2
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(non_mutated, aes(x=category, y=log2(total_hits+1), group=category)) +
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black") + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
  #stat_compare_means(label.y = max(non_mutated$total_hits)+3.5)
bp$layers[[3]]$aes_params$textsize <- 2
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(non_mutated, aes(x=category, y=(log10(sum_of_strength+0.01)), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('(log10(sum of all signal strengths+0.01))') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
  #stat_compare_means(label.y = max(log10(non_mutated$sum_of_strength+1))+1.5)
bp$layers[[2]]$aes_params$textsize <- 2
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(non_mutated, aes(x=category, y=log10(sum_of_strength+0.01), group=category)) + 
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black")  + ylab('log10(sum of all signal strengths+0.01)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
  #stat_compare_means(label.y = max(log10(non_mutated$sum_of_strength+1))+1.5)
bp$layers[[3]]$aes_params$textsize <- 2
print(bp + facet_grid(exp_lvl ~ .))

dev.off()

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_cposition.pdf'))
pdf(paste0(replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_cposition.pdf'))


bp <- ggplot(non_mutated, aes(x=category, y=log2(total_hits+1), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(non_mutated$total_hits)+3.5)
bp$layers[[2]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(non_mutated, aes(x=category, y=log2(total_hits+1), group=category)) +
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black") + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(non_mutated$total_hits)+3.5)
bp$layers[[3]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(non_mutated, aes(x=category, y=log10(sum_of_strength+0.01), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log10(sum of all signal strengths+0.01)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(non_mutated$sum_of_strength+1))+1.3)
bp$layers[[2]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(non_mutated, aes(x=category, y=log10(sum_of_strength+0.01), group=category)) + 
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black")  + ylab('log10(sum of all signal strengths+0.01)') + stat_compare_means(tip.length=0,comparison = my_comparison)# symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(non_mutated$sum_of_strength+1))+1.3)
bp$layers[[3]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))

dev.off()

#creating the sense dataframe
sense_minus<-non_mutated[which(non_mutated$strand=='-'),] %>% select(-hits2, -starting_position2, -ending_position2, -signal_strength2, -total_hits, -sum_of_strength)
sense_plus<-non_mutated[which(non_mutated$strand=='+'),] %>% select(-hits, -starting_position, -ending_position, -signal_strength, -total_hits, -sum_of_strength)
sense_plus<-sense_plus %>% rename(hits=hits2, starting_position=starting_position2, ending_position=ending_position2, signal_strength=signal_strength2)
sense<-rbind(sense_minus,sense_plus)

sense$sum_of_strength<-0
sense$sum_of_strength<-lapply(sense$signal_strength, sum)
sense$sum_of_strength<-lapply(sense$sum_of_strength, abs)
sense$hits<-as.numeric(sense$hits)
sense$sum_of_strength<-as.numeric(sense$sum_of_strength)

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_sense.pdf'))
pdf(paste0(replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_sense.pdf'))

bp <- ggplot(sense, aes(x=category, y=log2(hits+1), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(sense$hits)+3.5)
bp$layers[[2]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(sense, aes(x=category, y=log2(hits+1), group=category)) +
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black") + ylab('log2(number of signals)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(sense$hits)+3.5)
bp$layers[[3]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(sense, aes(x=category, y=log10(sum_of_strength+0.01), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log10(sum of all signal strengths+0.01)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(sense$sum_of_strength+1))+1.5)
bp$layers[[2]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(sense, aes(x=category, y=log10(sum_of_strength+0.01), group=category)) + 
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black")  + ylab('log10(sum of all signal strengths+0.01)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(sense$sum_of_strength+1))+1.5)
bp$layers[[3]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))

dev.off()

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_cposition_sense.pdf'))
pdf(paste0(replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_cposition_sense.pdf'))


bp <- ggplot(sense, aes(x=category, y=log2(hits+1), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(sense$hits)+3.5)
bp$layers[[2]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(sense, aes(x=category, y=log2(hits+1), group=category)) +
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black") + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(sense$hits)+3.5)
bp$layers[[3]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(sense, aes(x=category, y=log10(sum_of_strength+0.01), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log10(sum of all signal strengths+0.01)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(sense$sum_of_strength+1))+1.5)
bp$layers[[2]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(sense, aes(x=category, y=log10(sum_of_strength+0.01), group=category)) + 
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black")  + ylab('log10(sum of all signal strengths+0.01)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(sense$sum_of_strength+1))+1.5)
bp$layers[[3]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))

dev.off()


#creating the antisense dataframe
antisense_plus<-non_mutated[which(non_mutated$strand=='+'),] %>% select(-hits2, -starting_position2, -ending_position2, -signal_strength2, -total_hits, -sum_of_strength)
antisense_minus<-non_mutated[which(non_mutated$strand=='-'),] %>% select(-hits, -starting_position, -ending_position, -signal_strength, -total_hits, -sum_of_strength)
antisense_minus<-antisense_minus %>% rename(hits=hits2, starting_position=starting_position2, ending_position=ending_position2, signal_strength=signal_strength2)
antisense<-rbind(antisense_minus,antisense_plus)

antisense$sum_of_strength<-0
antisense$sum_of_strength<-lapply(antisense$signal_strength, sum)
antisense$sum_of_strength<-lapply(antisense$sum_of_strength, abs)
antisense$hits<-as.numeric(antisense$hits)
antisense$sum_of_strength<-as.numeric(antisense$sum_of_strength)

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_antisense.pdf'))
pdf(paste0(replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_antisense.pdf'))


bp <- ggplot(antisense, aes(x=category, y=log2(hits+1), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(antisense$hits)+3.5)
bp$layers[[2]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(antisense, aes(x=category, y=log2(hits+1), group=category)) +
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black") + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(antisense$hits)+3.5)
bp$layers[[3]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))



bp <- ggplot(antisense, aes(x=category, y=log10(sum_of_strength+0.001), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log10(sum of all signal strengths+0.001)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(antisense$sum_of_strength+1))+1.5)
bp$layers[[2]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))


bp <- ggplot(antisense, aes(x=category, y=log10(sum_of_strength+0.001), group=category)) + 
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black")  + ylab('log10(sum of all signal strengths+0.001)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(antisense$sum_of_strength+1))+1.5)
bp$layers[[3]]$aes_params$textsize <- 2.5
print(bp + facet_grid(exp_lvl ~ .))

dev.off()

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_cposition_antisense.pdf'))
pdf(paste0(replicate_name,'mutated_vs_unmutated_distance_',distance_cutoff,'_cposition_antisense.pdf'))


bp <- ggplot(antisense, aes(x=category, y=log2(hits+1), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(antisense$hits)+3.5)
bp$layers[[2]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(antisense, aes(x=category, y=log2(hits+1), group=category)) +
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black") + ylab('log2(number of signals+1)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(antisense$hits)+3.5)
bp$layers[[3]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(antisense, aes(x=category, y=log10(sum_of_strength+0.001), group=category)) + 
  geom_boxplot(aes(fill=category)) + ylab('log10(sum of all signal strengths+0.001)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(antisense$sum_of_strength+1))+1.5)
bp$layers[[2]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))


bp <- ggplot(antisense, aes(x=category, y=log10(sum_of_strength+0.001), group=category)) + 
  geom_violin(aes(fill=category)) + stat_summary(fun.y=median, geom="point", size=2, color="black")  + ylab('log10(sum of all signal strengths+0.001)') + stat_compare_means(tip.length=0,comparison = my_comparison)#, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) #+
#stat_compare_means(label.y = max(log10(antisense$sum_of_strength+1))+1.5)
bp$layers[[3]]$aes_params$textsize <- 1.3
print(bp + facet_grid(exp_lvl + base ~ .))

dev.off()

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'freq_of_TSS_around_sites.pdf'))
pdf(paste0(replicate_name,'freq_of_TSS_around_sites.pdf'))

#metaplot of TSSs around the mutations, sense and antisense
split_df_list<-split(non_mutated, list(non_mutated$exp_lvl, non_mutated$category), drop = FALSE)

par(mfrow=c(3,3))

for (m in 1:length(split_df_list)){
  sense_list<-matrix(NA, nrow=length(split_df_list[[m]]$base), ncol=1)
  antisense_list<-matrix(NA, nrow=length(split_df_list[[m]]$base), ncol=1)
  
  sense_count=0
  antisense_count=0
  
  for (l in 1:length(split_df_list[[m]]$base)){
    if (split_df_list[[m]]$strand[l]=='+'){
      if (split_df_list[[m]]$hits2[l]>0){
        sense_list[l]<-list(split_df_list[[m]]$ending_position2[[l]] - split_df_list[[m]]$end[l])
        sense_count=sense_count+as.integer(split_df_list[[m]]$hits2[l])
      }
      if (split_df_list[[m]]$hits[l]>0){
        antisense_list[l]<-list(split_df_list[[m]]$ending_position[[l]] - split_df_list[[m]]$end[l])
        antisense_count=antisense_count+as.integer(split_df_list[[m]]$hits[l])
      }
    }
    else if (split_df_list[[m]]$strand[l]=='-'){
      if (split_df_list[[m]]$hits[l]>0){
        sense_list[l]<-list(-split_df_list[[m]]$ending_position[[l]] + split_df_list[[m]]$end[l])
        sense_count=sense_count+as.integer(split_df_list[[m]]$hits[l])
      }
      if (split_df_list[[m]]$hits2[l]>0){
        antisense_list[l]<-list(-split_df_list[[m]]$ending_position2[[l]] + split_df_list[[m]]$end[l])
        antisense_count=antisense_count+as.integer(split_df_list[[m]]$hits2[l])
      }
    }
  }
  
  
  sense_table<-as.data.frame(table(unlist(sense_list)))
  sense_table$Freq<-as.numeric(sense_table$Freq)/(sense_count+antisense_count)
  
  antisense_table<-as.data.frame(table(unlist(antisense_list)))
  antisense_table$Freq<-as.numeric(antisense_table$Freq)/(sense_count+antisense_count)
  
  sense_table$Var1<-as.numeric(as.character(sense_table$Var1))
  antisense_table$Var1<-as.numeric(as.character(antisense_table$Var1))
  
  plot(sense_table$Freq~sense_table$Var1, type='l', col='blue', main=names(split_df_list)[m], ylab='Percentage',
       xlab='Distance from the site (nt)', ylim=c(min(antisense_table$Freq),max(sense_table$Freq)),
       lwd=1)
  lines(antisense_table$Freq~antisense_table$Var1, type='l', col='red', lwd=1)
  legend("right",c('Sense','Antisense'), fill=c('blue','red'), cex=0.7)
  abline(v=0, lty=2)
  mtext('C/G site',las=3, line=-7, padj=-0.5, cex=0.7)
}

dev.off()

# write.xlsx(non_mutated, file=paste0('../../10.Graphs_mutations/',replicate_name,'mutations_data_table_distance',distance_cutoff,'.xlsx'))
write.xlsx(non_mutated, file=paste0(replicate_name,'mutations_data_table_distance',distance_cutoff,'.xlsx'))

only_mutated<-non_mutated[which(non_mutated$category==2),]

rownames(only_mutated)<-NULL
for (i in 1:nrow(only_mutated)){
  only_mutated$freq[i]<-conc_mutations$freq[which(conc_mutations$chr==only_mutated$chr[i] & conc_mutations$end==only_mutated$end[i])][1]
}
only_mutated$freq<-as.numeric(only_mutated$freq)
boxplot(log10(only_mutated$freq)~only_mutated$total_hits)

# pdf(paste0('../../10.Graphs_mutations/',replicate_name,'frequency_of_mutations',distance_cutoff,'.pdf'))
pdf(paste0(replicate_name,'frequency_of_mutations',distance_cutoff,'.pdf'))

bp <- ggplot(only_mutated, aes(x=total_hits, y=log10(freq), group=total_hits)) + 
  geom_boxplot() + ylab('log10(frequency of mutation)') + xlab('number of signals')
print(bp)
print(bp + facet_grid(exp_lvl + base ~ .))

bp <- ggplot(only_mutated, aes(x=(sum_of_strength+1), y=(freq+1))) + 
  geom_point(aes(fill=sum_of_strength), alpha=0.1) + ylab('log10(frequency of mutation)') +
  ggpubr::color_palette("jco") + 
  facet_zoom(x = log10(sum_of_strength+1) < 0.1)+
  theme_bw() + xlab('log10(sum of RPM of all signals+1)') + theme(legend.position = 'none')
print(bp)
#print(bp + facet_grid(exp_lvl + base ~ .))

# NEW ANALYSIS
all_mutated<-non_mutated
all_mutated$freq <- as.numeric(all_mutated$freq)
all_mutated$bins<-bin(all_mutated$sum_of_strength, nbins=10, method='content')
all_mutated$bins_mutations<-bin(all_mutated$freq, nbins=10, method='content')

all_mutated$bins<-bin(all_mutated$sum_of_strength, nbins=5, method='length')
all_mutated$bins_mutations<-bin(all_mutated$freq, nbins=5, method='length')

all_mutated2<-non_mutated[which(non_mutated$category==2),]
all_mutated2$freq <- as.numeric(all_mutated2$freq)
all_mutated2$bins<-bin(all_mutated2$sum_of_strength, nbins=10, method='content')
all_mutated2$bins_mutations<-bin(all_mutated2$freq, nbins=10, method='content')

all_mutated2$bins <- as.character(all_mutated2$bins)

val <- sort(unique(all_mutated2$bins))
rpl <- unlist(str_split(val[1],","))[2]
val[1] <- paste0("[0,",rpl)

all_mutated2$bins[which(all_mutated2$bins==sort(unique(all_mutated2$bins))[1])] <- val[1]

all_mutated2$bins <- factor(all_mutated2$bins, levels = val)
bp <- ggplot(all_mutated2, aes(x=bins, y=log10(freq+1), group=bins)) + 
  geom_boxplot(outlier.shape = NA) + ylab('log10(frequency of mutation+1)') +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab('sum of RPM of all signals') + stat_compare_means() + 
  coord_cartesian(ylim = c(0,0.001))

pdf("RPM_vs_freq_Cs_and_Gs.pdf", width = 10, height = 7)
print(bp)
dev.off()
# bp <- ggplot(all_mutated2, aes(x=bins, y=log10(freq+1), group=bins)) + 
#   geom_violin() + ylab('log10(frequency of mutation + 1)') + 
#   theme_minimal() +
#   coord_cartesian(ylim = c(0,0.005)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   xlab('total RPM') + stat_compare_means() +
#   stat_summary(fun = "mean", geom = "point",
#                colour = "black") +
#   stat_summary(fun = "median", geom = "point",
#                colour = "blue")
# print(bp)

p_results <- pairwise.wilcox.test(all_mutated2$freq, all_mutated2$bins, p.adjust.method = "bonferroni")
write.table(p_results$p.value, file = "Wilcoxon_p_adjusted_bonferroni_RPMbins.tsv", quote = F, col.names = NA, row.names = T, sep = "\t")
stats <- all_mutated2%>%
  group_by(bins)%>%
  summarise(Mean=mean(freq), Median=median(freq), Std=sd(freq))

stats_table <- as.data.frame(matrix(nrow = nrow(stats), ncol = nrow(stats)))
for (i in 1:nrow(stats)){
  colnames(stats_table)[i] <- as.character(stats$bins[i])
  for (j in 1:nrow(stats)){
    stats_table[i,j] <- (stats$Mean[i]/stats$Mean[j])
    if (i==1){
      rownames(stats_table)[j] <- as.character(stats$bins[j])
    }
  }
}
write.table(stats_table, file = "FC_mean_RPMbins.tsv", sep = "\t", quote = F, col.names = NA, row.names = T)

stats_table_median <- as.data.frame(matrix(nrow = nrow(stats), ncol = nrow(stats)))
for (i in 1:nrow(stats)){
  colnames(stats_table_median)[i] <- as.character(stats$bins[i])
  for (j in 1:nrow(stats)){
    stats_table_median[i,j] <- (stats$Median[i]/stats$Median[j])
    if (i==1){
      rownames(stats_table_median)[j] <- as.character(stats$bins[j])
    }
  }
}
write.table(stats_table_median, file = "FC_median_RPMbins.tsv", sep = "\t", quote = F, col.names = NA, row.names = T)

# bp <- ggplot(all_mutated, aes(x=bins_mutations, y=log10(sum_of_strength+1), group=bins_mutations)) + 
#   geom_boxplot() + ylab('log10(sum of RPM of all signals + 0.01)') + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   xlab('frequency of mutation') + stat_compare_means()
# print(bp)
# print(bp + facet_grid(exp_lvl + base ~ .))

# separate non-mutated Cs/Gs that are in mutated genes (class 1) and in non-mutated genes (class 0)
all_mutated_0_freq <- all_mutated[which(all_mutated$freq==0),]
all_mutated_with_mutations <- all_mutated[which(all_mutated$freq>0),]
all_mutated_0_freq$category2_freq <- NA
all_mutated_0_freq$category2_freq[which(all_mutated_0_freq$category==0)] <- "0 (non-AID-off-target)"
all_mutated_0_freq$category2_freq[which(all_mutated_0_freq$category==1)] <- "0 (AID-off-target)"
all_mutated_with_mutations$category2_freq <- bin(all_mutated_with_mutations$freq, nbins=10, method='content')
all_mutated_with_mutations$category2_freq <- as.character(all_mutated_with_mutations$category2_freq)

val <- sort(unique(all_mutated_with_mutations$category2_freq))
rpl <- unlist(str_split(val[1],","))[2]
val[1] <- paste0("(0,",rpl)

all_mutated_with_mutations$category2_freq[which(all_mutated_with_mutations$category2_freq==sort(unique(all_mutated_with_mutations$category2_freq))[1])] <- val[1]

results <- rbind(all_mutated_0_freq,all_mutated_with_mutations)
results$category2_freq <- factor(results$category2_freq, levels = c("0 (non-AID-off-target)","0 (AID-off-target)",val))

bp <- ggplot(results, aes(x=category2_freq, y=log10(sum_of_strength+1), group=category2_freq)) + 
  geom_boxplot(outlier.shape = NA) + ylab('log10(sum of RPM of all signals + 1)') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab('frequency of mutation') + stat_compare_means() +
  coord_cartesian(ylim=c(0,0.78))

pdf("mutfreq_vs_RPM_Cs_and_Gs.pdf", width = 10, height = 7)
print(bp)

bp <- ggplot(results, aes(x=category2_freq, y=log10(sum_of_strength+1), group=category2_freq)) + 
  geom_violin() + ylab('log10(sum of RPM of all signals + 1)') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab('frequency of mutation') + stat_compare_means() +
  coord_cartesian(ylim=c(0,2.3)) + 
  stat_summary(fun = "mean", geom = "point",
               colour = "black") +
  stat_summary(fun = "median", geom = "point",
               colour = "blue")
print(bp)
dev.off()

# vioplot(log10(results$sum_of_strength+1)~results$category2_freq)

# plot(log10(all_mutated$freq+1),log10(all_mutated$sum_of_strength+1),
#      xlim = c(0,0.001))
# 
# plot(log10(only_mutated$freq+1),log10(only_mutated$sum_of_strength+1),
#      xlim = c(0,0.03))
# dev.off()

# sum_of_signals_comparisons<-as.data.frame(compare_means(sum_of_strength~bins_mutations,data=only_mutated))
# colnames(sum_of_signals_comparisons)[1]<-'variable'
# write.xlsx(sum_of_signals_comparisons, file=paste0('../../10.Graphs_mutations/',replicate_name,'sum_vs_freq_pairwise_comparisons_',distance_cutoff,'.xlsx'))

p_results <- pairwise.wilcox.test(results$sum_of_strength, results$category2_freq, p.adjust.method = "bonferroni")
write.table(p_results$p.value, file = "Wilcoxon_p_adjusted_bonferroni_freqofmutationbins.tsv", sep = "\t", quote = F, col.names = NA, row.names = T)

results_subset <- subset(results, select = c(chr,start,end,strand,base,freq,motif,gene,category,total_hits,sum_of_strength,bins,bins_mutations,category2_freq))
# write.table(results_subset, file = "results_subset.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

stats <- results_subset%>%
  group_by(category2_freq)%>%
  summarise(Mean=mean(sum_of_strength), Median=median(sum_of_strength), Std=sd(sum_of_strength))

stats_table <- as.data.frame(matrix(nrow = nrow(stats), ncol = nrow(stats)))
for (i in 1:nrow(stats)){
  colnames(stats_table)[i] <- as.character(stats$category2_freq[i])
  for (j in 1:nrow(stats)){
    stats_table[i,j] <- (stats$Mean[i]/stats$Mean[j])
    if (i==1){
      rownames(stats_table)[j] <- as.character(stats$category2_freq[j])
    }
  }
}
write.table(stats_table, file = "FC_mean_freqofmutationbins.tsv", sep = "\t", quote = F, col.names = NA, row.names = T)

stats_table_median <- as.data.frame(matrix(nrow = nrow(stats), ncol = nrow(stats)))
for (i in 1:nrow(stats)){
  colnames(stats_table_median)[i] <- as.character(stats$category2_freq[i])
  for (j in 1:nrow(stats)){
    stats_table_median[i,j] <- (stats$Median[i]/stats$Median[j])
    if (i==1){
      rownames(stats_table_median)[j] <- as.character(stats$category2_freq[j])
    }
  }
}

write.table(stats_table_median, file = "FC_median_freqofmutationbins.tsv", sep = "\t", quote = F, col.names = NA, row.names = T)

number_of_Cs_freq <- as.data.frame(table(results$category2_freq))
number_of_Cs_RPM <- as.data.frame(table(all_mutated2$bins))
names(number_of_Cs_freq) = names(number_of_Cs_RPM) = c("Category","Count")

write.table(number_of_Cs_freq, file = "number_of_Cs_freq.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(number_of_Cs_RPM, file = "number_of_Cs_RPM.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

save.image()
