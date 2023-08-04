# Cheng_2023_vitC
Scripts and files used in Cheng et al "Vitamin C activates young LINE-1 elements in mouse embryonic stem cells via H3K9me3 demethylation", bioRxiv


## ChIP

*ChIP-qPCR*

This folder contains H3K9me2 and H3K9me3 ChIP-qPCR data and associated script to generate plots.

*H3K9me2 ChIP-seq*

Data from Ebata et al (PMID: 28706564) were processed in Seqmonk to extract read counts (RPM) from full-length (>5kb) young L1s, as well 5kb-long regions downstream of the same elements. The associated script calculates the enrichment ratio between L1 and downstream regions and produces a boxplot.

*KDM4 ChIP-seq*

Data from Pedersen et al (PMID: 27266524, mouse KDM4A), Tomaz et al (PMID: 28087629, mouse KDM4C) and ENCODE (human KDM4A) were mapped allowing for multi-mapping reads, to generate average profiles over full-length young L1s. Scripts 'Kdm4_L1_profile.R' and 'hKDM4A_L1_profile.R' plot these profiles. 'Kdm4a_L1_subclasses.R' plots data on the number of mouse L1 elements overlapping with KDM4A peaks (either including multimapping reads, or just unique ones).


## oxBS

This folder contains all processed oxBS amplicon-seq data and associated script to generate plots.


## RNAseq

*mESC Genes*

Raw reads counts over genes were processed by DESeq2 ('DESeq_genes.R' script) to generate normalised gene expression values ('gene_expression_vsd.txt') and a list of differentially expressed genes ('de_genes_vitC.txt'). The same script generates an MA plot.

*mESC TEs*

RNA-seq data was processed using SQuIRE. The 'subFcounts' folder contains TE subfamily-level data, which was analysed using 'DESeq_subfamily.R' to identify differentially expressed subfamilies ('de_squire_subfamily.txt') and generate an MA plot. The 'flagSubset' folder contains data on individual elements, which was analysed using: 1) 'L1_flag_data.R' to plot the expression of independently transcribed L1 elements, 2) 'ETnERV3_flag_data.R' to identify upregulated proviral ETnERV3 elements ('up_ETnERV3.txt').

*hESC L1*

This folder contains data processed by TEtranscripts, TExP and SQuIRE for the analysis of L1 subfamily expression. Data from SQuIRE on individual elements is also included. All associated plots are generated using the 'hESC_L1.R' script. A file containing hKDM4A peaks ('ENCFF021QGZ.bed' from ENCODE) is also included to check for potential associations with expression changes.


## RT-qPCR

Multiple folders with processed data and sample info for each experiment are included. Nearly all plots using these data are generated using 'qPCR_plots.R'. There is a separate script for the time-course experiment ('qPCR_timeCourse.R'), and one for the hESC data in Supplementary Figure 4A ('qPCR_hESC2.R').


## TEbreak

This folder contains L1 insertion mapping data processed by TEbreak for both 13-day and 41-day experiments (the latter with two different read cut-offs). The 'TEBreak_summary.R' script identifies clonal insertions (output files in folder) as well as somatic insertions unique to one of the branches of the experiment.


