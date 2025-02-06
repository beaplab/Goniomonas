##############################################
######## DATASET PREPARATION #################
######## (scroll down for tree building) #####
##############################################



#manually removed other headings
##############################
#Here is where EukRibo dataset is written in file Cryptophyceae_EukRibo_full_seqs_clean.fasta is ready to go
##############################



#put seqs from various sources together
cat ./240_Cryptista_EukRibo_NCBI_PacBio.fas ./Cryptophyceae_EukRibo_full_seqs_clean.fasta ./all_Goniomonas.fasta > all_gonio_Jamy_EukRibo_NCBI_1.fasta

#mafft -1 
"/usr/bin/mafft"  --auto --inputorder "all_gonio_Jamy_EukRibo_NCBI.fasta" > "all_gonio_Jamy_EukRibo_NCBI_mafft.fasta" 


#visualise alignment and delete redundant sequences from the outgroup 
#calculate NJ tree via Jalview (by PID) saved to NJ_for_draft_mafft.nwk
#Manually deleteing closely related sequences from outgroup; saved into all_gonio_Jamy_EukRibo_NCBI_mafft_outgroup_edit.fasta


# {{{{ seqtk subseq ./all_gonio_Jamy_EukRibo_NCBI_mafft.fasta ./IDs_to_keep.txt > ./all_gonio_Jamy_EukRibo_NCBI_mafft_outgroup_edit.fasta



####################################
############Eukref run
####################################


# Sequences from goniomonas from all_gonio_Jamy_EukRibo_NCBI_mafft_outgroup_edit.fasta were selected (i.e.) no outgroup included, gaps deleted;
#resulting file saveed to goniomonads_nogaps.fasta:

#in /home/dzavadska/Data/Goniomonas/TO_UPLOAD/TREES/dataset_formation/eukref
# seqtk subseq ./all_gonio_Jamy_EukRibo_NCBI_mafft_outgroup_edit.fasta goniomonad_IDs_to_keep.fasta > ./goniomonads.fasta


sed 's/\-//g' goniomonads.fasta > goniomonads_nogaps.fasta
sed 's/\_.*//' goniomonads_nogaps.fasta > goniomonads_nogaps_renamed.fasta

#on nisaba

/home/dzavadska/eukref_datasets/usearch -sortbylength goniomonads_nogaps_renamed.fasta -fastaout goniomonads_nogaps_renamed.sorted.fasta -minseqlength 500 -notrunclabels
/home/dzavadska/eukref_datasets/usearch -cluster_smallmem goniomonads_nogaps_renamed.sorted.fasta -id 0.97 -centroids goniomonads_nogaps_renamed.clustered.fasta -uc goniomonads_nogaps_renamed.clusters -notrunclabels

nohup sudo python2 eukref_gbretrieve.py -i ../Goniomonas_eukref_FIN/goniomonads_nogaps_renamed.clustered.fasta -dbnt /scratch/data1/aauladell/databases/nt/nt_blastdb/nt_blast -dbsi Reference_DB.fas -n 100 -p 32 -g Goniomonadales -m megablast -idsi 75 -idnt 80 -td tax_d.bin

#new putatively goniomonas sequences not present in the input alignment
# EF526832 ==> CRY-4 lineage ALREADY INCLUDED INTO THE OUTGROUP
# EF526920  ==> sister to CRY-4 lineage ==> added to alignment
# EF674349 ==> Goniomonas related to LC683686 clade ==> added to alignment


# eukref-recovered sequences were added to ./all_gonio_Jamy_EukRibo_NCBI_mafft_outgroup_edit.fasta and resulting file is called eukref_goniomonads_outgroup.fasta

seqkit rmdup -D ./duplicates.fasta -s ./eukref_goniomonads_outgroup.fasta  > ./eukref_goniomonads_outgroup_nodup.fasta
#manually removed shorter duplicate of BEAP0078
#0075 eliminated as duplicate

#MAFFT v7.453 (2019/Nov/8)
"/usr/bin/mafft"  --auto --inputorder "eukref_goniomonads_outgroup_nodup.fasta" > "eukref_goniomonads_outgroup_mafft.fasta"
fasttree -nt ./eukref_goniomonads_outgroup_mafft.fasta > ./eukref_goniomonads_outgroup_mafft_fasttree

############################################################
######## ALIGNMENT AND TREE-BUILDING ITSELF ################
######## tree with big outgroup - used in supplement #######
############################################################


cp ./eukref_goniomonads_outgroup_nodup.fasta ./FIN.fasta

#manually added BEAP0271, BEAP0296, BEAP0351
seqkit rmdup -D ./duplicates.fasta -s ./FIN.fasta  > ./FIN_nodup.fasta
#no duplicates recovered at this point, continue with FIN.fasta
"/usr/bin/mafft"  --localpair  --maxiterate 1000 --inputorder "FIN.fasta" > "FIN_mafft.fasta" 

#trying a variety of trimming options
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt05.fasta -gt 0.5 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt02.fasta -gt 0.2 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt01.fasta -gt 0.1 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout.fasta -gappyout -fasta

#detecting redundant sequences which appear as a result of trimming
seqkit rmdup -D ./duplicates_gt05.fasta -s ./FIN_mafft_gappyout_gt05.fasta #4 duplicates
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN_mafft_gappyout_gt02.fasta #4 duplicates
seqkit rmdup -D ./duplicates_gt01.fasta -s ./FIN_mafft_gappyout_gt01.fasta #3 duplicates
seqkit rmdup -D ./duplicates_gappyout.fasta -s ./FIN_mafft_gappyout.fasta #13 duplicates

#detecting parsimony informative sites; col1: number of parsimony informative sites col2: total number of sites col3: percentage of parsimony informative sites
phykit parsimony_informative_sites ./FIN_mafft_gappyout_gt05.fasta #723	1703	42.4545
phykit parsimony_informative_sites ./FIN_mafft_gappyout_gt02.fasta #759	1781	42.6165
phykit parsimony_informative_sites ./FIN_mafft_gappyout_gt01.fasta #880	1958	44.9438
phykit parsimony_informative_sites ./FIN_mafft_gappyout.fasta #530	1320	40.1515

#Based on manual inspection and parsimony informative sites metrics, we continue with alignment FIN_mafft_gappyout_gt02.fasta 
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN_mafft_gappyout_gt02.fasta > ./FIN_mafft_gt02_FOR_TREE.fasta
#[INFO] 4 duplicated records removed : 2	EF526832.1/1-1390, EF526832_Cryptophyceae|CRY4-lineage|clone=NA1_3C10/1-1390
#2	SL4-0466_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-418_conseq_Otu466_6_soil_Skogaryd-peat/1-1424, FW3-0381_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-1293_conseq_Otu0381_19_freshwater_permafrost/1-1423
#2	FW3-0449_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-1344_conseq_Otu0449_14_freshwater_permafrost/1-1423, FW5-0385_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-305_conseq_Otu385_3_freshwater_Svartberget-fw/1-1422
#2	FW3-0548_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-1279_conseq_Otu0548_9_freshwater_permafrost/1-1420, FW5-0034_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-323_conseq_Otu034_913_freshwater_Svartberget-fw/1-1420



/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN_mafft_gt02_FOR_TREE.fasta --datatype nt --template raxml

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            GTR+I+G4    67838.6414        1.0000
#       AIC            GTR+I+G4    65726.9432        1.0000
#      AICc            GTR+I+G4    65939.9432        1.0000


raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -# 20 -s FIN_mafft_gt02_FOR_TREE.fasta -n FIN_23Oct2024
raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -b 12345 -# 1000 -s FIN_mafft_gt02_FOR_TREE.fasta -n FIN_23Oct2024_BS
raxmlHPC -m GTRGAMMAIX -p 12345 -f b -t RAxML_bestTree.FIN_23Oct2024 -z RAxML_bootstrap.FIN_23Oct2024_BS -n FIN_23Oct2024_BS_BIPART


/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN_mafft_gappyout_gt02.fasta --datatype nt -o mb --template mrbayes

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            GTR+I+G4    67861.1129        0.9917
#       AIC            GTR+I+G4    65705.5353        1.0000
#      AICc            GTR+I+G4    65928.5353        1.0000


#Preformatting for fkng MrBayes
#cp FIN_MAFFT_GONIO.fasta ./FIN_MAFFT_mb_rename.fasta


#manually rename HFCC_1666 into HFCC1666 otherwise gets remaned not nicely 
# convert fasts to nexus by running trimal again
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt02.nex -gt 0.2 -nexus

#cp ./FIN_mafft_gappyout.nex ./MrBayes_3/

#manually modify datatype to DNA
sed -i 's@/@VVV@g' FIN_mafft_gappyout_gt02.nex
sed -i 's@_.*.1-@@g' FIN_mafft_gappyout_gt02.nex
#manually remove duplicates, opening the .nex alignment in Aliview


mb
#mbprompt follows
#the number of gamma categories defaults to 4 in mrbayes
execute FIN_mafft_gappyout_gt02.nex
lset nst=6 rates=invgamma
mcmc ngen=2000000 samplefreq=10
sump
sumt

#stopped at generation 3010000, with Average standard deviation of split frequencies: 0.011541

#manually create csv file with two columns containing outputs from these two commands (dictionary of names)
grep ">" FIN_mafft_gappyout_gt02.fasta
grep ">" FIN_mafft_gappyout_gt02.fasta | sed 's@/@VVV@g' | sed 's@_.*.@@g' | sed 's@VVV.*.@@g' | sed 's@ 1781 bp@@g' | sed 's@>@@g'


#Run the tree renamer R script, and with the output from it:

#renaming tree outputs
sed 's@special_Separator.*@@g' tree_renamed_custom_mb.nwk > tree_renamed2_custom_mb.nwk
sed 's@VVV.*@@g' tree_renamed2_custom_mb.nwk > tree_renamed3_custom_mb.nwk



##################################################################
######## ALIGNMENT AND TREE-BUILDING ITSELF ######################
######## tree with ONLY CRY-1 outgroup - used in MAIN TEXT #######
##################################################################

# FIN_2.fasta is manually produced from FIN.fasta by removing all the outgroup sequences apart from the ones belonging to CRY-1.



"/usr/bin/mafft"  --localpair  --maxiterate 1000 --inputorder "FIN_2.fasta" > "FIN2_mafft.fasta" 


#trying a variety of trimming options
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt05.fasta -gt 0.5 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt02.fasta -gt 0.2 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt01.fasta -gt 0.1 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout.fasta -gappyout -fasta


#detecting redundant sequences which appear as a result of trimming
seqkit rmdup -D ./duplicates_gt05.fasta -s ./FIN2_mafft_gappyout_gt05.fasta #3 duplicates
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN2_mafft_gappyout_gt02.fasta #2 duplicates
seqkit rmdup -D ./duplicates_gt01.fasta -s ./FIN2_mafft_gappyout_gt01.fasta #1 duplicates
seqkit rmdup -D ./duplicates_gappyout.fasta -s ./FIN2_mafft_gappyout.fasta #11 duplicates


for i in ./duplicates* ; do echo "$i" ; wc -l "$i" ; done

#detecting parsimony informative sites; col1: number of parsimony informative sites col2: total number of sites col3: percentage of parsimony informative sites
phykit parsimony_informative_sites ./FIN2_mafft_gappyout_gt05.fasta #500	1717	29.1206
phykit parsimony_informative_sites ./FIN2_mafft_gappyout_gt02.fasta #629	1943	32.3726
phykit parsimony_informative_sites ./FIN2_mafft_gappyout_gt01.fasta #710	2146	33.0848
phykit parsimony_informative_sites ./FIN2_mafft_gappyout.fasta #309	1099	28.1165

#Based on manual inspection and parsimony informative sites metrics, we continue with alignment FIN2_mafft_gappyout_gt02.fasta 
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN2_mafft_gappyout_gt02.fasta > ./FIN2_mafft_gt02_FOR_TREE.fasta

/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN2_mafft_gt02_FOR_TREE.fasta --datatype nt --template raxml

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            HKY+I+G4    31561.7558        0.6294
#       AIC            GTR+I+G4    30442.8456        1.0000
#      AICc            GTR+I+G4    30488.8456        0.9999



raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -# 20 -s FIN2_mafft_gt02_FOR_TREE.fasta -n FIN_9Dec2024
raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -b 12345 -# 1000 -s FIN2_mafft_gt02_FOR_TREE.fasta -n FIN_9Dec2024_BS
raxmlHPC -m GTRGAMMAIX -p 12345 -f b -t RAxML_bestTree.FIN_9Dec2024 -z RAxML_bootstrap.FIN_9Dec2024_BS -n FIN_9Dec2024_BS_BIPART


#mb

/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN2_mafft_gappyout_gt02.fasta --datatype nt -o mb --template mrbayes

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            HKY+I+G4    31582.6340        0.6438
#       AIC            GTR+I+G4    30441.5867        1.0000
#      AICc            GTR+I+G4    30489.5867        0.9999



#Preformatting for fkng MrBayes
#cp FIN_MAFFT_GONIO.fasta ./FIN_MAFFT_mb_rename.fasta


#manually rename HFCC_1666 into HFCC1666 otherwise gets remaned not nicely 
# convert fasts to nexus by running trimal again
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt02.nex -gt 0.2 -nexus

#cp ./FIN_mafft_gappyout.nex ./MrBayes_3/

#manually modify datatype to DNA
sed -i 's@/@VVV@g' FIN2_mafft_gappyout_gt02.nex
sed -i 's@_.*.1-@@g' FIN2_mafft_gappyout_gt02.nex
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#manually remove duplicates, opening the .nex alignment in Aliview
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#for GTR

mb
#mbprompt follows
#the number of gamma categories defaults to 4 in mrbayes
execute FIN2_mafft_gappyout_gt02.nex
lset nst=6 rates=invgamma
mcmc ngen=1000000 samplefreq=10
sump
sumt

#Average standard deviation of split frequencies: 0.010011

#for HKY

mb
#mbprompt follows
#the number of gamma categories defaults to 4 in mrbayes
execute FIN2_mafft_gappyout_gt02.nex
lset nst=2 rates=invgamma
mcmc ngen=1000000 samplefreq=10
sump
sumt



#Run the tree renamer R script, and with the output from it:

#renaming tree outputs
sed 's@special_Separator.*@@g' gtr_tree2_0_renamed_custom_mb.nwk > gtr_tree2_0_renamed2_custom_mb.nex
sed 's@VVV.*@@g' gtr_tree2_0_renamed2_custom_mb.nex > gtr_tree2_0_renamed3_custom_mb.nex

sed 's@special_Separator.*@@g' hky_tree2_0_renamed_custom_mb.nwk > hky_tree2_0_renamed2_custom_mb.nex
sed 's@VVV.*@@g' hky_tree2_0_renamed2_custom_mb.nex > hky_tree2_0_renamed3_custom_mb.nex




