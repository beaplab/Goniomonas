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


#manually remove BEAP sequences from the same strains (they are shorter versions of original sequences, and have letters like a b c appended to bEAP ID) ======> eukref_goniomonads_outgroup_mafft_2_0.fasta 



#add sequences from phanphrasert et. al. - fetching and merging with previous file
for i in $(cat phannph_ids.txt) ; do efetch -db nuccore -id "$i" -format fasta >> phannph_seqs.fasta ; done

cat phannph_seqs.fasta eukref_goniomonads_outgroup_mafft_2_0.fasta >> eukref_goniomonads_outgroup_mafft_3_0.fasta 

############################################################
######## ALIGNMENT AND TREE-BUILDING ITSELF ################
######## tree with big outgroup - used in supplement #######
############################################################


cp ./eukref_goniomonads_outgroup_mafft_3_0.fasta ./FIN.fasta

#manually added BEAP0271, BEAP0296, BEAP0351
# seqkit rmdup -D ./duplicates.fasta -s ./FIN.fasta  > ./FIN_nodup.fasta
#no duplicates recovered at this point, continue with FIN.fasta
"/usr/bin/mafft" --localpair --maxiterate 1000 --inputorder "FIN.fasta" > "FIN_mafft.fasta" 

#trying a variety of trimming options
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt05.fasta -gt 0.5 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt02.fasta -gt 0.2 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt01.fasta -gt 0.1 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout.fasta -gappyout -fasta

#detecting redundant sequences which appear as a result of trimming
seqkit rmdup -D ./duplicates_gt05.fasta -s ./FIN_mafft_gappyout_gt05.fasta #6 duplicates
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN_mafft_gappyout_gt02.fasta #4 duplicates
seqkit rmdup -D ./duplicates_gt01.fasta -s ./FIN_mafft_gappyout_gt01.fasta #2 duplicates
seqkit rmdup -D ./duplicates_gappyout.fasta -s ./FIN_mafft_gappyout.fasta #19 duplicates

for i in ./duplicates* ; do echo "$i" ; wc -l "$i" ; done

#detecting parsimony informative sites; col1: number of parsimony informative sites col2: total number of sites col3: percentage of parsimony informative sites
phykit parsimony_informative_sites ./FIN_mafft_gappyout_gt05.fasta #725	1708	42.4473
phykit parsimony_informative_sites ./FIN_mafft_gappyout_gt02.fasta #755	1778	42.4634
phykit parsimony_informative_sites ./FIN_mafft_gappyout_gt01.fasta #903	1976	45.6984
phykit parsimony_informative_sites ./FIN_mafft_gappyout.fasta #434	1102	39.3829

#Based on manual inspection and parsimony informative sites metrics, we continue with alignment FIN_mafft_gappyout_gt02.fasta 
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN_mafft_gappyout_gt02.fasta > ./FIN_mafft_gt02_FOR_TREE.fasta
# FW3-0381_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-1293_conseq_Otu0381_19_freshwater_permafrost/1-1423
# FW5-0034_Cryptophyceae|goniomonads|g-Goniomonas|Gtrun-group|clone=c-323_conseq_Otu034_913_freshwater_Svartberget-fw/1-1420

/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN_mafft_gt02_FOR_TREE.fasta --datatype nt --template raxml

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            GTR+I+G4    69239.8222        1.0000
#       AIC            GTR+I+G4    66854.6108        1.0000
#      AICc            GTR+I+G4    67136.6108        1.0000



raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -# 20 -s FIN_mafft_gt02_FOR_TREE.fasta -n FIN_18Oct2025
raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -b 12345 -# 1000 -s FIN_mafft_gt02_FOR_TREE.fasta -n FIN_18Oct2025_BS
raxmlHPC -m GTRGAMMAIX -p 12345 -f b -t RAxML_bestTree.FIN_18Oct2025 -z RAxML_bootstrap.FIN_18Oct2025_BS -n FIN_18Oct2025_BS_BIPART


/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN_mafft_gappyout_gt02.fasta --datatype nt -o mb --template mrbayes

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            GTR+I+G4    69170.9355        0.9939
#       AIC            GTR+I+G4    66741.8582        1.0000
#      AICc            GTR+I+G4    67035.8582        1.0000



#Preformatting for MrBayes
#cp FIN_MAFFT_GONIO.fasta ./FIN_MAFFT_mb_rename.fasta


#manually rename HFCC_1666 into HFCC1666 otherwise gets remaned not nicely 
# convert fasts to nexus by running trimal again
/home/dzavadska/trimAl/source/trimal -in FIN_mafft.fasta -out FIN_mafft_gappyout_gt02.nex -gt 0.2 -nexus

#cp ./FIN_mafft_gappyout.nex ./MrBayes_3/

#manually modify datatype to DNA
sed -i 's@/@VVV@g' FIN_mafft_gappyout_gt02.nex
sed -i 's@_.*.1-@@g' FIN_mafft_gappyout_gt02.nex
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#manually remove duplicates, opening the .nex alignment in Aliview
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

mb
#mbprompt follows
#the number of gamma categories defaults to 4 in mrbayes
execute FIN_mafft_gappyout_gt02.nex
lset nst=6 rates=invgamma
mcmc ngen=2000000 samplefreq=10
sump
sumt

#   2000000 -- (-33256.857) (-33059.576) (-33097.499) [-33058.986] * [-33031.983] (-33116.544) (-33105.407) (-33080.768) -- 0:00:00
#   Average standard deviation of split frequencies: 0.014344



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

#check for duplicates: there should be none
#seqkit rmdup -D ./duplicates.fasta -s ./FIN_2.fasta  > ./FIN_2nodup.fasta

"/usr/bin/mafft" --localpair --maxiterate 1000 --inputorder "FIN_2.fasta" > "FIN2_mafft.fasta" 


#trying a variety of trimming options
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt05.fasta -gt 0.5 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt02.fasta -gt 0.2 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt01.fasta -gt 0.1 -fasta
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout.fasta -gappyout -fasta


#detecting redundant sequences which appear as a result of trimming
seqkit rmdup -D ./duplicates_gt05.fasta -s ./FIN2_mafft_gappyout_gt05.fasta #5 duplicates
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN2_mafft_gappyout_gt02.fasta #2 duplicates
seqkit rmdup -D ./duplicates_gt01.fasta -s ./FIN2_mafft_gappyout_gt01.fasta #1 duplicates
seqkit rmdup -D ./duplicates_gappyout.fasta -s ./FIN2_mafft_gappyout.fasta #16 duplicates


for i in ./duplicates* ; do echo "$i" ; wc -l "$i" ; done

#detecting parsimony informative sites; col1: number of parsimony informative sites col2: total number of sites col3: percentage of parsimony informative sites
phykit parsimony_informative_sites ./FIN2_mafft_gappyout_gt05.fasta #502	1717	29.237
phykit parsimony_informative_sites ./FIN2_mafft_gappyout_gt02.fasta #620	1913	32.4098
phykit parsimony_informative_sites ./FIN2_mafft_gappyout_gt01.fasta #712	2133	33.3802
phykit parsimony_informative_sites ./FIN2_mafft_gappyout.fasta #314	1101	28.5195

#Based on manual inspection and parsimony informative sites metrics, we continue with alignment FIN2_mafft_gappyout_gt02.fasta 
seqkit rmdup -D ./duplicates_gt02.fasta -s ./FIN2_mafft_gappyout_gt02.fasta > ./FIN2_mafft_gt02_FOR_TREE.fasta

/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN2_mafft_gt02_FOR_TREE.fasta --datatype nt --template raxml

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            HKY+I+G4    32222.6910        0.9030
#       AIC            GTR+I+G4    30844.7765        0.9999
#      AICc            GTR+I+G4    30919.7765        0.9994


raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -# 20 -s FIN2_mafft_gt02_FOR_TREE.fasta -n FIN2_18Oct2025
raxmlHPC -m GTRGAMMAIX -p 12345 -T 8 -b 12345 -# 1000 -s FIN2_mafft_gt02_FOR_TREE.fasta -n FIN2_18Oct2025_BS
raxmlHPC -m GTRGAMMAIX -p 12345 -f b -t RAxML_bestTree.FIN2_18Oct2025 -z RAxML_bootstrap.FIN2_18Oct2025_BS -n FIN2_18Oct2025_BS_BIPART


#mb

/home/dzavadska/Data/soft/modeltest/bin/modeltest-ng --input FIN2_mafft_gappyout_gt02.fasta --datatype nt -o mb --template mrbayes

#                         Model         Score        Weight
#----------------------------------------------------------
#       BIC            GTR+I+G4    32411.0272        0.7945
#       AIC            GTR+I+G4    30994.1380        1.0000
#      AICc            GTR+I+G4    31072.1380        0.9999




#Preformatting for  MrBayes
#cp FIN_MAFFT_GONIO.fasta ./FIN_MAFFT_mb_rename.fasta


#manually rename HFCC_1666 into HFCC1666 otherwise gets remaned not nicely 
# convert fasts to nexus by running trimal again
/home/dzavadska/trimAl/source/trimal -in FIN2_mafft.fasta -out FIN2_mafft_gappyout_gt02.nex -gt 0.2 -nexus


#cp ./FIN2_mafft_gappyout_gt02.nex ../mrbayes/

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
mcmc ngen=1500000 samplefreq=10
sump
sumt
  
# 1500000 -- (-15414.377) (-15407.570) (-15453.851) [-15381.736] * (-15418.645) (-15454.415) (-15405.726) [-15394.014] -- 0:00:00
#   Average standard deviation of split frequencies: 0.007153

#for HKY

mb
#mbprompt follows
#the number of gamma categories defaults to 4 in mrbayes
execute FIN2_mafft_gappyout_gt02.nex
lset nst=2 rates=invgamma
mcmc ngen=1000000 samplefreq=10
sump
sumt

#   1500000 -- (-15423.521) (-15469.211) (-15408.459) [-15386.869] * (-15429.307) (-15435.184) [-15409.435] (-15434.023) -- 0:00:00
#   Average standard deviation of split frequencies: 0.008679



#Run the tree renamer R script, and with the output from it:

#renaming tree outputs
sed -i 's@special_Separator.*@@g' tree2_0_renamed_custom_mb.nwk
sed -i 's@VVV.*@@g' tree2_0_renamed_custom_mb.nwk

sed 's@special_Separator.*@@g' hky_tree2_0_renamed_custom_mb.nwk > hky_tree2_0_renamed2_custom_mb.nex
sed 's@VVV.*@@g' hky_tree2_0_renamed2_custom_mb.nex > hky_tree2_0_renamed3_custom_mb.nex



