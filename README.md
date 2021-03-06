# pedigree
The folder was originally to be a log of how I created the genetically fixed pedigree using the original pedigree made in the field and the genomic kinship matrix from the first call of the reference panel. It has since expanded a bit but its purpose is overall to "log and maintain files needed or produced during the making of a pedigree for different purposes". 
There are more than 1 purposes so I will try to list them here and name the files implicated in each and the process behind them if needed.  

## Genetically fixed pedigree  
Focused around the `making_ped_with_genkinship.rmd` are several other files that were used or produced during the making of the genetically fixed pedigree.   

### Input files 
- `total_pedigree_052022.csv` is the pedigree created from the barn owl database Clutch and Indiv files on May 2022. The script to do this can be also added here (*TO DO*).   
- `swiss_beta_matrix.table` is a genomic kinship matrix that only has the Swiss samples of the reference panel. The GDS needed to produce this with SNPRelate is available here due to a later modification but the link is not made in code - although I get the overall kinship matrix in `making_a_maximally_unrelated_set.rmd`.   
- `ref_panel_k0.table` & `ref_panel_k1.table` are tables of K0 and K1 kinship produced like `swiss_beta_matrix.table` using SNPRelate from the gds available here.  
- `depth_stats_all.idepth` is an individual depth file output of VCFtools on a version of the refpanel (pre-BQSR) and is used to prioritize individuals based on depth.  
- `refpanel_metadata.csv` file that has the old-new names of individuals +  the population they come from and their sex (genetic).  

### Output files 
- `genetically_fixed_pedigree.csv` is the genetically fixed pedigree produced from the RMD.  
- `individuals_to_keep.csv` is a list of individuals retains after removing the 'twins' from the list of swiss in the reference panel.  
- `succint_ped.csv` a pedigree that only includes sequenced and one more generation and still captures most links in the dataset. Used to plot the next: 
- `succint_sequenced_ped.pdf` a plot of the succing pedigree where through the `making_a_set_of_families.r` script it is split in families and replotted in `succint_ped_per_family.pdf` for easier reading. 
- `making_ped_with_genkinship.pdf` & `making_ped_with_genkinship.html` pdf and html knits of the rmd

## Making a manual pedigree for LEPMAP3

Using `succint_ped.csv` and through `making_a_set_of_families.r` we can plot `succint_ped_per_family.pdf`  and then manually select trios as input for LepMap3. The pre-conversion file in readable format (before transposing & s/,/\t/g) is given in `manual_LM_pedigree.csv`. 

## Selecting unrelated individuals 

Selecting unrelated individuals is focused around `making_a_maximally_unrelated_set.rmd` and its `.html` knitted output and uses initially the `swiss_beta_matrix.table` for the swiss individuals but after the `ref_panel_snps_f1_masked_maf05_miss05_LDpruned.gds` for the EU individuals. It also outputs the `unrelated_0335_kinship.csv` which includes a set of 76 swiss individuals pruned at 0.0335 kinship.   
At some point I decided that we need a different cutoff of unrelated individuals and added unrelated_001_kinship_42_040722.csv and a pdf with a couple of plots for the beta of these guys. They show promising results in smc++.   

## WhatsHap phasing pedigree 

This pedigree was created by Sonia and includes trios to be used when phasing individuals in pedigree mode of WhatsHap ... it is in core similar to LepMap3 `manual_LM_pedigree.csv` but also includes small trios that the former does not (since a single trio is sufficient for read-based phasing but not for recombination inference).   *TO ADD FILE FROM SONIA*

## Sexing the 3k owls 

Tristan produced the sex of lcWGS individuals in `CoverageLowCov_Samples.txt` by comparing the coverage of an autosomal scaffold and a Z scaffold. Then in the `match_3kowls_sex_and_seq.r` script I match the sex of genomic data with the pedigree available sex to identify mistakes. I also use the `Metadata 3Kowls - 3kOwls_LowCoverage.csv` to match old-new names and output the plot in `sexxx_diff_3k.pdf`.

