#!/bin/bash
#SBATCH --account=your_compute_canada_account
#SBATCH --time=3-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

REF=reference/PO1735_Microstegium_vimineum.only.23.pseudochr.fasta

#Apparently the process is simpler by using the option -rmSNPs 1, however, this option didn't work when I ran these analyses
#The first step is calling snps present in the outgroup samples, because some SNPs segregate between species, only monomorphic sites in the outgroups can be considered ancestral
angsd -b list.outgroups -out outgroup.snps -minMapQ 30 -minQ 20 -baq 1 -ref $REF -doCounts 1 -snp_pval 1e-6 -domaf 1 -domajorminor 1 -minind 10 -P 8 -gl 2
#The second step is creating a reference file, which I name ancestral.ref.with.snps because both variant and invariant sites are included, make sure to use -explode 1 to include all positions in the reference file
angsd -b list.outgroups -doFasta 2 -out ancestral.ref.with.snps -minMapQ 30 -minQ 20 -baq 1 -ref $REF -doCounts 1 -domaf 1 -domajorminor 1 -gl 2 -minind 10 -explode 1 -P 8
#The next step is using the snps called in the first step to mask our ancestral reference, so only invariant sites in the outgroup are included
zcat outgroup.snps.mafs.gz |awk '{print$1"\t"$2}'|awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $2+1}' > snps.outgroup.bed
#I removed manually the header of the bed file
#Masking the reference
bedtools maskfasta -fi ancestral.ref.with.snps.fasta -bed fixed.bed -fo ancestral.reference.fasta
#indexing the ancestral reference
bwa index ancestral.reference.fasta
samtools faidx ancestral.reference.fasta
#Create a folder where the ancestral reference files will be
mkdir ancestral.ref
mv ancestral.reference.fasta* ancestral.ref/

#The ancestral reference will be used to polarize the snps and estimate the unfolded sfs, can be also used to compute abba-baba tests
