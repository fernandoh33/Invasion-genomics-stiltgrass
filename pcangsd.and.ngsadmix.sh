#!/bin/bash
#SBATCH --account=your-account
#SBATCH --time=2-16:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

BEAGLE=ld.estimation/all.chrs.vimineum.lowLD.maf.0.1.beagle.gz
REF=reference/PO1735_Microstegium_vimineum.only.23.pseudochr.fasta

#estimate pca using pcangsd
pcangsd -b $BEAGLE -t 16 -o out.pcangsd
#in R
samples.info=read.table("sample.info.txt",header=F) #sample.info.txt is a file with the name of the samples in a single column, with no header
out.pca=as.matrix(read.table("~/your_local_folder/out.pcangsd.cov"))
eigenvalues=eigen(out.pca)
plot(eigenvalues$vectors)
out.pca=cbind(samples.info, eigenvalues$vectors[,1:10])
colnames(out.pca)=c("sample","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

#ngsadmix
#generating 10 independent runs per K, needed to run Evanno test
for seed in 1234 5678 9101 1213 1415 1617 1819 2021 2223 2425; 
do for k in 1 2 3 4 5 6 7 8 9 10;
do NGSadmix -likes $BEAGLE -K $k -P 8 -seed $seed -o admixture/$seed.out.ngsadmix.k$k;
done;
done;
