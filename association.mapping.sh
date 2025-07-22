#SBATCH --time=5-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8

REF=reference/PO1735_Microstegium_vimineum.only.23.pseudochr.fasta

#for association mapping, I have used the hybrid model (-doAsso 5) implemented in angsd (https://www.popgen.dk/angsd/index.php/Association)
#this approach first do a score test and then a latent factor model with sites that passed a threshold (-hybridThres 0.05)
for chr in $(cat stiltgrass.chromosomes);
do angsd -bam bams.case.control.all.samples -nInd 300 -minInd 100 -doMajorMinor 1 -doPost 1 -doGlf 2 -doMaf 1 -out gwas/$chr.gwas.hybrid.model \
-gl 2 -minMapQ 30 -minQ 20 -minMaf 0.05 -SNP_pval 1e-6 -r $chr: -doAsso 5 -hybridThres 0.05 -cov covariates -yBin phenotype -baq 1 -ref $REF -fai $REF.fai;
done;

#bams.case.control.all.samples is a list of samples with known phenotype: absence or presence of awns
#phenotype is a single column file with phenotypes coded as 0 (absence) or 1 (presence)
#as covariates I included the first four principal components from the PCA (low LD dataset)
#IMPORTANT: the above three files must have the same number of rows (samples) in the same order

#windows-based analysis using the weighted-zscore developed by Booker et al. 2023 (10.1111/1755-0998.13768)
#prepare the input file with the following columns: chr 	position	snp	maf	window_10kb	var1 var2 var3
#var 1 var2 are pvalues for each snp and individual variable
for var in var1 var2 var3;
do python3 general_WZA_script.py --correlations input.csv --summary_stat $var --window window_10kb --MAF maf --output $var.csv --sep "," --min_snps 5;
done;
