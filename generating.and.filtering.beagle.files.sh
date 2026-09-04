#!/bin/bash
#SBATCH --account=your-account
#SBATCH --time=0-12:00
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

REF=reference/PO1735_Microstegium_vimineum.only.23.pseudochr.fasta
#generate beagle files, i.e., snp calling in the angsd way
for chr in $(cat stiltgrass.chrs);
do angsd -bam bams.vimineum -fai $REF.fai -ref $REF -nInd 431 -minInd 130 -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -out geno.for.ld/$chr.vimineum.for.ld.maf0.1.minInd0.5 -gl 2 -minMapQ 30 -minQ 20 -minMaf 0.1 -SNP_pval 1e-6 -nThreads 16 -baq 1 -only_proper_pairs 1 -remove_bads 1 -r $chr:;
done;

#further filtering of beagle files
#get headers from beagle and mafs files of any chromosome
zcat geno.for.ld/Scaffold_10.vimineum.for.ld.maf0.1.minInd0.5.beagle.gz | head -1 > header.beagle
zcat geno.for.ld/Scaffold_10.vimineum.for.ld.maf0.1.minInd0.5.mafs.gz | head -1 > header.mafs

for chr in $(cat stiltgrass.chrs);
do sed 's/:/_/g' ld.estimation/$chr.unlinked.positions > tmp.$chr.unlinked.positions;
zcat geno.for.ld/$chr.vimineum.for.ld.maf0.1.minInd0.5.beagle.gz | awk 'NR==FNR{a[$1]; next} $1 in a' tmp.$chr.unlinked.positions - > tmp.$chr.filtered.txt;
cat header.beagle tmp.$chr.filtered.txt |bgzip -c > ld.estimation/$chr.vimineum.lowLD.maf.0.1.beagle.gz;
done;

cat header.beagle tmp*.filtered.txt | bgzip -c > ld.estimation/all.chrs.vimineum.lowLD.maf.0.1.beagle.gz
rm tmp*
