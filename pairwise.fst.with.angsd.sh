#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-18:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

REF=reference/PO1735_Microstegium_vimineum.only.23.pseudochr.fasta
ANC=ancestral.ref/ancestral.reference.fasta
OUTDIR=pair.fst.only.snps
filters='-doMajorMinor 4 -doPost 1 -doGlf 2 -doMaf 1 -baq 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minMaf 0.01'

for chr in $(cat all.chrs);
do angsd -bam samples.north.q1 -nInd 18 -minInd 6 $filters -doSaf 1 -anc $REF -ref $REF -GL 2 -P 8 -r $chr -out $OUTDIR/north.modern.q1.$chr.only.snps.maf.0.01;
realSFS $OUTDIR/north.modern.q1.$chr.only.snps.maf.0.01.saf.idx -P 8 -fold 1 > $OUTDIR/north.modern.q1.$chr.only.snps.maf.0.01.sfs;
done;

for chr in $(cat all.chrs);
do angsd -bam samples.north.q3 -nInd 17 -minInd 6 $filters -doSaf 1 -anc $REF -ref $REF -GL 2 -P 8 -r $chr -out $OUTDIR/north.modern.q3.$chr.only.snps.maf.0.01;
realSFS $OUTDIR/north.modern.q3.$chr.only.snps.maf.0.01.saf.idx -P 8 -fold 1 > $OUTDIR/north.modern.q3.$chr.only.snps.maf.0.01.sfs;
done;

for chr in $(cat all.chrs);
do realSFS $OUTDIR/north.modern.q1.$chr.only.snps.maf.0.01.saf.idx $OUTDIR/north.modern.q3.$chr.only.snps.maf.0.01.saf.idx -P 8 -r $chr -fold 1 > $OUTDIR/$chr.north.modern.q1.vs.north.modern.q3.ml;
realSFS fst index $OUTDIR/north.modern.q1.$chr.only.snps.maf.0.01.saf.idx $OUTDIR/north.modern.q3.$chr.only.snps.maf.0.01.saf.idx -fold 1 -r $chr -sfs $OUTDIR/$chr.north.modern.q1.vs.north.modern.q3.ml -fstout $OUTDIR/$chr.north.modern.q1.vs.north.modern.q3;
realSFS fst stats2 $OUTDIR/$chr.north.modern.q1.vs.north.modern.q3.fst.idx -win 100000 -step 100000 > $OUTDIR/$chr.fst.100kb.only.snps.north.modern.q1.vs.north.modern.q3;
realSFS fst stats2 $OUTDIR/$chr.north.modern.q1.vs.north.modern.q3.fst.idx -win 10000 -step 10000 > $OUTDIR/$chr.fst.10kb.only.snps.north.modern.q1.vs.north.modern.q3;
realSFS fst stats $OUTDIR/$chr.north.modern.q1.vs.north.modern.q3.fst.idx > $OUTDIR/$chr.fst.global.only.snps.north.modern.q1.vs.north.modern.q3;
done;
