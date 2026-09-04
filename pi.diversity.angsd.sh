#!/bin/bash
#SBATCH --account=your-account
#SBATCH --time=0-18:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

REF=reference/PO1735_Microstegium_vimineum.only.23.pseudochr.fasta

#this code need to be run for each genetic group you want to get estimates of genetic diversity, here is an example for herbarium samples collected after 1960 in the north group 
for chr in $(cat all.chrs);
do angsd -bam samples.north.after.1960 -nInd 25 -noTrans 0 -minInd 8 -doSaf 1 -baq 1 -anc $REF -ref $REF -GL 2 -P 8 -minMapQ 30 -minQ 20 -uniqueOnly 1 -only_proper_pairs 1 -r $chr -out pi.north/north.after.$chr.w.invariants.0.3;
realSFS pi.north/north.after.$chr.w.invariants.0.3.saf.idx -P 8 -fold 1 > pi.north/north.after.$chr.sfs;
realSFS saf2theta pi.north/north.after.$chr.w.invariants.0.3.saf.idx -sfs pi.north/north.after.$chr.sfs -outname pi.north/north.after.$chr.w.invariants.0.3 -fold 1;
thetaStat do_stat pi.north/north.after.$chr.w.invariants.0.3.thetas.idx -win 100000 -step 100000 -outnames pi.north/north.after.$chr.theta.w100;
thetaStat do_stat pi.north/north.after.$chr.w.invariants.0.3.thetas.idx -win 10000 -step 10000 -outnames pi.north/north.after.$chr.theta.w10;
thetaStat do_stat pi.north/north.after.$chr.w.invariants.0.3.thetas.idx -outnames pi.north/north.after.$chr.theta.global;
done;

