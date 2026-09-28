#!/bin/bash
#SBATCH --account=your-account
#SBATCH --time=1-0:00
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=40

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

#create 100 folders with the required input files
for i in {1..100}
do mkdir run$i
cp north.q3_DAFpop0.obs run$i/
cp north.q3.tpl run$i/
cp north.q3.est run$i/
done


#run fastsimcoal, example for north.q3 group
#!/bin/bash
for i in {1..100}
do cd run$i
fsc27 -t north.q3.tpl -n 100000 -e north.q3.est -M -L 100 -c 40 -B 40 -d -q
cd ..
done
