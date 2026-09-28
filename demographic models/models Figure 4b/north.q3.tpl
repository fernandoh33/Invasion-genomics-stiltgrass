//Parameters for the coalescence simulation program : fastsimcoal.exe
1 samples to simulate :
//Population effective sizes (number of genes)
npresent
//Samples sizes and samples age 
34 0
//Growth rates	: negative growth implies population expansion
GR1
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
Tbot 0 0 0 RESIZE1 1 0 
Tendbot 0 0 0 RESIZE2 1 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1e-9 OUTEXP
