//Parameters for the coalescence simulation program : fastsimcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
north.q1
north.q3
south.q2
//Samples sizes and samples age 
36 0
34 0
36 0
//Growth rates	: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
Tdiv1 1 0 1 RESIZE1 0 0 
Tdiv2 2 0 1 RESIZE2 0 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1e-9 OUTEXP
