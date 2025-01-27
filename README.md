## Part II: Comparing shared ORFs to ANI per genome (in LGT-free simulated genomes)

---  
#### What it does: Calculates the number of orthologous genes per genome set at different pct_id thresholds.
#### Why? To plot how the detectibility of homology erodes as genomes diverge in sequence space.

#### High level steps:
1. Takes a set of simulated LGT-free genomes. (All simulations are copies of the same real genome but with different numbers of point-mutations added in silico).
2. Runs pyANI comparison of all genomes in the set to calculate Average Nucleotide Identity.
4. Uses prokka to call genes (i.e. ORFs) for each genome.
2. Creates a BLAST database of all ORFs in the genome set.
3. Runs a BLASTn query of EACH individual genome's ORFs against ALL genomes' ORFs in the set.
4. Summarizes each BLAST file per query genome, writes each summary per SAG to a file.
5. Combines all summaries into one file, combines with pyANI output.
6. Plots, like the above, to illustrate detected ORFs vs. cumulative point-mutation/GND.

#### Inputs and outputs:
* Inputs
  * This pipeline takes a zipped set of simulated genomes.
  * Expects a tarball named something like "AG-359-G18_a22.tar.gz"
  * The tarball contains genomes named like "AG-899-G06_a22_gnd001.fasta", "AG-899-G06_a22_gnd002.fasta", etc.
  * What do the names mean?
    * The rootname [AG-359-G18] is the name of a real bacterial genome from the GORG-tropics dataset.
    * The first suffix (either [a5] or [a22]) represents the alpha value used in the model that simulated point mutations.
    * The last suffix (e.g. [gnd001] vs. [gnd002]  represents the simulated genome's GND from the real genome.
* Outputs
  * For each input genome set, it produces a separate plot (see above).
  * Context: Our study examined 56 different genome sets, found here: [XXXX]. Each set contains 31 simulated genomes.
* Notes
  * Expects genome sequence names in default SPAdes format, like ">AG-359-G18_NODE_1"

---
#### Installation:



#### Test data
