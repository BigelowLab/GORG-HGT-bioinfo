#!/user/bin/nextflow
nextflow.enable.dsl=2

/*
I. Inputs:

Expects an input folder called ./input/ that contains several tarfiles (e.g. AG-359-G18_a5.tar.gz, AG-359-G18_a22.tar.gz, ).

Each of these tarfiles is a zipped "experiment" folder. (There are 58 experiments in all.)

II. What is an experiment?

It is a set of 27 mutated genomes that are all based on the same source genome (e.g. named 'AG-359-G18_a5') but with different levels of algorithmically generated "mutations"
    So experiment AG-359-G18_a5.tar.gz:
       Is a zipped folder containing 27 fastas named like:
            AG-359-G18_a5_aad0000.fasta
            AG-359-G18_a5_aad0005.fasta
            AG-359-G18_a5_aad0010.fasta
            etc.
        Where the 'aad' number reflects its average amino acid distance from the source genome (again, caused by simulated amino acid changes, as made by Shayestah)

    Note: The pipeline expects genomes with this specific namescheme: [experiment]_aad[number].fasta

III. So what does this pipeline do?

1. For a given experiment, tabulates the pairwise ANI (average nucleotide identity); i.e. between any two given genomes in that experiment.
2. Finds ORFs using Prokka
3. Uses BLAST between ORFs to see whether homology is still detectable between mutated genomes, as the mutation levels increase.
4. Integrates this homology info with ANI, resulting in a table that will let you address the above question.
5. So in total, get 58 tables, for the 58 experiments.


IV. What command did I use to run this nextflow job?
nextflow run 3b_siavash_pt2.nf -profile charlie

BTW: ask me (Greg Gavelis, ggavelis@gbigelow.org) for more information if you want to run nextflow yourself.
*/

params.dev = false // Lets user testrun nextflow command (by adding flag '--dev') which will have this pipeline run on JUST ONE SAG (again, as a test)
params.num_inputs = 2

params.outdir = "results"
params.prokka = "/mnt/scgc_nfs/ref/uniprot_swissprot_prokka.fasta"


workflow {
    CH_num_pair_AND_tar = Channel.fromPath("./input/*.tar.gz")
            .flatten()
            .map { file -> tuple(file.simpleName, file) }
            //.view()

    CH_num_pair_AND_tar = CH_num_pair_AND_tar.take ( params.dev ? params.num_inputs : -1)
    
    // Unzip each experiment (an experiment folder is based on a source genome with 27 permutations (i.e. 27 fastas.) Those permutations are genomes that have been algorithmically "mutated" to have various levels of amino acid distance 'aad' from the source genome.)
    // E.g. for experiment named like "AG-359-G18_a5", there will be fastas like "AG-359-G18_a5_aad00000", "AG-359-G18_a5_aad00005", etc., where 'aad' represents the average amino acid distance from the source genome.
    UNZIP(CH_num_pair_AND_tar)

    // Within an experiment, calculate the pairwise ANI between each of the 27 mutated genomes 
    PYANI(UNZIP.out.groupdir)
    TABULATE_ANI(PYANI.out)

    // From unzipped directory, emit each of the 27 fastas individually, so that prokka can run separately on each.
    // Why? ORFS must be named after their specific genome permutation e.g. ("AG-359-G18_a5_aad00000") so we can tell where they came from.
    CH_fastas = UNZIP.out.fastas.map{ it -> it[1]}.flatten()
    CH_ID_AAD_fasta = CH_fastas.map{ file -> tuple(file.getParent().toString().split('/')[12], file.simpleName, file) } 
    CH_ID_AAD_fasta = CH_ID_AAD_fasta.take ( params.dev ? params.num_inputs : -1)
    // Use prokka to find the ORFs in each genome -> .ffn files.
    PROKKA_v1_14_6(CH_ID_AAD_fasta)

    // Group the ORFs files back together (all 27 back into the same experiment)
    CH_ID_ORF = PROKKA_v1_14_6.out.map{ it -> tuple(it.simpleName.toString().replaceAll(/_aad(.+)/,''), it) }//derives the experiment ID by dropping filename characters starting at "_aad" (e.g. "AG-359-G18_a5_aad00000.fasta" -> "AG-359-G18_a5")
    CH_grouped_ORFS = CH_ID_ORF
                        .groupTuple() // groups all fastas together that share an experiment ID
                        //.view()

    // In each experiment, count up all 27 genomes' worth of ORFS (To ask: Does number of detected ORFs vary across the 27 genome permutations?)
    COUNT_ORFS(CH_grouped_ORFS)

    // In each experiment, BLAST those ORFs against eachother (To ask: Is homology still detectable?)
    BLAST_ORFS_TO_SELF(CH_grouped_ORFS)

    CH_FINAL_TABLES = BLAST_ORFS_TO_SELF.out.join(TABULATE_ANI.out).join(COUNT_ORFS.out)
    SUMMARIZE_EXPERIMENT(CH_FINAL_TABLES)
    }

process BLAST_ORFS_TO_SELF {
    tag "All v. all BLAST for ORFs from experiment ${ID} (across all 27 mutated genomes)"
    cpus = 24
    memory = { 10.GB * task.attempt }
    errorStrategy = 'finish'
    publishDir "results/${ID}/", mode: "copy"
    beforeScript 'module load anaconda/2.1.0'
    input: tuple val(ID), path(FILES)
    output:
        tuple val("${ID}"), path("${ID}_blast.tsv")
    
    shell:
    """
    cat *ffn > "!{ID}_orfs.fasta"

    makeblastdb -in "!{ID}_orfs.fasta" -dbtype nucl

    echo "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > "!{ID}_blast.tsv"

    blastn -query "!{ID}_orfs.fasta" \
    -db "!{ID}_orfs.fasta" \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -num_threads !{task.cpus} -max_target_seqs 10 \
    -out headerless_blast.tsv
    
    cat headerless_blast.tsv >> "!{ID}_blast.tsv"
    """ }

process COUNT_ORFS {
    tag "for experiment ${ID} make a table with the following columns: sag, gene_count, gene_count_500-1500bp (summarizing all 27 mutated genomes)"
    errorStrategy = 'finish'
    beforeScript 'module load anaconda3/2019.07'
    publishDir "results/${ID}/", mode: "copy"
    input: tuple val(ID), path(FILES)
    output:
        tuple val(ID), path("${ID}_orf_counts.csv")
    """
    #!/usr/bin/env python

    PATH_out = "${ID}_orf_counts.csv"
    
    from glob import glob
    import os
    import itertools
    
    length_threshold = [500, 1500]
    LIST_ffn = glob("*.ffn")
    
    def readfa(fh):
        for header, group in itertools.groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = next(group)
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq
    
    with open(PATH_out, 'w') as oh:
        print('sag','gene_count','gene_count_{}-{}bp'.format(length_threshold[0], length_threshold[1]), sep = ',', file = oh)
        for ffn in LIST_ffn:
            sag = os.path.basename(ffn).split(".")[0]
            gene_count = 0
            thresh_count = 0
            for name, seq in readfa(open(ffn)):
                gene_count += 1
                gene_length = len(seq)
                if length_threshold[0] <= gene_length and length_threshold[1] >= gene_length:
                    thresh_count += 1
            #print(sag, gene_count, thresh_count)
            print(sag, gene_count, thresh_count, sep = ',', file = oh)
    """ }

process PROKKA_v1_14_6 {
    tag "Output ORFs (as fnn file) of ${ADD} (which is a mutated genome from experiment ${ID})"
    cpus = 6
    memory = { 10.GB * task.attempt }
    errorStrategy = 'finish'
    publishDir "results/${ID}/orfs_from_prokka/", pattern: "${AAD}.ffn", mode: "copy"
    container='docker://quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_1'
    input: tuple val(ID), val(AAD), path(fasta)
    output: path("${AAD}.ffn")
    """
    prokka --outdir prokka_${AAD} --prefix ${AAD} --locustag ${AAD} --quiet --compliant --force --cpus ${task.cpus} ${fasta}
    mv prokka_${AAD}/*ffn ./${AAD}.ffn
    """ }

process PYANI {
    tag "calculate ANI between each genome permutation (all 27 mutated genomes from experiment ${ID})"
    cpus = 24
    memory = { 10.GB * task.attempt }
    errorStrategy = 'finish'
    beforeScript 'module load anaconda3/2019.07'
    publishDir "results/${ID}/", mode: "copy"
    input: tuple val(ID), path(unzip)
    output: tuple val(ID), path("${ID}_pyani")
    """
    average_nucleotide_identity.py -o "${ID}_pyani" \
    -i ${ID} -m ANIb --workers ${task.cpus}
    """ }

process SUMMARIZE_EXPERIMENT {
    tag "Collate results into single table for experiment ${ID}"
    errorStrategy = 'finish'
    publishDir "results/", mode: "copy"
    beforeScript 'module load anaconda3/2019.07'
    input: tuple val(ID), path(BLAST_tsv), path(ANI_csv), path(ORF_csv)
    output: tuple val(ID), path("${ID}_results.csv")
    """
    #!/usr/bin/env python

    PATH_blast = "${BLAST_tsv}"
    PATH_ani_table = "${ANI_csv}"
    PATH_orf_table = "${ORF_csv}"
    
    PATH_out = "${ID}_results.csv"
    
    TAB = "\t"
    
    import pandas as pd
    import numpy as np
    
    LIST_pid = [99.9, 99.8, 99.7, 99.6, 99.5, 99, 98.5, 98, 97.5, 97]
    len_threshold = [500, 1500]
    min_coverage = 0.03
    
    ### I. Load and filter BLAST results _______________________
    print("I. Loading BLAST results.")
    DF_blast = pd.read_csv(PATH_blast, sep=TAB)
    DF_blast.insert(0, "qsag", DF_blast["qseqid"].str.rsplit("_",n=1).str.get(0))
    DF_blast.insert(3, "ssag", DF_blast["sseqid"].str.rsplit("_",n=1).str.get(0))
    
    print(TAB+str(len(DF_blast["qsag"].unique())) + " SAGs were BLASTed")
    
    LIST_sags = DF_blast["qsag"].unique().tolist()
    
    DF = DF_blast
    print(TAB+str(len(DF_blast)) + " initial BLAST hits")
    
    # Exclude BLAST self-hits
    DF = DF[DF['qsag'] != DF['ssag']]
    print(TAB+"Dropping "+str(len(DF_blast)-len(DF))+" self-hits")
    
    print(TAB+"Dropping "+str(len(DF) - len(DF[(DF['length'] > (0.8 * DF['qlen'])) | (DF['length'] > 0.8 * DF['slen'])]))+ " hits where less than 80% of the query or 80% of the subject seq are aligned")
    DF = DF[(DF['length'] > (0.8 * DF['qlen'])) | (DF['length'] > 0.8 * DF['slen'])]
    
    ### II. DF_pairwise, a table of BLAST-inferred orthologs shared by pairs of SAGs ___________________________________
    
    DF_pairwise = pd.DataFrame(columns = ['qsag','ssag'])
        
    ni_alpha_mean = DF.groupby(['qsag','ssag'], as_index = False)['pident'].mean().rename(columns = {'pident':'mean_pident'}).merge(
                    DF.groupby(['qsag','ssag'])['pident'].std().reset_index().rename(columns = {'pident':'stdev_pident'})).merge(
                    DF.groupby(['qsag','ssag'], as_index = False)['pident'].count().rename(columns = {'pident':'total_hits'})).merge(
                    DF.groupby(['qsag','ssag'], as_index = False)['pident'].median().rename(columns = {'pident':'median_pident'}))
    
    print("II. Integrating ORF counts per SAG")
    DF_orf_count = pd.read_csv(PATH_orf_table)
    DF_pairwise = DF_pairwise.merge(ni_alpha_mean, how = 'outer').merge(
                DF_orf_count.rename(columns = {i:'s{}'.format(i) for i in DF_orf_count.columns}), how = 'left').merge(
                DF_orf_count.rename(columns = {i:'q{}'.format(i) for i in DF_orf_count.columns}), how = 'left')
    
    print(TAB+"Converting to 'DF_pairwise', which quantifies the orthologs detected across "+str(len(DF_pairwise))+" pairs of SAGs")
    
    ### How many orthologs persist across increasing percent identities? 
    
    for pid in LIST_pid:
            
        len_bound = DF[(DF['pident'] >= pid) & 
                ((DF['qlen'] >= len_threshold[0]) & (DF['qlen'] <= len_threshold[1])) & 
                ((DF['slen'] >= len_threshold[0]) & (DF['slen'] <= len_threshold[1]))]\
                .groupby(['qsag','ssag'], as_index=False)['sseqid']\
                .count().rename(columns = {'sseqid':'{}_pid_{}-{}bp_orthologs'.format(pid, len_threshold[0],len_threshold[1])})
        all_sizes = DF[DF['pident'] >= pid].groupby(['qsag','ssag'], as_index=False)['sseqid']\
                .count().rename(columns = {'sseqid':'{}_pid_orthologs'.format(pid, len_threshold[0],len_threshold[1])})
        
        DF_pairwise = DF_pairwise.merge(
                len_bound, how = 'outer').merge(
                all_sizes, how = 'outer')
    
    print(TAB+"Calculating orthologs across "+str(len(LIST_pid))+" percids ranging from " + str(LIST_pid[0])+ " to " +str(LIST_pid[-1]))
    
    # add a pairs column to have a unique identifier per genome pair:
    pairs = []
    for i, l in DF_pairwise.iterrows():
        lst = [l['qsag'], l['ssag']]
        lst.sort()
        pairs.append("_".join(lst))
    DF_pairwise['pair'] = pairs
    
    ortho_names = [i for i in DF_pairwise.columns if 'orthologs' in i]
    
    # Drop reciprocal ANI hits as they are redundant
    DF_pairwise = DF_pairwise.drop_duplicates(subset = 'pair')
    
    ### III. Combine ortholog info (DF_pairwise) with genome ANI (DF_ani) -> DF_final
    
    print("III. Loading ANI data.")
    DF_ani = pd.read_csv(PATH_ani_table)
    #DF_ani = pd.read_csv(PATH_ani_table, index_col=0)
    DF_ani = DF_ani.rename(columns = {'comp':'pair'})
    DF_ani = DF_ani[['Genome A','Genome B','ani', 'ani_aln_cov_ab','ani_aln_cov_ba', 'pair']].dropna() ### Keep only essential columns
    
    print(TAB+"Incorporating ANI data. Deriving GND.")
    DF_final = DF_pairwise.merge(DF_ani, on = 'pair', how = 'outer')
    DF_final['GND'] = [round(100 - i * 100, 2) for i in DF_final['ani']]
    
    # include all pairs for which ANI was measured in summary calculations:
    DF_final[ortho_names] = DF_final[ortho_names].fillna(0)
    
    tups = [(round(i-0.01, 2), round(i, 2),) for i in np.arange(1, 0.6, -0.01)]
    bins = pd.IntervalIndex.from_tuples(tups)
    DF_final['ani_bin'] = pd.cut(DF_final['ani'], bins)
    print(TAB+"Binning ANI across "+str(len(bins))+" increments.")
    
    # write to csv:
    DF_final.to_csv(PATH_out, index=False)
    print("Done.")
    """ }

process TABULATE_ANI {
    tag "Distill 4 pyANI tables into 1 for experiment ${ID}"
    errorStrategy = 'finish'
    publishDir "results/${ID}/", mode: "copy"
    beforeScript 'module load anaconda3/2019.07'
    input: tuple val(ID), path("DIR_ani")
    output: tuple val(ID), path("${ID}_ani.csv")
    """
    #!/usr/bin/env python

    ## What it does: Outputs a table with...
    ## * 1 row per genome-pair
    ## * 5 important columns: Genome A, Genome B, ani, ani_aln_cov_ab, ani_aln_cov_ba, comp

    PATH_out = "${ID}_ani.csv"

    PATH_pyani_identity = "${DIR_ani}/ANIb_percentage_identity.tab"
    PATH_pyani_lengths = "${DIR_ani}/ANIb_alignment_lengths.tab"
    PATH_pyani_coverage = "${DIR_ani}/ANIb_alignment_coverage.tab"
    PATH_pyani_simerrors = "${DIR_ani}/ANIb_similarity_errors.tab"
    
    import pandas as pd
    
    min_coverage = 0.03
    
    TAB = "\t"
    
    print("Loading pyani results.")
    load_matrix_table = lambda metric, table: pd.read_csv(table, sep=TAB).melt(id_vars='Unnamed: 0', var_name='Genome B', value_name=metric).rename(columns={'Unnamed: 0':'Genome A'})
    idf = load_matrix_table('ani', PATH_pyani_identity)
    ldf = load_matrix_table('ani_aln_len', PATH_pyani_lengths)
    cdf = load_matrix_table('ani_aln_cov', PATH_pyani_coverage)
    edf = load_matrix_table('ani_similarity_errs', PATH_pyani_simerrors)
    
    bdf = idf.merge(ldf).merge(cdf).merge(edf)
    bdf = bdf[bdf['Genome A'] != bdf['Genome B']]
    
    print(str(len(bdf))+" genome pairs were compared")
    
    ### Make a table "comp" to show what is being compared.
    
    comps = []
    
    for i, l in bdf.iterrows():
        lst = [l['Genome A'], l['Genome B']]
        lst.sort()
        string = "{}_{}".format(lst[0], lst[1])
        comps.append(string)
    
    bdf['comp'] = comps
    bdf = bdf.sort_values(by=['comp','ani'], ascending=False)
    
    dfa = bdf.drop_duplicates(subset=['comp'], keep='first')
    dfb = bdf.drop_duplicates(subset=['comp'], keep='last')
    
    dfa.rename(columns={'ani_aln_cov':'ani_aln_cov_ab', 'ani_aln_len':'ani_aln_len_ab', 
                           'ani_similarity_errs':'ani_similarity_errs_ab'}, inplace=True)
    dfb.rename(columns={'ani_aln_cov':'ani_aln_cov_ba', 'ani_aln_len':'ani_aln_len_ba', 
                           'ani_similarity_errs':'ani_similarity_errs_ba'}, inplace=True)
    
    dfab = dfa.merge(dfb[['comp', 'ani_aln_cov_ba','ani_aln_len_ba','ani_similarity_errs_ba']], on='comp', how='outer')
    print(str(len(bdf) - len(dfab))+" of the comparisons dropped as being redundant (reciprocal)")
    
    DF = dfab[(dfab['ani_aln_cov_ba'] > min_coverage) & (dfab['ani_aln_cov_ab'] > min_coverage)]
    print(str(len(DF)) + ' ANI estimates had over '+str(min_coverage)+' coverage for both compared genomes')

    DF.to_csv(PATH_out, index=False)
    """ }


process UNZIP {
    tag "Unzip all 27 mutated genomes from experiment ${ID}"
    beforeScript 'module load anaconda3/2019.07'
    input: tuple val(ID), path(targz)
    output:
        tuple val(ID), path("${ID}"), emit: groupdir
        tuple val(ID), path("${ID}/*fasta"), emit: fastas
    """
    echo "unpacking files"
    tar xvzfk ${targz}
    """ }
