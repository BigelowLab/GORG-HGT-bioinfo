#!/user/bin/nextflow
nextflow.enable.dsl=2

/*
I. Inputs:

Expects an input folder called ./gzipped_experiments/ that contains several tarfiles (e.g. AG-359-G18_a5.tar.gz, AG-359-G18_a22.tar.gz, ).

Each of these tarfiles is a zipped "experiment" folder. (There are 58 experiments in all.)

II. What is an experiment?

It is a set of 31 mutated genomes that are all based on the same source genome (e.g. named 'AG-359-G18_a5') but with different levels of algorithmically generated "mutations"
    So experiment AG-359-G18_a5.tar.gz:
       Is a zipped folder containing 31 fastas named like:
            AG-359-G18_a5_gnd0000.fasta
            AG-359-G18_a5_gnd0005.fasta
            AG-359-G18_a5_gnd0010.fasta
            etc.
        Where the 'gnd' number reflects its general nucleotide distance from the source genome (again, caused by simulated random mutations, as made by Siavash et al.)

    Note: The pipeline expects genomes with this specific namescheme: [experiment]_gnd[number].fasta

III. So what does this pipeline do?

1. For a given experiment, tabulates the pairwise ANI (average nucleotide identity); i.e. between any two given genomes in that experiment.
2. Finds ORFs using Prokka
3. Uses BLAST between ORFs to see whether homology is still detectable between mutated genomes, as the mutation levels increase.
4. Integrates this homology info with ANI, resulting in a table that will let you address the above question.
5. So in total, get 58 tables, for the 58 experiments.


IV. What command did I use to run this nextflow job?
nextflow run orfs_vs_ani.nf

BTW: ask me (Greg Gavelis, ggavelis@gbigelow.org) for more information if you want to run nextflow yourself.
*/

params.dev = false // Lets user testrun nextflow command (by adding flag '--dev') which will have this pipeline run on JUST ONE SAG (again, as a test)
params.num_inputs = 1

params.outdir = "results"

params.alpha10 = "/alpha/alpha10-c1-gl904-orf1.csv"
params.alpha6 = "/alpha/alpha6.74-c1-gl904-orf1.csv"
params.alpha4 = "/alpha/alpha4.59-c1-gl904-orf1.csv"

workflow {
    CH_num_pair_AND_tar = Channel.fromPath("./mid/2_gzipped_inputs/*.tar.gz")
            .flatten()                                                                                      // Emit each demultiplexed fastq as its own object. E.g. blahblah_R1.fastq.gz
            .map { file -> tuple(file.simpleName, file) }

    // Unzip each experiment (an experiment folder is based on a source genome with 31 permutations (i.e. 31 fastas.) Those permutations are genomes that have been algorithmically "mutated" to have various levels of general nucleotide distance 'gnd' from the source genome.)
    // E.g. for experiment named like "AG-359-G18_a5", there will be fastas like "AG-359-G18_a5_gnd00000", "AG-359-G18_a5_gnd00005", etc., where 'gnd' represents the average general nucleotide distance from the source genome.
    UNZIP(CH_num_pair_AND_tar)

    // Within an experiment, calculate the pairwise ANI between each of the 31 mutated genomes 
    PYANI(UNZIP.out.groupdir)
    TABULATE_ANI(PYANI.out)

    // From unzipped directory, emit each of the 31 fastas individually, so that prokka can run separately on each.
    // Why? ORFS must be named after their specific genome permutation e.g. ("AG-359-G18_a5_gnd00000") so we can tell where they came from.
    CH_fastas = UNZIP.out.fastas.map{ it -> it[1]}.flatten()
    CH_ID_AAD_fasta = CH_fastas.map{ file -> tuple(file.getParent().toString().split('/')[12], file.simpleName, file) } 
    PROKKA_v1_14_6(CH_ID_AAD_fasta) // Use prokka to find the ORFs in each genome -> .ffn files.

    // Group the ORFs files back together (all 31 back into the same experiment)
    CH_ID_ORF = PROKKA_v1_14_6.out.map{ it -> tuple(it.simpleName.toString().replaceAll(/_gnd(.+)/,''), it) }//derives the experiment ID by dropping filename characters starting at "_gnd" (e.g. "AG-359-G18_a5_gnd00000.fasta" -> "AG-359-G18_a5")
    CH_grouped_ORFS = CH_ID_ORF.groupTuple() // groups all fastas together that share an experiment ID

    // In each experiment, count up all 31 genomes' worth of ORFS (To ask: Does number of detected ORFs vary across the 31 genome permutations?)
    COUNT_ORFS(CH_grouped_ORFS)

    // In each experiment, BLAST those ORFs against eachother (To ask: Is homology still detectable?)
    BLAST_ORFS_TO_SELF(CH_grouped_ORFS)

    CH_FINAL_TABLES = BLAST_ORFS_TO_SELF.out.join(TABULATE_ANI.out).join(COUNT_ORFS.out)
    SUMMARIZE_EXPERIMENT(CH_FINAL_TABLES)
    CH_all_tables = SUMMARIZE_EXPERIMENT.out.csv.collect()
    CH_each_alpha_AND_table = CH_all_tables.flatten().map{ it -> tuple(it.simpleName.replaceAll('_results','').replaceAll(/(.+)_a/,''), it)}
    CH_grouped_by_alpha = CH_each_alpha_AND_table.groupTuple()

    CH_alpha_1SAGtable_vs_allSAGtablesSameAlpha = CH_each_alpha_AND_table.combine(CH_grouped_by_alpha, by:0)
    
    CH_SAGid_alpha_1SAGtable_vs_allSAGTablesSameAlpha = CH_alpha_1SAGtable_vs_allSAGtablesSameAlpha.map{ it -> tuple(it[1].simpleName.toString().replaceAll(/_(.+)/,''), it[0], it[1], it[2]) } // Adds 'SAGid' (e.g. "AM-359-G18") to tuple, for sake of naming output files

    // Allows --dev mode to run GGPLOT on just 1 sample emitted by this channel
    CH_SAGid_alpha_1SAGtable_vs_allSAGTablesSameAlpha = CH_SAGid_alpha_1SAGtable_vs_allSAGTablesSameAlpha.take ( params.dev ? params.num_inputs : -1)
    
    GGPLOT(CH_SAGid_alpha_1SAGtable_vs_allSAGTablesSameAlpha)
    }



process GGPLOT {
    errorStrategy = 'finish'
    beforeScript 'module load anaconda'
    publishDir "results/plots/${alpha}/divergence/", pattern: "${ID}_${alpha}_divergence.pdf", mode: "copy"
    publishDir "results/plots/${alpha}/shared-genes-percent/", pattern: "${ID}_${alpha}_shared-genes-percent.pdf", mode: "copy"

    conda '/mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/r_ggplot2'
    // env r_ggplot2 created like so: module load anaconda3; conda create --prefix /mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/r_ggplot2 conda-forge::r-ggpubr conda-forge::r-ggplot2 conda-forge::r-dplyr

    input: tuple val(ID), val(alpha), path(main_csv, stageAs: "?/*"), path(other_csvs)
    output: tuple val(ID), path("${ID}_${alpha}_divergence.pdf"), path("${ID}_${alpha}_shared-genes-percent.pdf")

    script:
    """
    #!/usr/bin/env Rscript

    require(ggplot2);require(scales); require(reshape2); library(ggpubr);library(dplyr)

    #a <- "5"
    #bd <- read.csv(paste0("./result_tables_2000/AG-359-G18_a", a, "_results.csv"))
    bd <- read.csv("${main_csv}")

    all_files  = list.files(path="./", pattern="*.csv", full.names=TRUE, recursive=FALSE)
    files <- all_files[ !grepl("${ID}", all_files) ] # Exclude the file that is redundant with the query SAG
    
    for (g in files) {
    data <- read.csv(g)
    bd <- rbind(bd, data)
    }

    head(bd)

    nrow(bd)
    
    bd_filtered <- bd[, c("qsag", "ssag", "mean_pident", "GND")]
    bd_filtered\$est_GND <- (100 - bd_filtered\$mean_pident)/100

    extract_base_genome <- function(string) {
      strsplit(string, "_")[[1]][1]
    }

    nrow(bd)
    head(bd)
    
    bd\$gndbin=cut(bd\$mean_pident/100,breaks=(670:2000)/2000, labels =round((1-(671:2000)/2000+0.0005)*100,digits=5))
    dcast(data=bd[,c("X99_pid_500.1500bp_orthologs","gndbin")],formula=gndbin~.)
    nrow(bd)
    head(bd)
    nrow(bd[bd\$X99_pid_500.1500bp_orthologs == 0,])
    
    (bd %>% mutate(SAG = sub("_a5_gnd.*\$","",ssag)) %>%
        filter(1-mean_pident/100 < 0.0042) %>% 
        filter(1-mean_pident/100 > 0.0037) %>% 
        filter(SAG == "AG-892-F15") 
    )
    
    #Completeness is 1 for simulated dataset
    bdrm <- bd
    bdrm\$CompA <- 1
    bdrm\$CompB <- 1
    ncol(bdrm)
    head(bdrm)
    
    bdrmi=bdrm[bdrm\$total_hits> 800-50*(100-bdrm\$mean_pident),]
    nrow(bdrmi)
    head(bdrmi)

    ### This is what we used
    alphas=quantile(with(bdrmi[bdrmi\$mean_pident!= 100 & bdrmi\$stdev_pident!=0 &  bdrmi\$mean_pident <95 & bdrmi\$mean_pident > 80  &
                              !is.na(bdrmi\$total_hits)   & bdrmi\$total_hits > 4,
                            ,], 
                         1/( (stdev_pident/100) /(1-mean_pident/100) )^2 ),c(0.1,0.5,0.9))
    
    alphas
    
    
    model=rbind(
      data.frame(read.csv('${params.alpha4}'),alpha=4.59,c=1,gl=904,orf=1),
      data.frame(read.csv('${params.alpha6}'),alpha=6.74,c=1,gl=904,orf=1),
      data.frame(read.csv('${params.alpha10}'),alpha=10,c=1,gl=904,orf=1)
    )
    
    head(model)
    model\$GND = round(model\$GND,5)
    tail(model)
    
    head(bdrmi)
    
    # bdrmim = melt(bdrmi,measure.vars = 10:23)
    bdrmim = melt(bdrmi,measure.vars = 11:30)
    nrow(bdrmim)
    bdrmim = bdrmim[grepl("500.1500bp" ,bdrmim\$variable),]
    bdrmim\$adjusted = with(bdrmim,value*(1/sgene_count_500.1500bp/CompB+1/qgene_count_500.1500bp/CompA)/2 )
    bdrmim\$similarity=sub("_.*","",bdrmim\$variable)
    head(bdrmim)
    levels(bdrmim\$variable)
    bdrmim\$SAG = sub("_a5_gnd.*\$","",bdrmim\$ssag)

    # View(bdrmim[bdrmim\$similarity=="X99" & bdrmim\$adjusted > 0.8 & bdrmim\$adjusted < 0.85 ,])
    
    #View(bdrmim[bdrmim\$similarity=="X99",] %>% 
    #  filter(1-mean_pident/100 < 0.006) %>% 
    #  filter(1-mean_pident/100 > 0.002) %>% 
    #  filter(SAG == "AG-892-F15") %>%
    #  filter(grepl(".*_gnd00000",qsag)|grepl(".*_gnd00000",ssag))
    #  )
    
    ggplot(aes(y=  adjusted ,
               color="Data", x=(1-mean_pident/100)),
           data=bdrmim[bdrmim\$similarity=="X99" ,])+
      #geom_hline(yintercept = 0.96)+
      #geom_vline(xintercept = c(0.0026,0.003))+
      geom_point(aes(color=ifelse(grepl(".*AG-891-K05",ssag),"AG-891-K05"," others")),alpha=0.75,size=0.95)+ #geom_smooth(se=F,aes(color="Data (fit)"),size=1)+
      geom_line(aes(y=`.`,color="Data (mean)",x=as.numeric(as.character(gndbin))/100),
                data=dcast(gndbin+variable~.,data=bdrmim[bdrmim\$similarity=="X99",c("gndbin","variable", "adjusted")],fun.aggregate = mean),
                size=1)+
      geom_ribbon(aes(ymin=`4.59`,y=`6.74`,ymax=`10`,x=GND/100,color="Model (80% CI alpha)"),
                  data=dcast(GND~alpha,data=model[,c("GND", "X99", "alpha")],value.var = "X99"),
                  size=0.2,alpha=0.3,fill="#EE60BB",show.legend = F)+
      geom_line(aes(y=X99/orf/c,x=GND/100, color="Model (median alpha)"),
                data=model[model\$alpha==6.74,],
                size=.81,alpha=0.68)+
      scale_x_continuous(name="GND",labels=percent)+
      scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
      scale_color_manual(name="",values = c("gray50","#1055EE","#FF60AA","#BB3333","#770010"))+
      scale_y_continuous(name="Shared genes (adjusted for incompleteness)",labels=percent,breaks=c(0,0.2,0.4,0.6,0.8,1))+
      theme_bw()+
      theme(legend.position = c(.8,.7),panel.grid.minor = element_blank())+
      coord_cartesian(xlim = c(0,0.04))
    ggsave("${ID}_${alpha}_shared-genes-percent.pdf",width=6.5,height = 4.5)


    ds = merge(dcast(variable+GND~alpha,data=melt(model[,c(1,3:12,14)],measure.vars = 2:11)[,c("GND","alpha","variable","value")],value.var = "value"),
               dcast(gndbin+similarity~"real",data=bdrmim[,c("gndbin","similarity", "adjusted")],fun.aggregate = mean),
               by.x=c("variable","GND"),by.y=c("similarity","gndbin"))
    ds=melt(ds,measure.vars = 3:5,variable.name = "alpha",value.name = "model")
    
    ds = merge(dcast(variable+GND~alpha,data=melt(model[,c(1,3:12,14)],measure.vars = 2:11)[,c("GND","alpha","variable","value")],value.var = "value"),
               dcast(gndbin+similarity~"real",data=bdrmim[,c("gndbin","similarity", "adjusted")],fun.aggregate = mean),
               by.x=c("variable","GND"),by.y=c("similarity","gndbin"))
    ds=melt(ds,measure.vars = 3:5,variable.name = "alpha",value.name = "model")
    
    
    ggplot(aes(y=real-model,x=GND/100,color=as.factor(100-as.numeric(sub("X","",variable)))),data=ds[ds\$alpha == 6.74 & ds\$variable %in% c("X99","X98.5","X98"),])+
      geom_line(size=0.8)+
      geom_point(size=0.8)+
      scale_color_manual(name="Gene ND threshold",values=c("red","#009900","blue"))+
      scale_x_continuous(lim=c(0,0.13),labels = function(x) x*100,breaks = (0:13)/100,"Genomic nucleotid difference, %")+
      #facet_wrap(~alpha,labeller = label_both,nrow=2)+
      theme_classic()+
      geom_hline(yintercept = 0,color="grey")+
      theme(legend.position = c(.8,.85))+
      coord_cartesian(ylim = c(-0.15,0.25))
    ggsave("${ID}_${alpha}_divergence.pdf",width=6.5,height = 4.5)
    """
}

process BLAST_ORFS_TO_SELF {
    tag "All v. all BLAST for ORFs from experiment ${ID} (across all 31 mutated genomes)"
    cpus = 24
    memory = { 10.GB * task.attempt }
    errorStrategy = 'finish'
    publishDir "results/raw/${ID}/", mode: "copy"
    beforeScript 'module load anaconda/2.1.0' // load the environment in which BLAST is installed
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
    -num_threads !{task.cpus} -max_target_seqs 2000 \
    -out headerless_blast.tsv
    
    cat headerless_blast.tsv >> "!{ID}_blast.tsv"
    """ }

process COUNT_ORFS {
    tag "for experiment ${ID} make a table with the following columns: sag, gene_count, gene_count_500-1500bp (summarizing all 31 mutated genomes)"
    errorStrategy = 'finish'
    beforeScript 'module load anaconda3/2019.07' // python environment
    publishDir "results/raw/${ID}/", mode: "copy"
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
    publishDir "results/raw/${ID}/orfs_from_prokka/", pattern: "${AAD}.ffn", mode: "copy"
    container='docker://quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_1'
    input: tuple val(ID), val(AAD), path(fasta)
    output: path("${AAD}.ffn")
    """
    prokka --outdir prokka_${AAD} --prefix ${AAD} --locustag ${AAD} --quiet --compliant --force --cpus ${task.cpus} ${fasta}
    mv prokka_${AAD}/*ffn ./${AAD}.ffn
    """ }

process PYANI {
    tag "calculate ANI between each genome permutation (all 31 mutated genomes from experiment ${ID})"
    cpus = 24
    memory = { 10.GB * task.attempt }
    errorStrategy = 'finish'
    beforeScript 'module load anaconda3/2019.07' // environment with pyani v. 0.2.9 installed
    publishDir "results/raw/${ID}/", mode: "copy"
    input: tuple val(ID), path(unzip)
    output: tuple val(ID), path("${ID}_pyani")
    """
    average_nucleotide_identity.py -o "${ID}_pyani" \
    -i ${ID} -m ANIb --workers ${task.cpus}
    """ }

process SUMMARIZE_EXPERIMENT {
    tag "Collate results into single table for experiment ${ID}"
    errorStrategy = 'finish'
    publishDir "results/raw/", mode: "copy"
    beforeScript 'module load anaconda3/2019.07'  // environment with pandas
    input: tuple val(ID), path(BLAST_tsv), path(ANI_csv), path(ORF_csv)
    output:
        tuple val(ID), path("${ID}_results.csv"), emit: tuple
        path("${ID}_results.csv"), emit: csv
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
    publishDir "results/raw/${ID}/", mode: "copy"
    beforeScript 'module load anaconda3/2019.07' // environment with pandas
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
    tag "Unzip all 31 mutated genomes from experiment ${ID}"
    beforeScript 'module load anaconda3/2019.07'
    input: tuple val(ID), path(targz)
    output:
        tuple val(ID), path("${ID}"), emit: groupdir
        tuple val(ID), path("${ID}/*fasta"), emit: fastas
    """
    echo "unpacking files"
    tar xvzfk ${targz}
    """ }
