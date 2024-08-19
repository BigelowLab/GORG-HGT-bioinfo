#PBS -N launch_nextflow_pipeline
#PBS -M ggavelis@bigelow.org
#PBS -m bea -j oe -o status_nextflow.pbs
#PBS -l ncpus=1,mem=50MB,walltime=120:00:00

cd /mnt/scgc_nfs/lab/ggavelis/analysis/siavash_hgt/aug9_second_round/replicate_on_pt1/

echo "launching nextflow pipeline"
module use /mnt/scgc_nfs/opt/modulefiles/common/
module load anaconda/2.1.0
module load singularity

export NXF_EXECUTOR=pbspro
export SINGULARITY_CACHEDIR="/mnt/scgc_nfs/containers"
export SINGULARITY_TMPDIR="/mnt/scgc_nfs/containers/temp"
export NXF_SINGULARITY_CACHEDIR="/mnt/scgc_nfs/containers"
export NXF_CACHE_DIR="/mnt/scgc/scgc_temp/nextflow"
export NXF_CONDA_ENABLED=true

export NXF_VER=21.10.4
NXF_VER=21.10.4 nextflow run 2b_hgt_pipeline_pt1.nf -profile charlie

echo "done"
