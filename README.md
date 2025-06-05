Processing PacBio data starting from FLNC reads.

The skeleton of the SLURM script to run this pipeline

```
#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --account={account}
#SBATCH --time=90:00:00
#SBATCH --mem=9000

module load java/21 apptainer

nextflow run -profile uva -params-file params.yml aakrosh/nf-pacbio-rna -resume
```

The params.yml file should be 

```
apptainer_cache_dir: '{cache_dir}'
clusterOptions: '--account={account}'
queue: 'standard-afton,standard-rivanna'
gtf: '{gtf}'
reference: '{reference}'
reference_faidx: '{reference_fai}''
scratch: '{tmp}'
genome: 'hg38'
platform: 'Revio'
flnc_bams:
  - id: "{sample1}"
    path: "{bam2}"
  - id: "{sample2}"
    path: "{bam2}
```
