# Alignment using Salmon

## Salmon Aligner
[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) is an aligner that uses transcripts for getting raw counts from RNA-Seq data. We will align using salmon to compare count data with our STAR alignments.

### Indexing the transcripts

**1\.** First go to your rnaseq_example directory, load the salmon module, and take a look at the options:

	cd /share/workshop/$USER/rnaseq_example
	module load salmon
	salmon -h

-----

**2\.** Look at the help docs for the salmon subcommands as well:

	salmon index -h
	salmon quant --help-reads

-----

**3\.** In order to align the data we first need to download and index the transcripts. Pull down a slurm script to do the downloading and indexing and take a look at it:

	wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/scripts/salmon_index.slurm
	less salmon_index.slurm

Press 'q' to exit.

> #!/bin/bash
>
> #SBATCH --job-name=salmon_index # Job name<br>
> #SBATCH --nodes=1<br>
> #SBATCH --ntasks=8<br>
> #SBATCH --time=60<br>
> #SBATCH --mem=15000 # Memory pool for all cores (see also --mem-per-cpu)<br>
> #SBATCH --partition=production<br>
> #SBATCH --reservation=workshop<br>
> #SBATCH --account=workshop<br>
> #SBATCH --output=slurmout/salmon-index_%A.out # File to which STDOUT will be written<br>
> #SBATCH --error=slurmout/salmon-index_%A.err # File to which STDERR will be written<br>
>
> start=\`date +%s\`<br>
> echo $HOSTNAME<br>
>
> outpath="References"<br>
>
> cd ${outpath}<br>
> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz<br>
> gunzip gencode.v29.transcripts.fa.gz<br>
>
> module load salmon<br>
> call="salmon index -i salmon_index -p 8 -t gencode.v29.transcripts.fa --gencode"<br>
>
> echo $call<br>
> eval $call<br>
>
> end=\`date +%s\`<br>
> runtime=$((end-start))<br>
> echo $runtime<br>

1. The script changes into the References directory.
1. It uses wget to download the transcript fasta file from GENCODE.
1. Uncompresses it using gunzip.
1. Run Salmon indexing, using the "gencode" flag to parse the GENCODE file properly, and outputting to a new directory called "salmon_index".

-----

**4\.** Submit the script:

	sbatch salmon_index.slurm

This step will take about 11 minutes to complete after the job actually starts running. Take a look at the [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html) while you wait. All of the output will be in the "salmon_index" directory.

**IF** for some reason it didn't finish, is corrupted, or you missed the session, you can copy over a completed copy.

	cp -r /share/biocore/workshops/2019_August_RNAseq/References/salmon_index /share/workshop/$USER/rnaseq_example/References/

-----

### Alignment

**5\.** Once your indexing is complete, then you can align all the samples. Go back to your rnaseq_example directory, download the slurm script for the alignment, and take a look at it:

	cd /share/workshop/$USER/rnaseq_example
	wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/scripts/salmon.slurm
	less salmon.slurm

Press 'q' to exit.

> #!/bin/bash<br>
>
> #SBATCH --array=1-16<br>
> #SBATCH --job-name=salmon # Job name<br>
> #SBATCH --nodes=1<br>
> #SBATCH --ntasks=8<br>
> #SBATCH --time=1440<br>
> #SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)<br>
> #SBATCH --partition=production<br>
> #SBATCH --reservation=workshop<br>
> #SBATCH --account=workshop<br>
> #SBATCH --output=slurmout/salmon_%A_%a.out # File to which STDOUT will be written<br>
> #SBATCH --error=slurmout/salmon_%A_%a.err # File to which STDERR will be written<br>
>
>
> start=\`date +%s\`<br>
> hostname<br>
>
> outdir="02-Salmon"<br>
> sampfile="samples.txt"<br>
> REF="References/salmon_index"<br>
> GTF="References/gencode.v29.primary_assembly.annotation.gtf"<br>
>
> SAMPLE=\`head -n ${SLURM_ARRAY_TASK_ID} $sampfile | tail -1\`<br>
> R1="01-HTS_Preproc/$SAMPLE/${SAMPLE}\_R1.fastq.gz"<br>
> R2="01-HTS_Preproc/$SAMPLE/${SAMPLE}\_R2.fastq.gz"<br>
>
> echo $SAMPLE<br>
>
> if [ ! -e $outdir ]; then<br>
>   mkdir $outdir<br>
> fi<br>
>
> module load salmon<br>
> call="salmon quant -p 8 -i $REF -l A \<br>
> --validateMappings -g $GTF \<br>
> -1 $R1 -2 $R2 \<br>
> -o $outdir/$SAMPLE"<br>
>
> echo $call<br>
> eval $call<br>
>
> end=\`date +%s\`<br>
> runtime=$((end-start))<br>
> echo Runtime: $runtime seconds<br>

1. The script specifies the output directory (02-Salmon), the samples file (samples.txt), the reference that we just indexed, and the annotation that we downloaded when we ran STAR.
1. It then defines the filenames for the forward and reverse reads (R1 and R2).
1. Creates the output directory.
1. Loads the salmon module and then runs salmon in "quant" mode to quantify (i.e. count) the reads aligning to transcripts. Salmon uses the annotation GTF file to roll up the transcript counts into gene-level counts.

-----

**6\.** Submit the script:

	sbatch salmon.slurm

This will run a job for every sample... each job should only take around 3 minutes to run.

**IF** for some reason it didn't finish, is corrupted, or you missed the session, you can copy over a completed copy.

	cp -r /share/biocore/workshops/2019_August_RNAseq/02-Salmon /share/workshop/$USER/rnaseq_example/

-----

**7\.** Once the jobs finish, take a look at one of the output files:

	cd 02-Salmon/SampleAC1
	less quant.genes.sf

These are the gene-level counts rolled up from the transcript counts. Take a look at the [Salmon output file format documentation](https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats) to understand each of the columns. Press 'q' to exit. We will be comparing the counts from Salmon to the counts from STAR, to see their differences, if any.