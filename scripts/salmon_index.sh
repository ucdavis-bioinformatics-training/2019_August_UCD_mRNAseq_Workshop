#!/bin/bash

## assumes star version 2.7.0e

start=`date +%s`
echo $HOSTNAME

outpath="References"
[[ -d ${outpath} ]] || mkdir ${outpath}

cd ${outpath}
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.pc_transcripts.fa.gz
gunzip gencode.v31.pc_transcripts.fa.gz
PC_FASTA="gencode.v31.pc_transcripts.fa"

module load salmon
call="salmon index -i salmon_gencode.v31.index -k 31 --gencode -p 8 -t ${PC_FASTA}

mkdir star.overlap100.gencodev31
cd star.overlap100.gencodev31

call="STAR
     --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir . \
     --sjdbOverhang 100 \
     --sjdbGTFfile ${GTF} \
     --genomeFastaFiles ${FASTA}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime