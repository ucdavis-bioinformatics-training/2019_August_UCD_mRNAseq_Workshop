## Links to Scripts

### HTStream

slurm script for preprocessing using slurm task array and htstream

[hts_preproc.slurm](./scripts/hts_preproc.slurm)

[hts_preproc_SE.slurm](./scripts/hts_preproc_SE.slurm)

[hts_preproc_TagSeq.slurm](./scripts/hts_preproc_TagSeq.slurm)

shell script for preprocessing using bash loop and htstream.

[hts_preproc.sh](./scripts/hts_preproc.sh)

[hts_preproc_SE.sh](./scripts/hts_preproc_SE.sh)

[hts_preproc_TagSeq.sh](./scripts/hts_preproc_TagSeq.sh)

R script to produce summary table, assumes exact htstream operations and order as described above.

[summary_stats.R](./scripts/summarize_stats.R)

[summary_stats_SE.R](./scripts/summarize_stats_SE.R)

[summary_stats_TagSeq.R](./scripts/summarize_stats_TagSeq.R)

### Star alignment

slurm script for indexing the genome

[star_index.slurm](./scripts/star_index.slurm)

shell script for indexing the genome

[star_index.sh](./scripts/star_index.sh)

slurm script for mapping using slurm task array and star

[star.slurm](./scripts/star.slurm)

[star_SE.slurm](./scripts/star_SE.slurm)

shell script for mapping using bash loop and star.

[star.sh](./scripts/star.sh)

[star_SE.sh](./scripts/star_SE.sh)

shell script to produce summary mapping table

[star_stats.sh](./scripts/star_stats.sh)

### Salmon alignment

slurm script for indexing the genome

[salmon_index.slurm](./scripts/salmon_index.slurm)

shell script for indexing the genome

[salmon_index.sh](./scripts/salmon_index.sh)

slurm script for mapping using slurm task array and star

[salmon.slurm](./scripts/salmon.slurm)

shell script for mapping using bash loop and star.

[salmon.sh](./scripts/salmon.sh)

R script to produce summary mapping table

[salmon_stats.R](./scripts/salmon_stats.R)

### Matt's R scripts from final day

Matt's final Rmd file for DE_analysis.Rmd

[Matt_DE_Analysis.Rmd](./differential_expression/Matt_DE_Analysis.Rmd)

Matt's final Rmd file for enrichment.Rmd

[Matt_enrichment.Rmd](./differential_expression/Matt_enrichment.Rmd)
