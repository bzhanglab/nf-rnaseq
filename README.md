# RNA-seq pipeline

### To run:

* run with docker locally

```console
nextflow run bzhanglab/nf-rnaseq -r main -profile docker \
   --case_id /path/to/my_case_id.txt \
   --start 1 --end 2 --run_version rnaseq-2020-09-22
```

* run with aws batch on AWS cloud

```console
nextflow run bzhanglab/nf-rnaseq -r main -profile awsbatch \
   --start 1 --end 2 \
   --case_id /path/to/my_case_id.txt \
   -bucket-dir s3://zhanglab-nextflow-workdir/workdir/2020-09-18 \
   --outdir s3://zhanglab-shiz/rnaseq-results \
   --run_version rnaseq-2020-09-22
```

In the above examples, if you want the results to be stored locally on your 
launch host, just provide `--run_version xyz`, and the results will be stored
at `results/xyz` under your launch directory. 

To run with awsbatch, you must specify an s3 path as `outdir`, e.g.
`--outdir s3://my-bucket/my-folder`.  In this case, results will be 
stored under `s3://my-bucket/my-folder/run_version`.


## Inputs

Required:

* CPTAC3 catalog file from https://github.com/ding-lab/CPTAC3.catalog/blob/master/CPTAC3.Catalog.dat (S3 path specified in the config file)

User inputs:

* Case ids file called `case_id.txt` contains all the case ids to be processed in a column.
  This file can be in the local directory where nextflow is launched or in S3.
  

## Data analysis steps

1. Build bwa index. The index will be used for all samples.

```console 
   tools/bwa-0.7.17/bwa index genome/GRCh38.d1.vd1.fa
```

2. Map RNAseq to genome using bwa. The mapped sam will be used for circRNA calling.
   
```console
   tools/bwa-0.7.17/bwa mem -t 16 -T 19 genome/GRCh38.d1.vd1.fa fastq/C3L-00006_R1.fastq fastq/C3L-00006_R2.fastq -o sam/bwa-pe-for-CIRI.sam 
```

3. circRNA calling using CIRI.

```console
    mkdir CIRI_results & perl tools/CIRI_v2.0.6/CIRI2.pl -T 16 -I sam/bwa-pe-for-CIRI.sam -O CIRI_results/results_CIRI.txt -F genome/GRCh38.d1.vd1.fa -A genome/UCSC_hg38_refseq_05_24_2019_with_gene_symbol_isoform_fixed.gtf 
```

4. Add gene names to CIRI outputs by the inhouse script.

```console
   perl 1_add_gene_name_to_CIRI_results.pl 
```

5. Add circRNA to gene annotation in the gft format by the inhouse script.

```console
   perl 2_add_linear_circular_isoform_to_gtf.pl 
```

6. Extract both linear and circRNA transcripts using RSEM.

```console
   tools/RSEM-1.3.1/rsem-extract-reference-transcripts genome/GRCh38.d1.vd1.linear.and.circrna 0 genome/UCSC_hg38_refseq_05_24_2019_with_gene_symbol_isoform_fixed_with_circular_RNA.gtf None 0 genome/GRCh38.d1.vd1.fa 
```

7. Transfer linear transcript of circRNA to psedo linear transcript. It also will remove transcript with length less than reads length and any transcritps with “N”. 

```console
   perl 3_circular_linear_to_psedo_linear.pl
```

8. Generate transcript and gene mapping table for RSEM index.

```console
   perl 4_gene_isoform_mapping_with_circular_rna.pl 
```

9. Build RSEM index using transcript from step 7 and mapping from step 8

```console
   tools/RSEM-1.3.1/rsem-prepare-reference -p 16 --transcript-to-gene-map genome/gene_isoform_mapping_with_circular_rna_for_RSEM.txt --bowtie2 --bowtie2-path tools/bowtie2-2.3.3/ genome/GRCh38.d1.vd1.linear.and.circrna.as.pseudo.linear.transcripts.fa RSEM_index/hg38 
```

10. run RSEM to do gene and transcript quantification.

```console
    tools/RSEM-1.3.1/rsem-calculate-expression --bowtie2 --bowtie2-path tools/bowtie2-2.3.3/ --paired-end -p 16 fastq/C3L-00006_R1.fastq fastq/C3L-00006_R2.fastq RSEM_index/hg38 RSEM_results/C3L-00006 
```

11. Summarize RSEM output.

```console
    perl 5-summary-gene-quantification.pl 
```

12. Combine and summarize

```console
    Rscript combine-and-summary.r
```