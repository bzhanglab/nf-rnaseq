#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir_run = "${params.outdir}/${params.run_version}"
params.outdir_ciri = "${params.outdir_run}/CIRI_results"
params.outdir_rsem = "${params.outdir_run}/RSEM_results"
params.outdir_rsem2 = "${params.outdir_run}/RSEM_results_2"


def helpMessage() {
    log.info rnaseqHeader()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile docker 
    or 
    nextflow run main.nf -profile awsbatch -bucket-dir s3://zhanglab-nextflow-workdir/workdir
    Mandatory arguments:
      -profile                      Configuration profile to use.
                                    Available: docker, awsbatch.
    
    Optional argument:
      --allow_dup                  Allow duplicated sample id 
    """.stripIndent()
}

header = """====================================================================================
   ___   _  __ ___    ____ ____ ____     ___   ____ ___   ____ __    ____ _  __ ____
  / _ \\ / |/ // _ |  / __// __// __ \\   / _ \\ /  _// _ \\ / __// /   /  _// |/ // __/
 / , _//    // __ | _\\ \\ / _/ / /_/ /  / ___/_/ / / ___// _/ / /__ _/ / /    // _/  
/_/|_|/_/|_//_/ |_|/___//___/ \\___\\_\\ /_/   /___//_/   /___//____//___//_/|_//___/ 
===================================================================================
"""

def rnaseqHeader() {
    return header.stripIndent()
}

// step 1
process genome_index {
  container "${params.container.bwa}" 
  publishDir  params.genome_basedir, mode: 'move'
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
  val(ref_prefix)
  path('ref.fa')
   
  output:
  path("${ref_prefix}.*")
  
  """
  mv ref.fa ${ref_prefix}
  bwa index ${ref_prefix}
  """
}


process generate_id_files_gdc {
   container "${params.container.r_tidyverse}"
   label 'r5_2xlarge'

   input:
     path('catalog.txt')
     path('case_id.txt')
     val start
     val end

   output:
     path 'case_*.tsv'

   script:
   """
   #!/usr/bin/env Rscript
   library(tidyverse)
   
   # if dup_id is 0,  there is no duplication
   get_file_info <- function(catalog_rnaseq, case, st, r1_name, dup_id,
                             matched_case, case_aliquot){
      if (dup_id == 0) {
        new_name <- paste0(case, "_", st)
        cur_aliquot <- case_aliquot[1, "aliquot"]
      } else {
        new_name <- paste0(case, "_", st, "_", dup_id)
        cur_aliquot <- case_aliquot[dup_id, "aliquot"] %>% pull(1)
      }
      r1_file_name <- catalog_rnaseq %>%
                  filter(.data[[names(.)[[1]]]] == r1_name) %>%
                  filter(.data[[names(.)[[6]]]] == cur_aliquot) %>%
                  select("filename") %>%
                  pull(1)
      r1_uuid <- catalog_rnaseq %>%
                  filter(.data[[names(.)[[1]]]] == r1_name) %>%
                  filter(.data[[names(.)[[6]]]] == cur_aliquot) %>%
                  select("UUID") %>%
                  pull(1)
      r1_md5 <- catalog_rnaseq %>%
                  filter(.data[[names(.)[[1]]]] == r1_name) %>%
                  filter(.data[[names(.)[[6]]]] == cur_aliquot) %>%
                  select("MD5") %>%
                  pull(1)
      r2_name <- matched_case %>%
                  filter(.data[[names(.)[[3]]]] == "R2") %>%
                  filter(.data[[names(.)[[4]]]] == st) %>%
                  pull(1) %>% 
                  unique()
      r2_file_name <- catalog_rnaseq %>%
                    filter(.data[[names(.)[[1]]]] == r2_name) %>%
                    filter(.data[[names(.)[[6]]]] == cur_aliquot) %>%
                    select("filename") %>%
                    pull(1)
      r2_uuid <- catalog_rnaseq %>%
                  filter(.data[[names(.)[[1]]]] == r2_name) %>%
                  filter(.data[[names(.)[[6]]]] == cur_aliquot) %>%
                  select("UUID") %>%
                  pull(1)
      r2_md5 <- catalog_rnaseq %>%
                  filter(.data[[names(.)[[1]]]] == r2_name) %>%
                  filter(.data[[names(.)[[6]]]] == cur_aliquot) %>%
                  select("MD5") %>%
                  pull(1)
      return(tibble(case=new_name, R1_filename=r1_file_name, R1_uuid=r1_uuid,
                  R1_md5=r1_md5, R2_filename=r2_file_name,
                  R2_uuid=r2_uuid, R2_md5=r2_md5))
   }

   case_id <- read_tsv("case_id.txt") %>%
              pull(1)
   case_start <- ifelse($start == -1, 1, $start)
   case_end <- ifelse($end == -1, length(case_id), $end)
   stopifnot(case_start >= 1, case_start <= length(case_id))
   stopifnot(case_end >= 1, case_end <= length(case_id))
   stopifnot(case_start <= case_end)

   case_id <- case_id[case_start:case_end]

   catalog <- read_tsv("catalog.txt")
   catalog_rnaseq <- catalog %>% 
                     filter(case %in% case_id) %>%
                     select(sample_name = 1, everything()) %>%
                     filter(experimental_strategy == "RNA-Seq") %>%
                     filter(data_format == "FASTQ") %>%
                     {
                        if("aliquot_annotation" %in% names(.)) filter(., is.na(aliquot_annotation)) else .
                     }


   sample_names <-  catalog_rnaseq %>% pull(1)
   all_case_tbl <- tibble(case=character(),
                          R1_filename=character(),
                          R1_uuid=character(),
                          R1_md5=character(),
                          R2_filename=character(),
                          R2_uuid=character(),
                          R2_md5=character()) 
   for (case in case_id) {
      tmp <- str_match(sample_names, paste0("(",case,")\\\\.RNA-Seq\\\\.(R[1-2]).([^\\\\.]+)(.*)"))
      matched_case <- as_tibble(tmp, .name_repair = "unique") %>% 
                    filter(.data[[names(.)[[5]]]] == "") %>%
                    drop_na()
      # e.g. "A", "T"
      sample_type <- unique(matched_case[[4]])
      for (st in sample_type){
          r1_name <- matched_case %>% 
                     filter(.data[[names(.)[[3]]]] == "R1") %>%
                     filter(.data[[names(.)[[4]]]] == st)  %>%
                     pull(1)
          case_aliquot <- catalog_rnaseq %>%
                          select(sample_name, case, aliquot) %>%
                          filter(sample_name == unique(r1_name)) %>%
                          arrange(aliquot)
          if (length(r1_name) > 1) {
            # if allow duplicates, 
            # give new names ordered by aliquot id (alphabetically)
            # e.g. C3N-01525_T_1, C3N-01525_T_2, ...
            # with aliquot ID:  CPT95870043, CPT95870046, ...
            if("${params.allow_dup}" == "false") {
              print("error!!! duplicates found!!! ")
              print(paste0(length(r1_name), " duplicates for ", unique(r1_name), " found"))
              quit(status=1)
            } else {
              # create multiple r1_name, one for each aliquot id
              temp_info <- tibble(case=character(), R1_filename=character(),
                            R1_uuid=character(), R1_md5=character(),
                            R2_filename=character(),
                            R2_uuid=character(), R2_md5=character())
              for (i in seq_len(nrow(case_aliquot))){
                 res <- get_file_info(catalog_rnaseq, case, st, 
                                      unique(r1_name), i,
                                      matched_case, case_aliquot) 
                 temp_info <- temp_info %>%
                              bind_rows(res)    
              }
            }
          }
          else{  # no duplication
            temp_info <- get_file_info(catalog_rnaseq, case, st, unique(r1_name),
                                        0, matched_case, case_aliquot)  
          }
          all_case_tbl <- all_case_tbl %>%
                          bind_rows(temp_info)
          
      }
   }

   # write_tsv(all_case_tbl, "all_case_info.tsv")
   # case_id   r1_uuid    r2_uuid
   for (i in seq(1, nrow(all_case_tbl))) {
      write_tsv(all_case_tbl[i, c("case", "R1_uuid", "R1_filename", "R1_md5", 
                "R2_uuid", "R2_filename", "R2_md5")], paste0('case_', i, '.tsv'),
                col_names = FALSE)
   }
   """
}


process generate_id_files {
   container "${params.container.r_tidyverse}"
   label 'r5_2xlarge'

   input:
     path('catalog.txt')
     val start
     val end

   output:
     path 'case_*.tsv'

   script:
   """
   #!/usr/bin/env Rscript
   library(tidyverse)
   
   catalog <- read_tsv("catalog.txt")
   case_id <- catalog %>% pull(1)
   case_start <- ifelse($start == -1, 1, $start)
   case_end <- ifelse($end == -1, length(case_id), $end)
   stopifnot(case_start >= 1, case_start <= length(case_id))
   stopifnot(case_end >= 1, case_end <= length(case_id))
   stopifnot(case_start <= case_end)

   case_id <- case_id[case_start:case_end]
   catalog_rnaseq <- catalog %>% 
                     filter(sample_name %in% case_id) 

   # case_id   R1_path    R2_path
   for (i in seq(1, nrow(catalog_rnaseq))) {
      write_tsv(catalog_rnaseq[i,], paste0('case_', i, '.tsv'),
                col_names = FALSE)
   }
   """
}


process extract_paths {
  container "${params.container.ubuntu}"
  cpus 8
  memory '60 GB'
  label 'r5_2xlarge'

  input:
    path(id_file)

  output:
    env case_id, emit: case_id
    env r1_path, emit: r1_path
    env r2_path, emit: r2_path

  """
   while IFS=\$'\\t' read -r -a myid
   do
      case_id="\${myid[0]}"
      r1_path="\${myid[1]}"
      r2_path="\${myid[2]}"
   done < $id_file
   """
}


process download_files {
  container "${params.container.gdc_client}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    path(token)
    path(id_file)

 output:
  tuple path('CASE_ID'),
        path('R1/*.fastq'),
        path('R2/*.fastq'),
        emit: res_ch

  """
   while IFS=\$'\\t' read -r -a myid
   do
      case_id="\${myid[0]}"
      gdc-client download "\${myid[1]}" -t "$token" -n ${task.cpus}
      mv "\${myid[1]}" R1
      md5_1="\${myid[3]}"
      md5_result=\$(md5sum R1/\${myid[2]} | cut -d " " -f1)
      if [[ "\$md5_result" != "\$md5_1" ]]
       then
        echo "\${myid[2]} md5 does not match"
        exit 1
       else
        echo "md5 check OK"
      fi
      # unzip 
      r1_fastq="R1/\${case_id}_R1.fastq.gz"
      mv "R1/\${myid[2]}" "\${r1_fastq}"
      gzip -d "\${r1_fastq}"
      gdc-client download "\${myid[4]}" -t "$token" -n ${task.cpus}
      mv "\${myid[4]}" R2
      md5_2="\${myid[6]}"
      md5_result=\$(md5sum R2/\${myid[5]} | cut -d " " -f1)
      if [[ "\$md5_result" != "\$md5_2" ]]
      then
        echo "\${myid[5]} md5 does not match"
        exit 1
       else
        echo "md5 check OK"
      fi
      # unzip 
      r2_fastq="R2/\${case_id}_R2.fastq.gz"
      mv "R2/\${myid[5]}" "\${r2_fastq}"
      gzip -d "\${r2_fastq}"
      echo "\${case_id}" > CASE_ID
   done < $id_file
   """
}


process stage_files {
  container "${params.container.ubuntu}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    val(case_id)
    val(r1_path)
    val(r2_path)
    path(r1_path)
    path(r2_path)

 output:
  tuple path('CASE_ID'),
        path('R1/*.gz'),
        path('R2/*.gz'),
        emit: res_ch

  """
  echo "${case_id}" > CASE_ID
  mkdir R1 R2
  r1_name=\$(basename ${r1_path})
  mv \${r1_name} R1
  r2_name=\$(basename ${r2_path})
  mv \${r2_name} R2
  """
}


// step 2 
process fastq_to_sam {
  container "${params.container.bwa}"
  cpus 8
  memory '60 GB'
  label 'r5_2xlarge_500g'

  input:
    tuple path('CASE_ID'),
          path('r1.fastq.gz'),
          path('r2.fastq.gz')
    val(genome_ref_prefix)
    path('*')
    path('*')

  output:
    tuple path('CASE_ID'), 
          path('bwa-pe-for-CIRI.sam'),
          emit: res_ch

  """
   bwa mem -t ${task.cpus} -T 19 ${genome_ref_prefix} \
      r1.fastq.gz r2.fastq.gz -o bwa-pe-for-CIRI.sam 
  """

}


// step 3
process ciri_calling {
  publishDir params.outdir_ciri, 
             pattern: '*/results_CIRI.txt',
             mode: 'copy', overwrite: true 
  container "${params.container.ciri}"
  label 'r5_4xlarge_500g'
  cpus 16
  memory '124 GB'

  input:
    tuple path('CASE_ID'), path('bwa-pe-for-CIRI.sam')
    val(genome_ref_prefix)
    path('*')
    path('*')
    path('anno.gtf')

  output:
    tuple path('CASE_ID'), path("*/results_CIRI.txt"), emit: res_ch

  """
  case_id=\$(cat CASE_ID) 
  mkdir \${case_id} 
  perl /usr/src/CIRI2.pl -T ${task.cpus} \
  -I bwa-pe-for-CIRI.sam \
  -O \${case_id}/results_CIRI.txt \
  -F ${genome_ref_prefix} \
  -A anno.gtf 
  """
}


// step 4
process add_gene_name_to_CIRI_results {
  publishDir params.outdir_ciri,
             pattern: '*/results_CIRI_added_gene_name.txt',
             mode: 'copy', overwrite: true
  container "${params.container.rna_seq_misc_ydou}"
  label 'r5_2xlarge'

  input:
    tuple path('CASE_ID'), path('results_CIRI.txt') 
    path('ref_file.gtf')
  
  output:
    tuple path('CASE_ID'), path("*/results_CIRI_added_gene_name.txt"), emit: res_ch
  
  """
  case_id=\$(cat CASE_ID) 
  mkdir \${case_id} 
  1_add_gene_name_to_CIRI_results.pl -I results_CIRI.txt -F ref_file.gtf \
  -O \${case_id}/results_CIRI_added_gene_name.txt
  """

}


// step 5
process add_linear_circular_isoform_to_gtf{
  container "${params.container.rna_seq_misc_ydou}"
  label 'r5_2xlarge'

  input:
    tuple path('CASE_ID'), path('results_CIRI_added_gene_name.txt')
    path('ref_file.gtf')
  
  output:
    tuple path('CASE_ID'),
          path("*/genome_with_circular_RNA.gtf"), emit: res_ch
  
  """
  case_id=\$(cat CASE_ID) 
  mkdir \${case_id}
  2_add_linear_circular_isoform_to_gtf.pl -I results_CIRI_added_gene_name.txt -F ref_file.gtf \
  -O \${case_id}/genome_with_circular_RNA.gtf
  """
}


// step 6
process extract_linear_and_circRNA_transcripts {
  // publishDir params.outdir_genome,  
  //            pattern: "*/${genome_ref_prefix}.linear.and.circrna.*",
  //            mode: 'copy', overwrite: true
  container "${params.container.rsem}"
  label 'r5_2xlarge'

  input:
    tuple path('CASE_ID'), path('anno.gtf')
    val(genome_ref_prefix)
    path('*')
    path('*')

  output:
    tuple path('CASE_ID'),
          path("*/${genome_ref_prefix}.linear.and.circrna.transcripts.fa"), emit: res_ch
          path "*/${genome_ref_prefix}.linear.and.circrna.*", emit: output_ch
     
"""
case_id=\$(cat CASE_ID) 
mkdir \${case_id}
/opt/RSEM/rsem-extract-reference-transcripts \${case_id}/${genome_ref_prefix}.linear.and.circrna \
0 anno.gtf None 0 ${genome_ref_prefix}
"""
}


// step 7
process circular_linear_to_psedo_linear{
  // publishDir params.outdir_genome,
  //            pattern: '*/linear.and.circrna.as.pseudo.linear.transcripts.fa',
  //            mode: 'copy', overwrite: true
  container  "${params.container.rna_seq_misc_ydou}"
  label 'r5_2xlarge'

  input:
    tuple path('CASE_ID'), path('linear.and.circrna.transcripts.fa')

  output:
    tuple path('CASE_ID'), path("*/linear.and.circrna.as.pseudo.linear.transcripts.fa"), emit: res_ch

  """
  case_id=\$(cat CASE_ID) 
  mkdir \${case_id}
  3_circular_linear_to_psedo_linear.pl -I linear.and.circrna.transcripts.fa \
  -O \${case_id}/linear.and.circrna.as.pseudo.linear.transcripts.fa
  """
}


// step 8
process gene_isoform_mapping_with_circular_rna{
  // publishDir params.outdir_genome,
  //            pattern: '*/gene_isoform_mapping_with_circular_rna_for_RSEM.txt',
  //            mode: 'copy',
  //            overwrite: true
  container  "${params.container.rna_seq_misc_ydou}"
  label 'r5_2xlarge'

  input:
    tuple path('CASE_ID'), path('genome_with_circular_RNA.gtf')

  output:
    tuple path('CASE_ID'), 
          path('*/gene_isoform_mapping_with_circular_rna_for_RSEM.txt'),
          emit: res_ch

 """
  case_id=\$(cat CASE_ID) 
  mkdir \${case_id}
  4_gene_isoform_mapping_with_circular_rna.pl \
  -I genome_with_circular_RNA.gtf \
  -O \${case_id}/gene_isoform_mapping_with_circular_rna_for_RSEM.txt
 """

}


// step 9
process build_rsem_index {
  // publishDir params.outdir_genome, mode: 'copy', overwrite: true
  container "${params.container.rsem}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    tuple val(case_id), 
          path('linear.and.circrna.as.pseudo.linear.transcripts.fa'),
          path('gene_isoform_mapping_with_circular_rna_for_RSEM.txt')
    
  output: 
    tuple val(case_id),
          path("${case_id}/hg38*"),
          emit: res_ch

"""
mkdir ${case_id}
/opt/RSEM/rsem-prepare-reference -p ${task.cpus} \
--transcript-to-gene-map gene_isoform_mapping_with_circular_rna_for_RSEM.txt \
--bowtie2 --bowtie2-path /opt/bowtie2 \
linear.and.circrna.as.pseudo.linear.transcripts.fa ${case_id}/hg38 
"""
}


// step 10
process gene_and_transcript_quantification {
  publishDir params.outdir_rsem2, 
             pattern: "${case_id}/${case_id}.stat",
             mode: 'copy',
             overwrite: true
  publishDir params.outdir_rsem2, 
             pattern: "${case_id}/${case_id}.*.results",
             mode: 'copy',
             overwrite: true
  container "${params.container.rsem}"
  label 'r5_2xlarge_500g'
  cpus 8
  memory '60 GB'

  input:
    tuple val(case_id),
          path('r1.fastq.gz'),
          path('r2.fastq.gz'),
          path('*')
    
  output: 
    tuple val(case_id),
          path("${case_id}/${case_id}*"),
          emit: res_ch

"""
mkdir ${case_id}
/opt/RSEM/rsem-calculate-expression --bowtie2 --bowtie2-path \
/opt/bowtie2/ --paired-end -p ${task.cpus} \
r1.fastq.gz r2.fastq.gz hg38 ${case_id}/${case_id}
"""
}

// step 11
process summarize_gene_quantification {
 publishDir params.outdir_rsem,
             pattern: "${case_id}",
             mode: 'copy',
             overwrite: true
  container  "${params.container.rna_seq_misc_ydou}"
  label 'r5_2xlarge'
  cpus 4
  memory '30 GB'

  input:
    tuple val(case_id), path('*')

  output:
    path "${case_id}", emit: res_ch

  """
  mkdir ${case_id}
  5_summary-gene-quantification.pl
  mv *.txt ${case_id}
  """
}


// step 12
process combine_and_summary {
  publishDir params.outdir_rsem,
             pattern: "RSEM_results_summarized",
             mode: 'copy',
             overwrite: true
  container  "${params.container.rna_seq_combine_summary}"
  label 'r5_2xlarge'
  cpus 4
  memory '30 GB'

  input:
    path('*')

  output:
  // all results will be generated in the following folder
    path('RSEM_results_summarized')

  """
  Rscript /usr/src/combine_and_summary.r 
  """
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

assert params.run_version
if (params.data_source == 'gdc') {
  assert params.case_id
}

if (params.start > 0 && params.end > 0) {
  assert params.start <= params.end
}

if (workflow.profile.contains('awsbatch')) {
  if (!params.outdir.startsWith('s3:')) exit 1, "outdir not on S3 - specify S3 path to run on AWSBatch!"
}

workflow {
   log.info rnaseqHeader()
   // step 1
   if (params.run_indexing) {
      genome_index(params.genome_ref_prefix, params.genome_ref)
   }
   if (params.data_source == 'gdc') {
     generate_id_files_gdc(params.catalog_file, params.case_id,
                    params.start, params.end)
     download_files(params.gdc_token, generate_id_files_gdc.out.flatten())
     file_ch = download_files.out.res_ch
   } else { // s3 or local files will be staged
     generate_id_files(params.catalog_file, params.start, params.end)
     extract_paths(generate_id_files.out.flatten())
     stage_files(extract_paths.out.case_id,
                 extract_paths.out.r1_path,
                 extract_paths.out.r2_path,
                 extract_paths.out.r1_path,
                 extract_paths.out.r2_path
                 )
     file_ch = stage_files.out.res_ch
   }
   // step 2
   fastq_to_sam(file_ch,
        params.genome_ref_prefix,
        params.genome_ref,
        Channel.fromPath(params.genome_ref_index).collect())
   // step 3
   ciri_calling(fastq_to_sam.out.res_ch,
        params.genome_ref_prefix,
        params.genome_ref,
        Channel.fromPath(params.genome_ref_index).collect(),
        params.genome_ref_anno)
   // step 4
  add_gene_name_to_CIRI_results(ciri_calling.out.res_ch,
        params.genome_ref_anno) 
   // step 5
  add_linear_circular_isoform_to_gtf(add_gene_name_to_CIRI_results.out.res_ch,
        params.genome_ref_anno) 
   // step 6
  extract_linear_and_circRNA_transcripts(add_linear_circular_isoform_to_gtf.out.res_ch,
        params.genome_ref_prefix,
        params.genome_ref,
        Channel.fromPath(params.genome_ref_index).collect())
   // step 7
  circular_linear_to_psedo_linear(extract_linear_and_circRNA_transcripts.out.res_ch)
   // step 8
  gene_isoform_mapping_with_circular_rna(add_linear_circular_isoform_to_gtf.out.res_ch)
   // step 9
  ch_1 = circular_linear_to_psedo_linear.out.res_ch
  ch_2 = gene_isoform_mapping_with_circular_rna.out.res_ch
  ch1_new = ch_1.map{it -> [it[0].text.trim(), it[1]]}
  ch2_new = ch_2.map{it -> [it[0].text.trim(), it[1]]}
  new_ch = ch1_new.combine(ch2_new, by:0)
  build_rsem_index(new_ch)
   // step 10
  fastq_ch_new = file_ch.map{it -> [it[0].text.trim(), it[1], it[2]]}
  rsem_idx_ch = build_rsem_index.out.res_ch
  new_input_ch = fastq_ch_new.combine(rsem_idx_ch, by:0)
  gene_and_transcript_quantification(new_input_ch)
   // step 11
  summarize_gene_quantification(gene_and_transcript_quantification.out.res_ch)
   // step 12
  combine_and_summary(summarize_gene_quantification.out.res_ch.collect())
}

report_txt = """$header
 name = ${workflow.manifest.name}
 author = ${workflow.manifest.author}
 homePage = ${workflow.manifest.homePage}
 description = ${workflow.manifest.description}
 mainScript = ${workflow.manifest.mainScript}
 nextflowVersion = ${workflow.manifest.nextflowVersion}
 version = ${workflow.manifest.version}
""".stripIndent()


workflow.onComplete {
  def output_d = new File("results/pipeline_info/")
  if (!output_d.exists()) {
      output_d.mkdirs()
  }
  println "Pipeline completed at: $workflow.complete"
  def output_tf = new File(output_d, "pipeline_report.txt")
  output_tf.withWriter { w -> w << report_txt }
}