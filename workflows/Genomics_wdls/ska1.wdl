version 1.0

workflow SKA_1 {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
  }

  input {
  Array[File] read1_clean
  Array[File] read2_clean
  Array[String] samplename
  File ref_genome
  String strain_name
  Float? minor_allele_freq
  Int? kmer_size
  Int? file_Coverage_cutoff
  Int? total_Coverage_cutoff
  Float? identity_cutoff
  Int? snp_cutoff
  }

  scatter (i in range(length(samplename))) {
    call SKA1_build {
      input:
      fq1 = read1_clean[i],
      fq2 = read2_clean[i],
      name = samplename[i],
      ref = ref_genome,
      minor_freq = minor_allele_freq,
      kmers = kmer_size,
      file_cutoff= file_Coverage_cutoff,
      total_cutoff = total_Coverage_cutoff
    }
  }


  call SKA1_distance {
    input:
      skf_files = SKA1_build.skf_file,
      skf_summary = SKA1_build.skf_summary,
      skf_vcf = SKA1_build.skf_vcf,
      strain = strain_name,
      params = SKA1_build.build_parameters,
      snp_cutoff = snp_cutoff,
      identity_cutoff = identity_cutoff
  }

  output {
    File skf_summary = SKA1_distance.summaries
    File ska_vcfs = SKA1_distance.vcfs
    File ska_distance = SKA1_distance.distance_matrix
    File ska_clusters = SKA1_distance.clusters

  }
}

#### Tasks ####

task SKA1_build {  
  input {
  File fq1
  File fq2
  String name
  File ref
  Float? minor_freq
  Int? kmers
  Int? file_cutoff
  Int? total_cutoff
  }

  # Tweakable parameters currently set to SKA defaults
  Float MAF = select_first([minor_freq,0.2])
  Int kmers_actual = select_first([kmers,15])
  Int file_cov = select_first([file_cutoff,4])
  Int total_cov = select_first([total_cutoff,2])
  String skf_summary = "~{name}_k15_summary.txt"
  String input_parameters = "~{kmers_actual}_~{MAF}_~{total_cov}_~{file_cov}"
  
  command <<<
            
            ska fastq -m ~{MAF} -k ~{kmers_actual} -c ~{total_cov} -C ~{file_cov} -o ~{name} ~{fq1} ~{fq2}
            ska summary ~{name}.skf > ~{skf_summary}
            ska annotate -r ~{ref} -o ~{name} ~{name}.skf
        
  >>>

  output {
    File skf_file = glob("*.skf")[0]
    File skf_summary = skf_summary
    File skf_vcf = "~{name}.vcf"
    String build_parameters = input_parameters
  }

  runtime {
      docker: "staphb/ska:latest"
      mem: "3 GB"
      cpu: 2
      preemptible: 0
      maxRetries: 3
  }
}


task SKA1_distance {  
  input {
    String strain
    Array[File] skf_files
    Array[File] skf_summary
    Array[File] skf_vcf
    Array[String] params
    Float? identity_cutoff
    Int? snp_cutoff
   }

  Float identity_cutoff_actual = select_first([identity_cutoff,0.9])
  Int snp_cutoff_actual = select_first([snp_cutoff,20])
  String skf_filelist = "all_skf_files.txt"
  String user_params = params[0]
  String skf_distances_named = "~{strain}_~{user_params}"


command <<<

        # Generate distances and clusters files

            skf_array=(~{sep=" " skf_files})
            for i in ${skf_array[@]}; do echo $i >> ~{skf_filelist}; done
            ska distance -f ~{skf_filelist} -i ~{identity_cutoff_actual} -s ~{snp_cutoff_actual} -o ~{skf_distances_named}

        # Generate vcf tarball

            touch all_vcf_files.txt
            mkdir vcf_files
            vcf_array=(~{sep=" " skf_vcf})
            for i in ${vcf_array[@]}; do cp $i vcf_files; done
            tar -czf ~{skf_distances_named}_vcf.tar.gz vcf_files


        # Generate summaries tarball

            touch all_summaries_files.txt
            summary_array=(~{sep=" " skf_summary})
            for i in ${summary_array[@]}; do cat $i >> all_summaries_files.txt; done
            sed '1!{/^Sample/d;}' all_summaries_files.txt > all_summaries_files_clean.txt
            tar -czf ~{skf_distances_named}_summaries.tar.gz -T all_summaries_files_clean.txt
  >>>


output {
        File distance_matrix = "~{skf_distances_named}.distances.tsv"
        File clusters = "~{skf_distances_named}.clusters.tsv"
        File vcfs = "~{skf_distances_named}_vcf.tar.gz"
        File summaries = "~{skf_distances_named}_summaries.tar.gz"

  }

   runtime {
        docker:"staphb/ska:latest"
        memory: "150 GB"
        disks: "local-disk 200 HDD"
  }

  }
