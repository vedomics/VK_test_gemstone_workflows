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
  String strain_name
  Float? minor_allele_freq
  Int? kmer_size
  Int? file_Coverage_cutoff
  Int? total_Coverage_cutoff
  Float? identity_cutoff
  Int? snp_cutoff
  #Get genome positions
  Boolean generate_vcf = false
  File? ref_genome

  }

  scatter (i in range(length(samplename))) {
    call SKA1_build {
      input:
      fq1 = read1_clean[i],
      fq2 = read2_clean[i],
      name = samplename[i],
      minor_freq = minor_allele_freq,
      kmers = kmer_size,
      file_cutoff= file_Coverage_cutoff,
      total_cutoff = total_Coverage_cutoff
    }

    if (generate_vcf) {
      call SKA1_annotate {
        input:
          name = samplename[i],
          skf_file = SKA1_build.skf_file[i],
          ref = ref_genome
      }
    }
  }
  # If this fails, then accept that you'll just have to do a loop inside the command

  if (generate_vcf){
  call SKA1_vcf {
     input:
          vcf_file = SKA1_annotate.skf_vcf,
          strain = strain_name,
          params = SKA1_build.build_parameters

  }
}

  call SKA1_distance {
    input:
      skf_files = SKA1_build.skf_file,
      skf_summary = SKA1_build.skf_summary,
      strain = strain_name,
      params = SKA1_build.build_parameters,
      snp_cutoff = snp_cutoff,
      identity_cutoff = identity_cutoff
  }




  output {
    File skf_summary = SKA1_distance.summaries
    File? ska_vcfs = SKA1_vcf.vcfs
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
        
  >>>

  output {
    File skf_file = glob("*.skf")[0]
    File skf_summary = skf_summary
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

task SKA1_annotate {

   input {
      String name
      File skf_file
      File? ref
  
  }

      # Tweakable parameters currently set to SKA defaults
     
  
  command <<<

            ska annotate -r ~{ref} -o ~{name} ~{name}.skf
        
  >>>

  output {

    File skf_vcf = "~{name}.vcf"

  }

  runtime {
      docker: "staphb/ska:latest"
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


        # Generate summaries tarball

            touch all_summaries_files.txt
            summary_array=(~{sep=" " skf_summary})
            for i in ${summary_array[@]}; do cat $i >> all_summaries_files.txt; done
            sed '1!{/^Sample/d;}' all_summaries_files.txt > all_summaries_files_clean.txt
            mv all_summaries_files_clean.txt ~{skf_distances_named}_summaries.txt

  >>>


output {
        File distance_matrix = "~{skf_distances_named}.distances.tsv"
        File clusters = "~{skf_distances_named}.clusters.tsv"
        File summaries = " ~{skf_distances_named}_summaries.txt"

  }

   runtime {
        docker:"staphb/ska:latest"
        memory: "150 GB"
        disks: "local-disk 200 HDD"
  }

  }



task SKA1_vcf {

 input {
    String strain
    Array[File] vcf_file
    Array[String] params
   }

  String user_params = params[0]
  String skf_distances_named = "~{strain}_~{user_params}"

  command <<<


        # Generate vcf tarball

            touch all_vcf_files.txt
            mkdir vcf_files
            vcf_array=(~{sep=" " vcf_file})
            for i in ${vcf_array[@]}; do cp $i vcf_files; done
            tar -czf ~{skf_distances_named}_vcf.tar.gz vcf_files


  >>>

output {

        File vcfs = "~{skf_distances_named}_vcf.tar.gz"

  }

   runtime {
        docker:"staphb/ska:latest"
        memory: "150 GB"
        disks: "local-disk 200 HDD"
  }


  }
