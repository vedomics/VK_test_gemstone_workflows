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
  Boolean generate_tree = false
  Float? min_kmer_freq

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

if (generate_tree) {
  call SKA_align {
    input:
        skf_files = SKA1_build.skf_file,
        strain = strain_name,
        kmer_freq = min_kmer_freq,
        params = SKA1_build.build_parameters
  }
  call build_tree {
    input:
      aligned_skf = SKA_align.aligned_skf,
      strain = strain_name,
      params = SKA1_build.build_parameters
  }
}

  if (generate_vcf){
   call SKA1_annotate {
        input:
          names = samplename,
          skf_files = SKA1_build.skf_file,
          ref = ref_genome,
          strain = strain_name,
          params = SKA1_build.build_parameters
      }
}

  output {

    File skf_summary = SKA1_distance.summaries
    File? ska_vcfs = SKA1_annotate.vcfs
    File? ska_tree = build_tree.treefile
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
      mem: "1 GB"
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


# Optional tasks

task SKA1_annotate {

   input {
      Array[String] names
      Array[File] skf_files
      File? ref
      Array[String] params
      String strain
  
  }

    String user_params = params[0]
    String skf_distances_named = "~{strain}_~{user_params}"

  command <<<

            skf_array=(~{sep=" " skf_files})
            names_array=(~{sep=" " names})

            mkdir vcf_files

            for index in ${!skf_array[*]}; do 
              ska annotate -r ~{ref} -o ${names_array[$index]} ${skf_array[$index]}
              head -7 ${names_array[$index]}.vcf > ${names_array[$index]}_filt.vcf && grep 'NS5' ${names_array[$index]}.vcf >> ${names_array[$index]}_filt.vcf
              sed -i '' '8d' ${names_array[$index]}_filt.vcf
              mv ${names_array[$index]}_filt.vcf vcf_files/
            done

             # Generate vcf tarball

            tar -czf ~{skf_distances_named}_vcf.tar.gz vcf_files

        
  >>>

  output {

    File vcfs = "~{skf_distances_named}_vcf.tar.gz"

  }

  runtime {
        docker:"staphb/ska:latest"
        memory: "100 GB"
        disks: "local-disk 200 HDD"
  }

}

task SKA_align {
   input {
      Array[File] skf_files
      Array[String] params
      String strain
      Float? kmer_freq
  
  }

    String user_params = params[0]
    String skf_distances_named = "~{strain}_~{user_params}"
    Float kfreq = select_first([kmer_freq,0.9])

  command <<<

            skf_array=(~{sep=" " skf_files})
            ska merge -o ~{skf_distances_named}_merged ${skf_array[@]}
            ska align -v -p ~{kfreq} -o ~{skf_distances_named} ~{skf_distances_named}_merged.skf

  >>>

  output {

    File aligned_skf = "~{skf_distances_named}_variants.aln"

  }

  runtime {
        docker:"staphb/ska:latest"
        memory: "50 GB"
        disks: "local-disk 200 HDD"
  }

}

task build_tree {
  
input {
      Array[String] params
      String strain
      File aligned_skf
  
  }

    String user_params = params[0]
    String skf_distances_named = "~{strain}_~{user_params}"


  command <<<

    VeryFastTree -nt -gamma -gtr -threads 4 ~{aligned_skf} > ~{skf_distances_named}.tre


  >>> 

  output {

    File treefile = "~{skf_distances_named}.tre"

  }

  runtime {
        docker:"vkhadka/veryfasttree:v4.0.5"
        memory: "5 GB"
        cpu: 4
        disks: "local-disk 200 HDD"
  }

}
