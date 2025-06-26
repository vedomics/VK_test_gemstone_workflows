version 1.0

workflow SKA_compare_samples {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
     }
  input {
    Array[String] samplenames
    Array[File] assembly_or_chromosome
    String straingst_strain
  }
  call ska2_build_to_distance {
    input:
        samplenames = samplenames,
        assembly_or_chromosome = assembly_or_chromosome,
        strain = straingst_strain
    }
  
  output {
    File ska2_skf_distances_file = ska2_build_to_distance.skf_distances_file
    File ska2_descriptor_stats = ska2_build_to_distance.ska_nk_out
    File ska2_skf_file = ska2_build_to_distance.skf_file
    File ska2_snps_vcf = ska2_build_to_distance.snps_vcf
    String ska2_strain = ska2_build_to_distance.strain_name
  }
}


# TASKS #

task ska2_build_to_distance {

  input {
        Array[String] samplenames
        Array[File] assembly_or_chromosome
        String strain 
    }

  String skf_filelist = "all_skf_files.txt"
  
  command <<<

    fasta_array=(~{sep=" " assembly_or_chromosome})
    printf "%s\n" "${fasta_array[@]}" > fastas.txt

    names_array=(~{sep=" " samplenames})
    printf "%s\n" "${names_array[@]}" > names.txt

    paste names.txt fastas.txt > ~{strain}.distance.txt

    touch ~{strain}_ska_nk_out.txt
    touch seqs.skf
    touch ~{strain}_skalo_out_snps.vcf
    touch strain.txt

  >>>


  output {
        File skf_distances_file = "~{strain}.distance.txt"
        File ska_nk_out = "names.txt"
        File skf_file = "fastas.txt"
        #File ska_nk_out = "~{strain}_ska_nk_out.txt"
        #File skf_file = "seqs.skf"
        File snps_vcf = "~{strain}_skalo_out_snps.vcf"
        String strain_name = read_string("strain.txt")

  }
  
  runtime {
        docker:"ubuntu:latest"
        memory: "1 GB" # REMEMEBR TO CHANGE BACK TO 150
        disks: "local-disk 200 HDD"
  }
  
}
