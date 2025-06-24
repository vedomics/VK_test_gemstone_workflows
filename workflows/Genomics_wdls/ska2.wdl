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

  command <<<

   # create text file with filenames

      fasta_array=(~{sep=" " assembly_or_chromosome})
      for i in ${fasta_array[@]}; do echo $i >> fastas.txt ; done

      names_array=(~{sep=" " samplenames})
      printf "%s\n" "${names[@]}" > names.txt

      paste names.txt fastas.txt > ska_input_file.txt

      # Run SKA BUILD - generates skf file with all isolates  
      ska_build -o seqs -f ska_input_file.txt

      # Run SKA nk - generates characteristics of each isolate. Need to parse in subsequent analysis

      ska nk seqs.skf > ~{strain_name}_ska_nk_out.txt

      # Run SKA distance

      ska distance -o ~{strain}_distance seqs.skf

      # Run SKA lo

      ska lo seqs.skf ~{strain}_skalo_out

      # Report strain info if provided

      echo ~{strain} > strain.txt


  >>>


  output {
        File skf_distances_file = "~{strain}.distance.txt"
        File ska_nk_out = "~{strain}_ska_nk_out.txt"
        File skf_file = "seqs.skf"
        File snps_vcf = "~{strain}_skalo_out_snps.vcf"
        String strain_name = read_string("strain.txt")

  }
  
  runtime {
      docker:"staphb/ska:latest"
        memory: "150 GB"
        disks: "local-disk 200 HDD"
  }
  
}
