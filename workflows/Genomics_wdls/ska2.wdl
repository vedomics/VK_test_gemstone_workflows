version 1.0

workflow SKA_compare_samples {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
     }
  input {
    Array[String] samples
    Array[File] assembly_or_chromosome
    String straingst_strain
    File reference
    Float? min_freq
    Int? kmer_size
  }
  call ska2_build_to_distance {
    input:
        samplenames = samples,
        assembly_or_chromosome = assembly_or_chromosome,
        strain = straingst_strain,
        ref = reference,
        minfreq = min_freq,
        kmers = kmer_size
    }
  
  output {
    File ska2_skf_distances_file = ska2_build_to_distance.skf_distances_file
    File ska2_descriptor_stats = ska2_build_to_distance.ska_nk_out
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
        File ref 
        Float? minfreq
        Int? kmers
    }

  String skf_filelist = "all_skf_files.txt"
  Float minfreq_actual = select_first([minfreq,0.9])
  Int kmers_actual = select_first([kmers,31])

  command <<<

            printf "%s\n" "~{sep="\n" samplenames}" >> names.txt
            printf "%s\n" "~{sep="\n" assembly_or_chromosome}" >> fastas.txt

            paste names.txt fastas.txt > ska_input_file.txt

             # Run SKA BUILD - generates skf file with all isolates  
            
            ska build -o seqs -k ~{kmers_actual} -f ska_input_file.txt

            # Run SKA nk - generates characteristics of each isolate. Need to parse in subsequent analysis

            ska nk seqs.skf > ~{strain}_ska_nk_out.txt

            # Run SKA distance

            ska distance -o ~{strain}_distance.txt --min-freq ~{minfreq_actual} seqs.skf

            # Run SKA lo

            ska lo seqs.skf ~{strain}_skalo_out -r ~{ref}

            # Report strain info if provided

            echo ~{strain} > strain.txt

  >>>

  output {
        File skf_distances_file = "~{strain}_distance.txt"
        File ska_nk_out = "~{strain}_ska_nk_out.txt"
        File snps_vcf = "~{strain}_skalo_out_snps.vcf"
        String strain_name = "~{strain}"
  }
  
  runtime {
        docker:"vkhadka/ska2:v0.4.1"
        memory: "150 GB" 
        disks: "local-disk 200 HDD"
        shell: "/bin/bash"
  }
  
}
