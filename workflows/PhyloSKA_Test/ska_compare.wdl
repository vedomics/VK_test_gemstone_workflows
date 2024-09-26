version 1.0

workflow SKA_compare_samples {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
	}
  input {
    Array[File] skf_report
    Array[String] straingst_strain
    Float? Identity_cutoff
    Int? SNP_cutoff
    String? Filtering_Params_Used_to_Generate_Set
	}
  call ska_distance_matrix {
    input:
    	skf_report = skf_report,
    	strains = straingst_strain,
    	filter_params = Filtering_Params_Used_to_Generate_Set,
    	id_cutoff = Identity_cutoff,
    	snp_cutoff = SNP_cutoff
    }
  output {
    String filtering_params = ska_distance_matrix.params
    File skf_distances_file = ska_distance_matrix.skf_distances_file
    File skf_clusters_file = ska_distance_matrix.skf_clusters_file
    File skf_dot_file = ska_distance_matrix.skf_dot_file
  }
}

# Tasks #


task ska_distance_matrix {
	input {
	    Array[File] skf_report
	    Array[String] strains
	    Int? snp_cutoff
	    Float? id_cutoff
        String? filter_params
	}

	Float identity_cutoff_actual = select_first([id_cutoff,0.9])
	Int snp_cutoff_actual = select_first([snp_cutoff,20])
	String filter_params_actual = select_first([filter_params, "NA"])
	String params_file = "params.txt"
	String skf_filelist = "all_skf_files.txt"
    String strain_name = strains[0]

	command <<<

            echo ~{filter_params_actual} > ~{params_file}
            skf_array=(~{sep=" " skf_report})
            for i in ${skf_array[@]}; do echo $i >> ~{skf_filelist}; done
            ska distance -f ~{skf_filelist} -i ~{identity_cutoff_actual} -s ~{snp_cutoff_actual} -o ~{strain_name}

	>>>

	output {
        File params = "params.txt"
        File skf_distances_file = "~{strain_name}.distances.tsv"
        File skf_clusters_file = "~{strain_name}.clusters.tsv"
        File skf_dot_file = "~{strain_name}.dot"

	}
	runtime {
	    docker:"staphb/ska:latest"
        memory: "64 GB"
        disks: "local-disk 100 HDD"
	}
}
