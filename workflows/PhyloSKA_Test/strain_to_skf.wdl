version 1.0


workflow SKA_dists {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
	}
  input {
    File straingst_report
	File fastq_read1
	File fastq_read2
	String samplename
	}
  call read_straingst_report {
    input:
    	straingst_report = straingst_report
    }
  call SKA_fastq {
    input:
        fq1 = fastq_read1,
        fq2 = fastq_read2,
        samplename = samplename
  }

  output {
    String straingst_strain = read_straingst_report.straingst_top_strain
    File skf_file = SKA_fastq.skf_file
    File skf_summary = SKA_fastq.skf_summary
  }
}

# Tasks #

task read_straingst_report {	
	input {
		File straingst_report
	}
	command <<<
		python3 /app/read_tsv.py ~{straingst_report}
	>>>
	output {
		String straingst_top_strain = read_string("STRAIN_REF")
	}
	runtime{
		docker: "vkhadka/reader-test:multi_V1"
	}
}

task SKA_fastq {	
	input {
		File fq1
		File fq2
		String samplename
	}
	# Tweakable parameters currently set to SKA defaults
	Int kmer_size = 15
	Float minor_allele_freq = 0.3
	# If "File" type is used Cromwell attempts to localize it, which fails because it doesn't exist yet.
 	String skf_path = "~{samplename}_k15"
	String skf_summary = "~{samplename}_k15_summary"
	command <<<
		ska fastq -m ~{minor_allele_freq} -k ~{kmer_size} -o ~{skf_path} ~{fq1} ~{fq2}
		ska summary ~{skf_path}.skf > ~{skf_summary}
	>>>
	output {
		File skf_file = glob("*.skf")[0]
		File skf_summary = skf_summary
	}
	runtime {
	    docker: "staphb/ska:latest"
	}
}


