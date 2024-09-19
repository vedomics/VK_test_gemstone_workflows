version 1.0

workflow read_strain {
 meta {
    author: "Veda Khadka"
	}
	input {
		File straingst_report
		String QC_check
		Float? coverage_cutoff
	}
  	call read_straingst_report {
	    input:
	    	straingst_report = straingst_report,
	    	covg_cutoff = coverage_cutoff
	  }
  output {
  	String straingst_top_strain = read_straingst_report.straingst_top_strain
  }
}

# Tasks #

task read_straingst_report {	
	input {
		File straingst_report
		Float? covg_cutoff
	}
	Float covg_cutoff_actual = select_first([covg_cutoff,0.8])
	command <<<
		python3 /app/read_tsv.py ~{straingst_report} ~{covg_cutoff_actual}
	>>>
	output {
		String straingst_top_strain = read_string("STRAIN_REF")
	}
	runtime{
		docker: "vkhadka/reader-test:multi_V1"
	}
}
