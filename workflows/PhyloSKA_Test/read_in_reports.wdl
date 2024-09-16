version 1.0

workflow read_strain {
 meta {
    author: "Veda Khadka"
	}
	input {
		File straingst_report
	}
  call read_straingst_report {
    input:
    	straingst_report = straingst_report
  }
  output {
  	String straingst_strain = read_straingst_report.straingst_top_strain
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
