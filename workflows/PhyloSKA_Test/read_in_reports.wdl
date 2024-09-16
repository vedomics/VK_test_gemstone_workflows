version 1.0

workflow read_strain {
 meta {
    author: "Veda Khadka"
	}
	input {
		File straingst_report
	}
  call straingst_strain {
    input:
    	straingst_report = straingst_report
  }
  output {
  	String straingst_strain = reader.straingst_top_strain
  }
}

# Tasks #

task straingst_strain {	
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
