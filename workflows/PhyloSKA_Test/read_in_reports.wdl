version 1.0

workflow baby {
 meta {
    author: "Veda Suffers"
	}
	input {
		File straingst_report
	}
  call reader {
    input:
    	straingst_report = straingst_report
  }
  output {
  	String straingst_strain = reader.straingst_top_strain
  }
}

# Tasks #

task reader {	
	input {
		File straingst_report
	}
	command <<<
		python3 /app/read_tsv.py ~{straingst_report}
	>>>
	output {
		String straingst_top_strain = stdout()
	}
	runtime{
		docker: "vkhadka/reader-test@sha256:f3fc630bc4929cc8bc0cbbf06902c76b98729c940dac2b51a5753d293ffd8c7d"
	}
}
