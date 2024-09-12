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
		python3 ../scripts/read_tsv.py ~{straingst_report}
	>>>
	output {
		String straingst_top_strain = read_string("STRAIN_REF")
	}
	runtime{
		docker: "3.12.6-bookworm⁠"
	}
}