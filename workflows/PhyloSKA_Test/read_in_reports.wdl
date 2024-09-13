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
	python3 <<< EOF
	import sys

	strains_to_cov = {}
			
	report_file = "~{straingst_report}"

	with open(report_file) as infile:
		next(infile) #skip headers 
		for line in infile:
			l=line.split('\t')
			strains_to_cov[l[1]] = l[5]

	output_strain = max(strains_to_cov, key = strains_to_cov.get)
	out_cov = strains_to_cov[max(strains_to_cov, key = strains_to_cov.get)]

	if out_cov > 0.8:
		print(output_strain)
	else:
		print("Insufficient_COV")

	EOF
	>>>
	output {
		String straingst_top_strain = stdout()
	}
	runtime{
		docker: "python:latest"
	}
}