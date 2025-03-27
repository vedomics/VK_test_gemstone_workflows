version 1.0

workflow mummer {
	meta {
		author: "Veda Khadka"
		email: "vkhadka@broadinstitute.org"
	}
	input {
		File assembly_1
		File assembly_2
		String samplename_1
		String samplename_2
	}

	call mummer {
		input:
			assembly_1 = assembly_1,
			assembly_2 = assembly_2,
			sample1 = samplename_1,
			sample2 = samplename_2
	}


	output {
		File mummer_report = mummer.report
	    File mummer_rpdf = mummer.rplot
	    File mummer_fpdf = mummer.fplot
	    File mummer_delta = mummer.delta
	}
}

#TASKS#

task mummer {
  input {
    File assembly_1
    File assembly_2
    String sample1
    String sample2
  }

   command <<<

    mummer -mum -b -c ~{assembly_1} ~{assembly_2} > mummmer.mums
    mummerplot -postscript -p mummer mummmer.mums
    dnadiff ~{assembly_1} ~{assembly_2}


  >>>
  output {
    File report = "out.report"
    File rplot = "mummer.rplot"
    File fplot = "mummer.fplot"
    File delta = "out.delta"
    
  }
  runtime {
    docker: "staphb/mummer:4.0.1"
    memory: "2 GB"
    cpu: 1
    preemptible: 0
    maxRetries: 3
  }
}