version 1.0

workflow tutorial {
	meta {
		author: "Veda Khadka"
		email: "vkhadka@broadinstitute.org"
	}
	input {
		File fastqfile
		String samplename
	}

	call fastp_trim {
		input:
			examplefile = fastqfile,
			samplename = samplename
	}

	call readcount {
	input:
			clean_fq1 = fastp_trim.clean_fq1,
			clean_fq2 = fastp_trim.clean_fq2
	}

	output {
		File read1_trimmed = fastp_trim.clean_fq1
	    File read2_trimmed = fastp_trim.clean_fq2
	    File fastp_stats = fastp_trim.fastp_stats
	    Float read1_clean_reads = readcount.reads1
	    Float read2_clean_reads = readcount.reads2
	}
}

# Tasks # 

task fastp_trim {
  input {
    File examplefile
    String samplename
  }
  
  command <<<

    fastp -i ~{examplefile} -o ~{samplename}_1.fastq -O ~{samplename}_2.fastq --interleaved_in

  >>>
  output {
    File clean_fq1 = "~{samplename}_1.fastq"
    File clean_fq2 = "~{samplename}_2.fastq"
    File fastp_stats = "fastp.html"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/staphb/fastp:0.23.2"
    memory: "2 GB"
    cpu: 1
    preemptible: 0
    maxRetries: 3
  }
}

task readcount {
  input {
    File clean_fq1
	File clean_fq2
	}

  command <<<

    read1_reads="$(grep -c "^@" ~{clean_fq1})"
    read2_reads="$(grep -c "^@" ~{clean_fq2})"

  >>>
  output {
    Float reads1 = read_string("read1_reads")
    Float reads2 = read_string("read2_reads")
  }
  runtime {
    docker: "ubuntu:latest"
    memory: "2 GB"
    cpu: 1
    preemptible: 0
    maxRetries: 3
  }
}

