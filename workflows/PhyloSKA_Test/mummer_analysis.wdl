version 1.0

workflow mummer {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
  }
  input {
    Array[File] comparison_assemblies
    Array[String] samplenames
  }

  call mummer {
    input:
    samples = samplenames,
    files = comparison_assemblies

  }


  output {
    File mummer_report = mummer.report
      File mummer_rpdf = mummer.rplot
      File mummer_fpdf = mummer.fplot
      File mummer_delta = mummer.delta
      String sample1 = mummer.sample_1
      String sample2 = mummer.sample_2

  }
}

#TASKS#

task mummer {
  input {
    Array[File] files
    Array[String] samples
  }

  String sample1 = samples[0]
  String sample2 = samples[1]

   command <<<

    mummer -mum -b -c ~{sep=" " files} > mummmer.mums
    mummerplot -postscript -p mummer mummmer.mums
    dnadiff ~{sep=" " files}


  >>>
  output {
    File report = "out.report"
    File rplot = "mummer.rplot"
    File fplot = "mummer.fplot"
    File delta = "out.delta"
    String sample_1 = read_string(~{sample1})
    String sample_2 = read_string(~{sample2})
    
  }
  runtime {
    docker: "staphb/mummer:4.0.1"
    memory: "2 GB"
    cpu: 1
    preemptible: 0
    maxRetries: 3
  }
}