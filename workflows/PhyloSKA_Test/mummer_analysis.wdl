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
      File mummer_plot = mummer.plot
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
    s1="$( echo ~{sample1})"
    s2="$( echo ~{sample2})"


  >>>
  output {
    File report = "out.report"
    File plot = "mummer.ps"
    File delta = "out.delta"
    String sample_1 = read_string("s1")
    String sample_2 = read_string("s2")
    
  }
  runtime {
    docker: "staphb/mummer:4.0.1"
    memory: "100 GB"
    cpu: 1
    preemptible: 0
    maxRetries: 3
  }
}