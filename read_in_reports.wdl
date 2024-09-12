version 1.0

workflow baby {
  call reader {
  }
  input:
   File straingst_report
}

# Tasks #

task reader {
}