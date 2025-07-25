version 1.0

workflow test {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
     }

  call hello
}

task hello {
  command {
    echo 'Hello world!'
  }

  output {
    File response = stdout()
  }

  runtime {
    docker: "quay.io/ga4gh-dream/dockstore-tool-helloworld:1.0.2"
  }
}


