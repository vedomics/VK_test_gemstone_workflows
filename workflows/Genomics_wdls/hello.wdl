version 1.0

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

workflow test {
  call hello
}
