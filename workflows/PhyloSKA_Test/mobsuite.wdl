version 1.0

workflow mobsuite {
  meta {
    author: "Veda Khadka"
    email: "vkhadka@broadinstitute.org"
    description: "Run MOBsuite on isolates, adapted from MT"
  }
  
  input {
    File assembly
    String samplename
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
    }
  
  call mob_recon {
    input:
      assembly = assembly,
      samplename = samplename
    }

  output {
    File mob_recon_results = mob_recon.mob_recon_results
    File mob_typer_results = mob_recon.mob_typer_results
    File chromosome_fasta = mob_recon.chromosome_fasta
    File plasmid_gz = mob_recon.plasmid_gz
    String mob_recon_version = mob_recon.mob_recon_version
   }
  }

### Tasks ###

task mob_recon {
  input {
    File assembly
    String samplename
    String docker = "kbessonov/mob_suite:3.0.3"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
  }
  command <<<
    # version capture
    mob_recon --version | cut -d ' ' -f2 > VERSION.txt

    mkdir mob_recon

    # unzip assembly FASTA
    recompress=false
    assembly=~{assembly}
    if [[ ~{assembly} == *.gz ]]; then 
      gzip -d ~{assembly}
      recompress=true
      assembly=${assembly: 0:-3}
    fi

    # run mob-recon
    mob_recon --infile $assembly --outdir mob_recon/~{samplename}

    # If the assembly FASTA was originally compressed, compress it again
    if $recompress; then
      gzip $assembly
    fi

    # Compress output FASTAs
    for fasta in mob_recon/~{samplename}/*.fasta; do 
      gzip $fasta
    done

    # Gzip all plasmid fastas together with plasmid name appended to each contig header

    touch mob_recon/plasmids.fasta
    for plasmid in mob_recon/~{samplename}/plasmid*.fasta.gz; do
      NAME=$(basename ${plasmid} ".fasta.gz")
      zcat $plasmid| sed "/^>.*/ s/$/_$NAME/" >> mob_recon/plasmids.fasta
    done

    tar -czvf plasmid_fastas.tar.gz mob_recon/plasmids.fasta


    # If the mobtyper report file does not exist, create it so that 
    # the output link does not point to an endless void
    if [ ! -f mob_recon/~{samplename}/mobtyper_results.txt ]; then
      touch mob_recon/~{samplename}/mobtyper_results.txt
    fi

  >>>
  output {
    File mob_recon_results = "mob_recon/~{samplename}/contig_report.txt"
    File mob_typer_results = "mob_recon/~{samplename}/mobtyper_results.txt"
    File chromosome_fasta = "mob_recon/~{samplename}/chromosome.fasta.gz"
    File plasmid_gz = "plasmid_fastas.tar.gz"
    String mob_recon_version = read_string("VERSION.txt")
    String mob_recon_docker = "~{docker}"
  }
  runtime {
    docker: "kbessonov/mob_suite:3.0.3"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}