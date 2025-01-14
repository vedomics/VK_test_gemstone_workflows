version 1.0

workflow bakta {
 
  meta {
    author: "Theiagen edited by VK"
    email: "vkhadka@broadinstitute.org"
  }

  input {
    File assembly
    File bakta_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_db_2022-08-29.tar.gz"
    String samplename
    Int cpu = 8
    Int memory = 16
    String docker = "quay.io/biocontainers/bakta:1.5.1--pyhdfd78af_0"
    Int disk_size = 100
    # Parameters 
    #  proteins: Fasta file of trusted protein sequences for CDS annotation
    #  prodigal_tf: Prodigal training file to use for CDS prediction
    # bakta_opts: any additional bakta arguments
    Boolean proteins = false
    Boolean compliant = false
    File? prodigal_tf
    String? bakta_opts
  }
  
  call bakta {
    input:
      assembly = assembly,
      samplename = samplename
    }

  output {
    File bakta_embl = bakta.bakta_embl
    File bakta_faa = bakta.bakta_faa
    File bakta_ffn = bakta.bakta_ffn
    File bakta_fna = bakta.bakta_fna
    File bakta_gbff = bakta.bakta_gbff
    File bakta_gff3 = bakta.bakta_gff3
    File bakta_hypotheticals_faa = bakta.bakta_hypotheticals_faa
    File bakta_hypotheticals_tsv = bakta.bakta_hypotheticals_tsv
    File bakta_tsv = bakta.bakta_tsv
    File bakta_txt = bakta.bakta_txt
  }

}

task bakta {
  input {
    File assembly
    File bakta_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_db_2022-08-29.tar.gz"
    String samplename
    Int cpu = 8
    Int memory = 16
    String docker = "quay.io/biocontainers/bakta:1.5.1--pyhdfd78af_0"
    Int disk_size = 100
    # Parameters 
    #  proteins: Fasta file of trusted protein sequences for CDS annotation
    #  prodigal_tf: Prodigal training file to use for CDS prediction
    # bakta_opts: any additional bakta arguments
    Boolean proteins = false
    Boolean compliant = false
    File? prodigal_tf
    String? bakta_opts
  }
  command <<<
  date | tee DATE
  
  # Extract Bakta DB
  mkdir db
  time tar xzvf ~{bakta_db} --strip-components=1 -C ./db

  # Install amrfinderplus db
  amrfinder_update --database db/amrfinderplus-db
  amrfinder --database_version | tee AMRFINDER_DATABASE_VERSION

  bakta \
    ~{bakta_opts} \
    --db db/ \
    --threads ~{cpu} \
    --prefix ~{samplename} \
    --output ~{samplename} \
    ~{true='--compliant' false='' compliant} \
    ~{true='--proteins' false='' proteins} \
    ~{'--prodigal-tf ' + prodigal_tf} \
    ~{assembly}
  
  # rename gff3 to gff for compatibility with downstream analysis (pirate)
  mv "~{samplename}/~{samplename}.gff3" "~{samplename}/~{samplename}.gff"
  
  >>>
  output {
    File bakta_embl = "~{samplename}.embl"
    File bakta_faa = "~{samplename}.faa"
    File bakta_ffn = "~{samplename}.ffn"
    File bakta_fna = "~{samplename}.fna"
    File bakta_gbff = "~{samplename}.gbff"
    File bakta_gff3 = "~{samplename}.gff"
    File bakta_hypotheticals_faa = "~{samplename}.hypotheticals.faa"
    File bakta_hypotheticals_tsv = "~{samplename}.hypotheticals.tsv"
    File bakta_tsv = "~{samplename}.tsv"
    File bakta_txt = "~{samplename}.txt"

  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: "quay.io/biocontainers/bakta:1.5.1--pyhdfd78af_0"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}