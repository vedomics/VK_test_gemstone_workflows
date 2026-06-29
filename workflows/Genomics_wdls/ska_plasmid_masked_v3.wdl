version 1.0

workflow SKA_plasmid_masked_v3 {
  meta {
    author: "Veda Khadka & Claude Code"
    email: "vkhadka@broadinstitute.org"
    description: "SKA analysis with plasmid masking using pre-existing plasmid fastas - per-isolate plasmid merging"
  }

  input {
    Array[File] read1_clean
    Array[File] read2_clean
    File plasmid_fastas
    Array[String] samplename
    String strain_name
    Float? minor_allele_freq
    Int? kmer_size
    Int? file_Coverage_cutoff
    Int? total_Coverage_cutoff
    Float? identity_cutoff
    Int? snp_cutoff
    Boolean generate_vcf = false
    File? ref_genome
    Boolean generate_tree = false
    Float? min_kmer_freq = 0.9
  }

  # Extract plasmids from single tarball first
  call extract_all_plasmids {
    input:
      plasmid_tarball = plasmid_fastas,
      sample_names = samplename
  }

  scatter (i in range(length(samplename))) {
    call SKA1_build {
      input:
        fq1 = read1_clean[i],
        fq2 = read2_clean[i],
        name = samplename[i],
        minor_freq = minor_allele_freq,
        kmers = kmer_size,
        file_cutoff = file_Coverage_cutoff,
        total_cutoff = total_Coverage_cutoff
    }

    # Convert plasmid fasta to skf for this sample
    call plasmid_to_skf {
      input:
        plasmid_fasta = extract_all_plasmids.sample_plasmid_fastas[i],
        sample_id = samplename[i]
    }

    # Mask each genome with its own plasmids
    call mask_genomes {
      input:
        unmasked_skf = SKA1_build.skf_file,
        plasmid_skf = plasmid_to_skf.plasmid_skf,
        sample_id = samplename[i]
    }
  }

  call SKA1_distance_masked {
    input:
      skf_files = mask_genomes.masked_skf,
      skf_summary = SKA1_build.skf_summary,
      strain = strain_name,
      params = SKA1_build.build_parameters,
      snp_cutoff = snp_cutoff,
      identity_cutoff = identity_cutoff
  }

  if (generate_tree) {
    call SKA_align_masked {
      input:
        skf_files = mask_genomes.masked_skf,
        strain = strain_name,
        kmer_freq = min_kmer_freq,
        params = SKA1_build.build_parameters
    }
    call build_tree_masked {
      input:
        aligned_skf = SKA_align_masked.aligned_skf,
        strain = strain_name,
        params = SKA1_build.build_parameters
    }
  }

  if (generate_vcf) {
    call SKA1_annotate_masked {
      input:
        names = samplename,
        skf_files = mask_genomes.masked_skf,
        ref = ref_genome,
        strain = strain_name,
        params = SKA1_build.build_parameters
    }
  }

  output {
    File skf_summary_masked = SKA1_distance_masked.summaries
    File? ska_vcfs_masked = SKA1_annotate_masked.vcfs
    File? ska_tree_masked = build_tree_masked.treefile
    File ska_distance_masked = SKA1_distance_masked.distance_matrix
    File ska_clusters_masked = SKA1_distance_masked.clusters
  }
}

#### Tasks ####

task SKA1_build {
  input {
    File fq1
    File fq2
    String name
    Float? minor_freq
    Int? kmers
    Int? file_cutoff
    Int? total_cutoff
  }

  Float MAF = select_first([minor_freq, 0.2])
  Int kmers_actual = select_first([kmers, 15])
  Int file_cov = select_first([file_cutoff, 4])
  Int total_cov = select_first([total_cutoff, 2])
  String skf_summary = "~{name}_k15_summary.txt"
  String input_parameters = "~{kmers_actual}_~{MAF}_~{total_cov}_~{file_cov}"

  command <<<
    ska fastq -m ~{MAF} -k ~{kmers_actual} -c ~{total_cov} -C ~{file_cov} -o ~{name} ~{fq1} ~{fq2}
    ska summary ~{name}.skf > ~{skf_summary}
  >>>

  output {
    File skf_file = glob("*.skf")[0]
    File skf_summary = skf_summary
    String build_parameters = input_parameters
  }

  runtime {
    docker: "staphb/ska:latest"
    memory: "8 GB"
    preemptible: 0
    maxRetries: 1
    cpu: 4
  }
}

task extract_all_plasmids {
  input {
    File plasmid_tarball
    Array[String] sample_names
  }

  command <<<
    # Extract the tar.gz file
    tar -xzf ~{plasmid_tarball}

    # Find the combined plasmid fasta file
    PLASMID_FILE=$(find . -name "plasmids.fasta" | head -1)

    if [ -z "$PLASMID_FILE" ]; then
      echo "ERROR: plasmids.fasta not found in tarball"
      # Create empty files for all samples
      while IFS= read -r sample; do
        touch "${sample}_plasmids.fasta"
      done < ~{write_lines(sample_names)}
    else
      # Split plasmids by sample ID (extracted from header suffix plasmid_SAMPLEID)
      python3 <<PYTHON
import re
import sys

# Read sample names
with open('~{write_lines(sample_names)}', 'r') as f:
    samples = [line.strip() for line in f]

# Create empty fasta files for all samples first
sample_files = {}
for sample in samples:
    fname = f"{sample}_plasmids.fasta"
    sample_files[sample] = open(fname, 'w')

# Parse the combined plasmid fasta and assign to samples
current_sample = None
current_header = None
current_seq = []

with open("$PLASMID_FILE", 'r') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            # Write previous sequence if exists
            if current_sample and current_header:
                sample_files[current_sample].write(current_header + '\\n')
                sample_files[current_sample].write(''.join(current_seq) + '\\n')

            # Extract sample ID from header (looks for plasmid_SAMPLEID at end)
            match = re.search(r'plasmid_([A-Za-z0-9_-]+)', line)
            if match:
                sample_id = match.group(1)
                if sample_id in sample_files:
                    current_sample = sample_id
                    current_header = line
                    current_seq = []
                else:
                    current_sample = None
            else:
                current_sample = None
        else:
            if current_sample:
                current_seq.append(line)

    # Write last sequence
    if current_sample and current_header:
        sample_files[current_sample].write(current_header + '\\n')
        sample_files[current_sample].write(''.join(current_seq) + '\\n')

# Close all files
for f in sample_files.values():
    f.close()
PYTHON
    fi

    # Verify all sample files exist (create empty if missing)
    while IFS= read -r sample; do
      if [ ! -f "${sample}_plasmids.fasta" ]; then
        touch "${sample}_plasmids.fasta"
      fi
    done < ~{write_lines(sample_names)}
  >>>

  output {
    Array[File] sample_plasmid_fastas = glob("*_plasmids.fasta")
  }

  runtime {
    docker: "python:3.9-slim"
    memory: "4 GB"
    cpu: 1
    preemptible: 0
  }
}

task plasmid_to_skf {
  input {
    File plasmid_fasta
    String sample_id
  }

  String skf_name = "~{sample_id}_plasmids"

  command <<<
    # Only convert if file is not empty
    if [ -s ~{plasmid_fasta} ]; then
      ska fasta -o ~{skf_name} ~{plasmid_fasta}
    else
      # Create empty skf placeholder
      touch ~{skf_name}.skf
    fi
  >>>

  output {
    File plasmid_skf = "~{skf_name}.skf"
  }

  runtime {
    docker: "staphb/ska:latest"
    memory: "2 GB"
    cpu: 1
    preemptible: 0
  }
}

task mask_genomes {
  input {
    File unmasked_skf
    File plasmid_skf
    String sample_id
  }

  String masked_name = "~{sample_id}_masked"

  command <<<
    cp ~{unmasked_skf} ~{masked_name}.skf

    # Only weed if plasmid file is not empty
    if [ -s ~{plasmid_skf} ]; then
      ska weed -i ~{plasmid_skf} ~{masked_name}.skf
    else
      # No plasmids found, rename to match expected output
      mv ~{masked_name}.skf ~{masked_name}.weeded.skf
    fi
  >>>

  output {
    File masked_skf = "~{masked_name}.weeded.skf"
  }

  runtime {
    docker: "staphb/ska:latest"
    memory: "4 GB"
    cpu: 1
    preemptible: 0
  }
}

task SKA1_distance_masked {
  input {
    String strain
    Array[File] skf_summary
    Array[File] skf_files
    Array[String] params
    Float? identity_cutoff
    Int? snp_cutoff
  }

  Float identity_cutoff_actual = select_first([identity_cutoff, 0.9])
  Int snp_cutoff_actual = select_first([snp_cutoff, 20])
  String skf_filelist = "all_masked_skf_files.txt"
  String user_params = params[0]
  String skf_distances_named = "~{strain}_masked_~{user_params}"

  command <<<
    skf_array=(~{sep=" " skf_files})
    for i in ${skf_array[@]}; do echo $i >> ~{skf_filelist}; done
    ska distance -f ~{skf_filelist} -i ~{identity_cutoff_actual} -s ~{snp_cutoff_actual} -o ~{skf_distances_named}

    touch all_summaries_files.txt
    summary_array=(~{sep=" " skf_summary})
    for i in ${summary_array[@]}; do cat $i >> all_summaries_files.txt; done
    sed '1!{/^Sample/d;}' all_summaries_files.txt > all_summaries_files_clean.txt
    mv all_summaries_files_clean.txt ~{skf_distances_named}_summaries.txt
  >>>

  output {
    File distance_matrix = "~{skf_distances_named}.distances.tsv"
    File clusters = "~{skf_distances_named}.clusters.tsv"
    File summaries = "~{skf_distances_named}_summaries.txt"
  }

  runtime {
    docker: "staphb/ska:latest"
    memory: "250 GB"
    disks: "local-disk 200 HDD"
  }
}

task SKA_align_masked {
  input {
    Array[File] skf_files
    Array[String] params
    String strain
    Float? kmer_freq
  }

  String user_params = params[0]
  String skf_distances_named = "~{strain}_masked_~{user_params}"
  Float kfreq = select_first([kmer_freq, 0.9])

  command <<<
    skf_array=(~{sep=" " skf_files})
    ska merge -o ~{skf_distances_named}_merged ${skf_array[@]}
    ska align -v -p ~{kfreq} -o ~{skf_distances_named} ~{skf_distances_named}_merged.skf
  >>>

  output {
    File aligned_skf = "~{skf_distances_named}_variants.aln"
  }

  runtime {
    docker: "staphb/ska:latest"
    memory: "50 GB"
    disks: "local-disk 200 HDD"
  }
}

task build_tree_masked {
  input {
    Array[String] params
    String strain
    File aligned_skf
  }

  String user_params = params[0]
  String skf_distances_named = "~{strain}_masked_~{user_params}"

  command <<<
    VeryFastTree -nt -gamma -gtr -threads 4 ~{aligned_skf} > ~{skf_distances_named}.txt
  >>>

  output {
    File treefile = "~{skf_distances_named}.txt"
  }

  runtime {
    docker: "vkhadka/veryfasttree:v4.0.5"
    memory: "5 GB"
    cpu: 4
    disks: "local-disk 200 HDD"
  }
}

task SKA1_annotate_masked {
  input {
    Array[String] names
    Array[File] skf_files
    File? ref
    Array[String] params
    String strain
  }

  String user_params = params[0]
  String skf_distances_named = "~{strain}_masked_~{user_params}"

  command <<<
    skf_array=(~{sep=" " skf_files})
    names_array=(~{sep=" " names})

    mkdir vcf_files

    for index in ${!skf_array[*]}; do
      ska annotate -r ~{ref} -o ${names_array[$index]} ${skf_array[$index]}
      head -7 ${names_array[$index]}.vcf > ${names_array[$index]}_filt.vcf && grep 'NS5' ${names_array[$index]}.vcf >> ${names_array[$index]}_filt.vcf
      sed -i '' '8d' ${names_array[$index]}_filt.vcf
      mv ${names_array[$index]}_filt.vcf vcf_files/
    done

    tar -czf ~{skf_distances_named}_vcf.tar.gz vcf_files
  >>>

  output {
    File vcfs = "~{skf_distances_named}_vcf.tar.gz"
  }

  runtime {
    docker: "staphb/ska:latest"
    memory: "100 GB"
    disks: "local-disk 200 HDD"
  }
}
