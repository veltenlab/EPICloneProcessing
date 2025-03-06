configfile: "config.yaml"

names = config['name']
samples = config['sample']
dirs = config['path']
batch = config['batch']

## Here, you will have to provide the path to the installation of the on-premise Mission Bio pipeline, which will create a conda environment for you. Please specify it here.
TAPESTRI=""

rule all:
    input:
        expand("barcode/{name}_read1.with.barcodes.fastq.gz", name=names), 
        expand("barcode/{name}_read2.with.barcodes.fastq.gz", name=names), 
        expand("bam/{name}_aligned.bam", name=names), 
        expand("bam/{name}_aligned_fixed.bam", name=names),
        expand("distribution/{name}.cellfinder.barcode.distribution.txt", name=names),
        expand("tsv/{name}.barcode.cell.distribution.tsv", name=names),
        expand("tsv/{name}_doublet_scores_DoubletDetection.csv", name=names),
        expand("plots/{name}_reads_per_cell.png", name=names),
        expand("larry/{name}_LARRY_sorted_and_filtered_barcodes.fastq.gz", name=names),
        "larry/clones_{batch}"

rule barcodecode:
    input:
        first = lambda wildcards: f"{config['path']}/{config['samples'][wildcards.name]}_R1_001.fastq.gz", 
        second = lambda wildcards: f"{config['path']}/{config['samples'][wildcards.name]}_R2_001.fastq.gz"
    output:
        first = "barcode/{name}_read1.with.barcodes.fastq.gz",
        second = "barcode/{name}_read2.with.barcodes.fastq.gz",
        barcode = "barcode/{name}_raw.barcodes.txt",
        metrics = "barcode/{name}_01_barcode.metrics.json",
        sam_header = "tmp/{name}.sam.header"
    params:
        threads = 4
    log:
        "barcode/{name}_barcodecode.log"
    shell:
        "conda activate $TAPESTRI ; "
        "$TAPESTRI/bin/tapestri barcode extract --input-fwd {input.first} --input-rev {input.second} --output-fwd {output.first} --output-rev {output.second} --output-barcodes {output.barcode} --output-metrics {output.metrics} --chemistry dna.V2 --no-cut-adapters --decode --cores {params.threads} --barcode-suffix -1 ; "
        r"sed -n 's/^\([^\t]*\).*/@RG\tID:\1\tSM:\1/p' {output.barcode} > {output.sam_header} ; "
        "2> {log}"

rule alignment:
    input:
        fastq1 = "barcode/{name}_read1.with.barcodes.fastq.gz",
        fastq2 = "barcode/{name}_read2.with.barcodes.fastq.gz",
        sam_header = "tmp/{name}.sam.header"
    output:
        "bam/{name}_aligned.bam"
    params:
        tmp_ref = "tmp/panel_reference_{name}.fasta",
        tmp_sai = "tmp/tmp_{name}.sai",
        tmp_sam = "tmp/tmp_{name}.sam",
        tmp_bam = "tmp/tmp_{name}.bam",
        threads = 6
    log:
        "bam/{name}_bwa.log"
    conda:
        "genomics.yaml"
    shell:
        "Rscript scripts/generate_costum_reference_genome.R --amplicon {config[panel_bed]} --output {params.tmp_ref} --assembly {config[assembly]} ; "
        "conda activate $TAPESTRI ; "
        "bwa index {params.tmp_ref} ; "
        "bwa mem -C -M -t {params.threads} -H {input.sam_header} {params.tmp_ref} {input.fastq1} {input.fastq2} > {params.tmp_sam} ; "
        "samtools view -@ {params.threads} -S -b {params.tmp_sam} > {params.tmp_bam} ; "
        "samtools view -@ {params.threads} -bh -q 30 -F 4 -F 8 -F 0X0100 {params.tmp_bam} > {output} ; "
        "rm -rf {params.tmp_sam} ; "
        "rm -rf {params.tmp_bam} ; "
        "rm -rf {params.tmp_ref} ; "
        "2> {log}"

rule picard:
    input:
        "bam/{name}_aligned.bam",
    output:
        "bam/{name}_aligned_fixed.bam"
    log:
        "bam/{name}_picard.log"
    shell:
        "conda activate $TAPESTRI ; "
        "picard FixMateInformation -Xmx4g -Djava.io.tmpdir=tmp I={input} O={output} SORT_ORDER=coordinate ; "
        "2> {log}"

rule distribution:
    input:
        "bam/{name}_aligned_fixed.bam",
    output:
        distribution = "distribution/{name}.cellfinder.barcode.distribution.txt",
    params:
        cf = "distribution/{name}.mapped.target.count.txt",
        minus_5 = "tmp/{name}_panel_minus_5.bed",
        foo = "tmp/{name}_tmp_count.txt",
        tmp = "tmp/{name}_foo.txt",
        threads = 6
    log:
        "distribution/{name}_picard.log"
    shell:
        "conda activate $TAPESTRI ; "
        r"""awk '{{print $1 "\t" ($2 + 5) "\t" ($3 - 5) "\t" $4}}' {config[panel_bed]} > {params.minus_5} ; """
        r"""samtools view -@ {params.threads} {input} | grep -oP 'RG:Z:\K.+?\b' | sort | uniq -c | sort -n -r > {params.foo} ; """
        r"""awk '{{$1=$1}}1' {params.foo} > {params.cf} ; """
        "perl scripts/distribution_barcodes_amplicon_local_alignment.pl --amplicon-input {params.minus_5} --bam-input {input} --output {params.tmp} ; "
        "tr -d $'\r' < {params.tmp} > {output.distribution} ; "
        "2> {log}"
        
rule cellfinder:
    input:
        "distribution/{name}.cellfinder.barcode.distribution.txt",
    output:
        "tsv/{name}.barcode.cell.distribution.tsv"
    conda:
        "argparse.yml"
    log:
        "tsv/{name}_cellfinder.log"
    shell:
        "Rscript scripts/cellfinder_HSCs_selected.R -f {input} -a {config[ampli_file]} -c {config[cellfinder_cutoff]} -o {output} ; "
        "2> {log}"

rule doubletdetection:
    input:
        "tsv/{name}.barcode.cell.distribution.tsv",
    output:
        "tsv/{name}_doublet_scores_DoubletDetection.csv"
    log:
        "tsv/{name}_doublet_detection.log"
    shell:
        "conda activate $TAPESTRI ; "
        "python3 scripts/DoubletDetection.py --input {input} --output {output} ; "
        "2> {log}"

rule reads_per_cell:
    input:
        cells="distribution/{name}.cellfinder.barcode.distribution.txt",
        doublets="tsv/{name}_doublet_scores_DoubletDetection.csv"
    output:
        "plots/{name}_reads_per_cell.png"
    log:
        "plots/{name}_read_per_cell.log"
    params:
        folder = "plots/",
        n = "{name}"
    conda:
        "genomics.yaml"
    shell:
        "Rscript scripts/plot_reads_per_cell.R --file {input.cells} --doublet {input.doublets} --amplicon {config[ampli_file]} --output {params.folder} --name {params.n}"

rule extract_barcode:
    input: 
        fastq="barcode/{name}_read1.with.barcodes.fastq.gz",
        barcodes="tsv/{name}.barcode.cell.distribution.tsv"
    output:
        "larry/{name}_LARRY_sorted_and_filtered_barcodes.fastq.gz"
    params:
        folder = "larry/",
        name = "{name}",
        barcode_extension=lambda wildcards: f"{config['barcode_extension'][wildcards.name]}"
    log:
        "larry/{name}_larry.log"
    conda:
        "larry.yaml"
    shell:
        """
        python scripts/extract_LARRY_barcodes_v2_plasmid_batch.py --fastq {input.fastq} --barcode {input.barcodes} --output {params.folder} --name {params.name} --barcode_extension {params.barcode_extension};
        """

rule clonal_annotation:
    input:
        extracted="larry/{name}_LARRY_sorted_and_filtered_barcodes.fastq.gz",
        barcodes="tsv/{name}.barcode.cell.distribution.tsv"
    output:
        directory("larry/clones_{batch}")
    params:
        all_larry="all_extracted.fastq.gz",
        all_barcodes="all_barcodes.tsv"        
    log: 
        "larry/clones_{batch}/annotation.log"
    conda:
        "larry.yaml"
    shell:
        "zcat {input.extracted} | gzip -c > {params.all_extracted} ; "
        "cat {input.barcodes} | gzip -c > {params.all_barcodes} ; "
        "python scripts/LARRY_clonal_annotation.py --fastq {params.all_extracted} --barcode {params.all_barcodes} --output {output} --nreads {config[nreads_larry]} --nhamming {config[hamming_larry]} ; "
        "rm -rf tmp/"
