configfile: "config.yaml"
configfile: "samples.yaml"
samplename = config["samples"]

##### Target rules #####

rule all:
  input:
    "barcode_metrics.txt",
    expand("/umireads/{samplename}_consesnsus_unmapped.bam", samplename=config["samples"])

#First step is to extract the barcodes from the raw .bcl
rule extract_illumina_barcodes:
    input: basecalls_dir="basecallsdir", barcodes_file="barcode_file.txt"
    output: barcodes_dir="barcodes", metrics="barcode_metrics.txt"
    shell:
        "module load picard/2.8.0 "
        "mkdir -p {output.barcodes_dir}; "
        "java -Xmx4g -jar /nfs/sw/picard-tools/picard-tools-2.8.0/picard.jar ExtractIlluminaBarcodes BASECALLS_DIR={input.basecalls_dir} "
        "OUTPUT_DIR={output.barcodes_dir} "
        "LANE=1 "
        "BARCODE_FILE={input.barcodes_file} "
        "READ_STRUCTURE=100T8B9M8B100T "
        "METRICS_FILE={output.metrics} "

#This step does the demux based on the 2 sample indeces
rule illumina_basecalls_to_sam_demux:
    input:
      basecalls_dir="basecallsdir",
      barcodes_dir="barcodes",
      library_params="library_params.txt"
    output:
      demuxed=expand("/demuxed/{samplename}_unmapped.bam", samplename=config["samples"])
    shell:
        "module load picard/2.8.0"
        "java -Xmx8g -jar /nfs/sw/picard-tools/picard-tools-2.8.0/picard.jar IlluminaBasecallsToSam "
        "BASECALLS_DIR={input.basecalls_dir} "
        "BARCODES_DIR={input.barcodes_dir}  "
        "LANE=1 "
        "READ_STRUCTURE=100T8B9M8B100T "
        "RUN_BARCODE=CLL_pilot "
        "LIBRARY_PARAMS={input.library_params} "
        "TMP_DIR=/mnt/tmp "
        "MOLECULAR_INDEX_TAG=RX "
        "ADAPTERS_TO_CHECK=INDEXED "
        "READ_GROUP_ID=BN573-S1 "
        "NUM_PROCESSORS=8 "


#This command consists of three steps:
#1) Convert BAM to FASTQ, 2)Align reads using BWA-MEM, 3)Include UMI tags from unmapped BAM in the mapped BAM
rule unmapped_bam_to_mapped_bam_with_umi:
    input:
      unmapped_sample="/demuxed/{samplename}_unmapped.bam",
    output:
      mapped_sample="/demuxed/{samplename}_mapped.bam"
    shell:
      "module load picard/2.8.0"
      "java -Xmx4g -jar /nfs/sw/picard-tools/picard-tools-2.8.0/picard.jar SamToFastq "
      "I={input.unmapped_sample}"
      "F=/dev/stdout INTERLEAVE=true "
          "| bwa mem -p -t 8 hg38.fa /dev/stdin "
          "| java -Xmx4g -jar /nfs/sw/picard-tools/picard-tools-2.8.0/picard.jar MergeBamAlignment "
          "/demuxed/{samplename}_unmapped.bam ALIGNED=/dev/stdin "
          "O={output.mapped_sample} R=hg38.fa "
          "SORT_ORDER=coordinate MAX_GAPS=-1 "
          "ORIENTATIONS=FR"

#The reads are grouped into familes that share the same UMI
rule GroupReadsByUmi:
    input:
      mapped_sample="/demuxed/{samplename}_mapped.bam"
    output:
      grouped_sample="/umireads/{samplename}_grouped.bam"
    shell:
      "module load java/1.8"
      "wget https://github.com/fulcrumgenomics/fgbio/releases/download/0.7.0/fgbio-0.7.0.jar"
      "java -Xmx1g -jar fgbio-0.7.0.jar GroupReadsByUmi "
          "--{input.mapped_sample} --output={output.grouped_sample} "
          "--strategy=adjaceny --edits=1 --min-map-q=20 "
          "--assign-tag=MI "



#Consensus reads will be generate using fgbio's CallMolecularConsensusReads with default parameters
rule CallMolecularConsensusReads:
    input:
        grouped_sample="/umireads/{samplename}_grouped.bam"
    output:
        consensus_unmapped="/umireads/{samplename}_consesnsus_unmapped.bam"
    shell:
        "java -Xmx4g -jar fgbio-0.7.0.jar CallMolecularConsensusReads "
            "--input={input.grouped_sample} "
            "--output={output.consensus_unmapped} "
            "--min-reads=1 "
            "--rejects=/umireads/{samplename}_consesnsus_rejected.bam "
            "--min-input-base-quality=30 "
            "--read-group-id={samplename}"
