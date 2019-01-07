configfile: "config.yaml"
configfile: "samples.yaml"


##### Target rules #####

rule all:
  input:
    "barcode_metrics.txt",
    expand("demuxed/{samplename}.bam", samplename=config["samples"])

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

rule illumina_basecalls_to_sam:
    input:
      basecalls_dir="basecallsdir",
      barcodes_dir="barcodes",
      library_params="library_params.txt"
    output:
      demuxed=expand("demuxed/{samplename}.bam", samplename=config["samples"])
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




  #Samplename Example   cllid_0381_121407_P_g01.bam
  #Samplename Structure cllid_####_MMDDYY_T_g##
