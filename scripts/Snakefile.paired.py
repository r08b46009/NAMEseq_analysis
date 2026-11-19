import sys
from os.path import join
INDEX = config["index"]
INPUT_DIR = config["input_path"]
OUTPUT_DIR = config["output_path"].rstrip("/") + "/" # or to be replaced by your dir
NUM_THREADS = int(config["num_threads"])
SAMPLES = config["basename_list"]
ref_genome = config["fasta"]
DATE = config["date"]

Read = ["R1", "R2"]
Strand = ["fwd", "rev"]

# ruleorder: dedupe>trimming>bowtie2_align>mapped_fastq_filtering>sam_to_original_pe>samtools_calmd>sam_tag_filtering>samtools_sort_index>igvtools
rule all:
    input:
        # expand(OUTPUT_DIR + "deduped_fastq/{sample}_{read}_dedupe.fastq.gz", sample=SAMPLES, read=Read),
        # expand(OUTPUT_DIR + "trimmed_fastq/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=Read),
        # expand(OUTPUT_DIR + "trimmed_fastq/{sample}_{read}_trimmed2.fastq", sample=SAMPLES, read=Read),
        # expand(OUTPUT_DIR + "converted_sam/{sample}_trimmed2.sam", sample=SAMPLES),
        # expand(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}.sam", sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}.sorted.bam", sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}_baq.sam", sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_{strand}_baq_filtered.sam",date=DATE, sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",date=DATE, sample=SAMPLES, strand=Strand),
        expand(OUTPUT_DIR + "preprocessed_data_{date}/{sample}_adenine.csv",date = DATE, sample=SAMPLES),
        # expand(OUTPUT_DIR + "preprocessed_data_{date}/{sample}_cytosine.csv",date=DATE, sample=SAMPLES),

####


rule dedupe:
    input:
        r1 = INPUT_DIR + "/{sample}_R1.fastq.gz",
        r2 = INPUT_DIR + "/{sample}_R2.fastq.gz",
    output:
        r1 = OUTPUT_DIR + "deduped_fastq_{date}/{sample}_R1_dedupe.fastq.gz",
        r2 = OUTPUT_DIR + "deduped_fastq_{date}/{sample}_R2_dedupe.fastq.gz",
    log:
        OUTPUT_DIR + "logs/dedupe_{date}/{sample}.log"
    shell:
        """
        clumpify.sh -Xmx20G in1={input.r1} in2={input.r2} out={output.r1} out2={output.r2} dedupe 1>{log} 2>&1
        """
        
rule trimming:
    input:
        r1 = OUTPUT_DIR + "deduped_fastq_{date}/{sample}_R1_dedupe.fastq.gz",
        r2 = OUTPUT_DIR + "deduped_fastq_{date}/{sample}_R2_dedupe.fastq.gz",
    output:
        r1 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R1_trimmed.fastq.gz",
        r2 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R2_trimmed.fastq.gz",
    log:
        OUTPUT_DIR + "logs/cutadpt_{date}/{sample}.log"
    params:
        TRIM_OPTS = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 10"
    threads: NUM_THREADS
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} -u 9 -U 9 1>{log}
        """

rule trimming2:
    input:
        r1 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R1_trimmed.fastq.gz",
        r2 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R2_trimmed.fastq.gz",
    output:
        r1 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R1_trimmed2.fastq",
        r2 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R2_trimmed2.fastq",
    log:
        OUTPUT_DIR + "logs/cutadpt2_{date}/{sample}.log"
    params:
        TRIM_OPTS = "-m 10 -u -9 -U -9"
    threads: NUM_THREADS
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} 1>{log}
        """

rule bowtie2_align:
    input:
        r1 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R1_trimmed2.fastq",
        r2 = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_R2_trimmed2.fastq",
    output:
             OUTPUT_DIR + "converted_sam_{date}/{sample}_trimmed2.sam"
    params:
        index = INDEX,
    log:
        OUTPUT_DIR + "logs/bowtie_{date}/{sample}.log"
    threads: NUM_THREADS
    shell:
        "bowtie2 -x {params} -p {threads} -1 {input.r1} -2 {input.r2} -S {output} > {log} 2>&1 --seed 777777777 --verbose -k 40 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"


rule mapped_fastq_filtering:# (correct)
    input:
        sam =  OUTPUT_DIR + "converted_sam_{date}/{sample}_trimmed2.sam",
        fastq = OUTPUT_DIR + "trimmed_fastq_{date}/{sample}_{read}_trimmed2.fastq"
    output:
        fastq = OUTPUT_DIR + "mapped_fastq_{date}/{sample}_{read}_mapped.fastq",
        readid = OUTPUT_DIR + "mapped_fastq_{date}/{sample}_{read}_mapped.readid"
    shell:
        """
        samtools view -F 4 {input.sam} | awk '!seen[$1]++ {{print $1}}' > {output.readid}
        seqkit grep -f {output.readid} {input.fastq} > {output.fastq}
        """
rule sam_to_original_pe:
    input:
        sam =  OUTPUT_DIR + "converted_sam_{date}/{sample}_trimmed2.sam",
        r1 = OUTPUT_DIR + "mapped_fastq_{date}/{sample}_R1_mapped.fastq",
        r2 = OUTPUT_DIR + "mapped_fastq_{date}/{sample}_R2_mapped.fastq",
    output:
        OUTPUT_DIR + "original_sam_{date}/{sample}_original_reads_fwd.sam",
        OUTPUT_DIR + "original_sam_{date}/{sample}_original_reads_rev.sam"
    shell:
        """
        echo "Executing: python3 scripts/sam_to_original_pe.py {input.sam} {input.r1} {input.r2} {output[0]} {output[1]}"
        python3 scripts/sam_to_original_pe.py {input.sam} {input.r1} {input.r2} {output[0]} {output[1]} 
        """

rule samtools_calmd:
    input:
        OUTPUT_DIR + "original_sam_{date}/{sample}_original_reads_{strand}.sam",
        ref_genome = ref_genome,
    output:
        OUTPUT_DIR + "original_sam_{date}/{sample}_original_reads_{strand}.sorted.bam",
        OUTPUT_DIR + "original_sam_{date}/{sample}_original_reads_{strand}_baq.sam"
    threads: NUM_THREADS
    shell:
        """
        samtools sort {input[0]} > {output[0]}
        samtools calmd -@ {threads} {output[0]} {input.ref_genome} > {output[1]} 2>/dev/null
        """
        
        
        
####
rule sam_tag_filtering:
    input:
        fwd = OUTPUT_DIR + "original_sam_{date}/{sample}_original_reads_fwd_baq.sam",
        rev = OUTPUT_DIR + "original_sam_{date}/{sample}_original_reads_rev_baq.sam"
    output:
        fwd = OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_fwd_baq_filtered.sam",
        rev = OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_rev_baq_filtered.sam"

    shell:
        """
        python3 scripts/sam_tag_filtering.py {input.fwd} {output.fwd}
        python3 scripts/sam_tag_filtering.py {input.rev} {output.rev}
        """

rule samtools_sort_index:
    input:
        OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_{strand}_baq_filtered.sam"
    output:
        OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",
        OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_{strand}_baq_filtered.sorted.bam.bai"
    threads: NUM_THREADS
    shell:
        """
        samtools sort -@ {threads} {input} > {output[0]}
        samtools index {output[0]}
        """

rule igvtools:
    input:
        bam = OUTPUT_DIR + "filtered_sam_{date}/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",
        ref_genome = ref_genome,
    output:
        OUTPUT_DIR + "readcount_{date}/{sample}_original_reads_{strand}_baq_filtered.wig"
    shell:
        """
        igvtools count --bases -w 1 -e 0 {input.bam} {output} {input.ref_genome} --includeDuplicates
        """
rule data_proprocess:
    input:
        fwd = OUTPUT_DIR + "readcount_{date}/{sample}_original_reads_fwd_baq_filtered.wig",
        rev = OUTPUT_DIR + "readcount_{date}/{sample}_original_reads_rev_baq_filtered.wig",
        ref_genome = ref_genome,
    output:
        OUTPUT_DIR + "preprocessed_data_{date}/{sample}_adenine.csv"
    shell:
        """
        python3 scripts/readcount_to_csv.py {input.fwd} {input.rev} {input.ref_genome} {output}
        """
        
# # rule data_proprocess_cytosine:
# #     input:
# #         fwd = OUTPUT_DIR + "readcount_{date}/{sample}_original_reads_fwd_baq_filtered.wig",
# #         rev = OUTPUT_DIR + "readcount_{date}/{sample}_original_reads_rev_baq_filtered.wig",
# #         ref_genome = ref_genome,
# #     output:
# #         OUTPUT_DIR + "preprocessed_data_{date}/{sample}_cytosine.csv"
# #     shell:
# #         """
# #         python3 scripts/readcount_to_csv-cytosine.py {input.fwd} {input.rev} {input.ref_genome} {output}
# #         """
