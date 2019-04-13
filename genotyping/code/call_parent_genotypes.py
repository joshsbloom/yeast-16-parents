#https://snakemake.bitbucket.io/snakemake-tutorial.html

samples = glob.glob("data/*fqL.gz")
samples = [os.path.basename(sample).split(".fqL.gz")[0] for sample in samples]

PICARD_PATH="/media/theboocock/data/james/variant_calling/software/picard.jar"
GATK_PATH="/media/theboocock/data//james/variant_calling/software/GenomeAnalysisTK.jar"

REFERENCE="/media/theboocock/data/james/variant_calling/data/reference/sacCer3.fasta"
rule all:
    input:
        "outputs/vcf/all.vcf"

rule bwa_map:
    input:
        read_one="data/{sample}.fqL.gz",
        read_two="data/{sample}.fqR.gz"
    output:
        "outputs/{sample}.bam",
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:lib1\\tPU:unit1"
    log:
        "logs/bwa_map/{sample}.log"
    shell:
        "(bwa mem -R '{params.rg}' -M {REFERENCE} {input.read_one} {input.read_two} | samtools view -Sb - > {output}) 2> {log}" 

rule picard_sort:
    input:
        "outputs/{sample}.bam"
    output:
        "outputs/{sample}.sorted.bam"
    params:
        order="SORT_ORDER=coordinate",
        tmpdir="-Djava.io.tmpdir=tmpdir"
    log:
        "logs/picard_sort/{sample}.log"
    shell:
        "java -Xmx2g {params.tmpdir} -jar {PICARD_PATH} SortSam INPUT={input} OUTPUT={output} {params.order} 2> {log}"

rule picard_mark_duplicates:
    input:
        "outputs/{sample}.sorted.bam"
    output:
        bam_out="outputs/md/{sample}.sorted.md.bam",
        metrics_out="outputs/md/{sample}.metrics.txt"
    log:
        "logs/picard_mark_duplicates/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir"
    shell:
        "java -Xmx2g {params.tmpdir} -jar {PICARD_PATH} MarkDuplicates INPUT={input} OUTPUT={output.bam_out} METRICS_FILE={output.metrics_out} 2> {log}"

rule picard_index:
    input:
        "outputs/md/{sample}.sorted.md.bam"
    output:
        "outputs/md/{sample}.sorted.md.bam.bai"
    log:
        "logs/picard_index/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir"
    shell:
        "java -Xmx2g {params.tmpdir} -jar {PICARD_PATH} BuildBamIndex INPUT={input} OUTPUT={output} 2> {log}"

rule haplotype_caller:
    input:
        bam="outputs/md/{sample}.sorted.md.bam",
        bam_index="outputs/md/{sample}.sorted.md.bam.bai"
    output:
        "outputs/gvcfs/{sample}.g.vcf"
    log:
        "logs/haplotype_caller/${sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir"
    shell:
        "java -Xmx8g {params.tmpdir} -jar {GATK_PATH} " 
        "-T HaplotypeCaller " 
        "-R {REFERENCE} "
        "-I {input.bam} "
        " -dcov 200 "
        "--emitRefConfidence GVCF "
        "--genotyping_mode DISCOVERY " 
        " --sample_ploidy 1 "
        "-o {output} " 

def _gatk_multi_arg(flag, files):
    cmd_string = flag + flag.join(files)
    return cmd_string 

rule haplotype_joint:
    input:
        expand("outputs/gvcfs/{sample}.g.vcf",sample=samples)
    output:
        "outputs/vcf/all.vcf"
    log:
        "logs/haplotype_caller/joint_calling.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir"
    run:
        bams = _gatk_multi_arg(" --variant ", input)
        shell("java -Xmx64g {params.tmpdir} -jar {GATK_PATH}"
    " -T GenotypeGVCFs "
    " -R {REFERENCE} "
        " -nt 8 "  
    "{bams}"
    " -o {output} ") 
