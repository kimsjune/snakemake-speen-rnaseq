GTF = "../index/genes/gencode.vM22.annotation.gtf"
mm10Genome = "../index/genomes/mm10.fa"
GENOTYPE = ["WT_B","RLR_B","SPL"]
TREATMENT = ["D","G"]
REPLICATES = ["1","2","3"]

ALL_SAMPLES = expand("{geno}_{treatment}{rep}", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES)

rule all:
    input:
        fastqc = expand("fastqc/{geno}_{treatment}{rep}.html", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES),
        star =  expand("htseqcount/{geno}_{treatment}{rep}.csv", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES),
        RepEnrich2= expand("RepEnrich2/table/{geno}_{treatment}{rep}_fraction_counts.txt", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES),



rule testfastqc:
    input: expand("fastqc/{geno}_{treatment}{rep}.html", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES)

rule fastqc:
    output:
        html = "fastqc/{sample}.html",
        zip  = "fastqc/{sample}.zip"
    input:  "../gz/gz_spleen/{sample}.fastq.gz"
    threads: 8
    shell:
       r"""fastqc -o . {input} -t {threads}
           mv {wildcards.sample}_fastqc.html {output.html}
           mv {wildcards.sample}_fastqc.zip  {output.zip}
        """
# limitGenomeGenerate needed to avoid error
# adjust runThreadN
rule STARindex:
    output: 
        folder= directory("../index/genes/gencode.vM22.annotation_STARindex"),
    input: 
        gtf = GTF,
        fa = mm10Genome,
    params:
        length = 75
    threads: 8
    shell:
        "STAR  --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.folder} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.length} --limitGenomeGenerateRAM 8000000000"


rule STARalignOnePass:
    output: 
        bam="{sample}_STARalignOnePass/Aligned.sortedByCoord.out.bam",
        log="{sample}_STARalignOnePass/Log.final.out"
        # needs to be this way because STAR adds extra string to the end of the file that can't be changed
    input: 
        genome_dir = "../index/genes/gencode.vM22.annotation_STARindex",
        fq = "../gz/gz_spleen/{sample}.fastq.gz"
    threads: 4
    shell:
        r"""
        STAR --runThreadN {threads} --genomeDir {input.genome_dir} --outSAMtype BAM SortedByCoordinate
        --readFilesCommand zcat --readFilesIn {input.fq} --outSAMstrandField
        intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {output.bam} >& {output.log} 
        """


#rule samtoolsViewSortIndex:
#   output: "STARalign/{sample}_sorted.bam"
#   input: "STARalign/{sample}.sam"
#   shell:
#       r"""
#       samtools view -b -@ {threads} {input} | samtools sort -@ {threads} - -o {output}
#        """

rule STARalignTwoPass:


rule htseqCount:
    output: "htseqcount/{sample}.csv"
    input:
        bam = "{sample}_STARalignOnePass/Aligned.sortedByCoord.out.bam",
        gtf = GTF
    threads: 8
    shell: "htseq-count -n {threads} -f bam -s no -q {input.bam} {input.gtf} > {output}"

rule testSTAR: 
    input: expand("htseqcount/{geno}_{treatment}{rep}.csv", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES)

rule multiqc:
    output:
        mqc_out = directory('multiqc_out'),
        mqc_in  = directory('multiqc_in'),
    input:
        star = expand("{geno}_{treatment}{rep}_STARalignOnePass/Log.final.out", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES),
        # although multiqc casn automatically find the log files, they should be explicitly stated for snakemake?
        fastqc = expand("fastqc/{geno}_{treatment}{rep}.zip", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES),
        bt2 = expand("bowtie2Align/{geno}_{treatment}{rep}.log", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES),
    shell:
      r"""mkdir {output.mqc_in}
          ln -snr -t {output.mqc_in} {input}
          multiqc {output.mqc_in} -o {output.mqc_out}
       """
rule bowtie2Index:
    output: "../index/genomes/bt2/mm10"
    log: "../index/genomes/bt2/mm10Log.txt"
    input:  "../index/genomes/mm10.fa"
    shell:
        r"""
        bowtie2-build {input} {output} >& {log}
        """


rule bowtie2Align:
    output: "bowtie2Align/{sample}.sam"
    log: "bowtie2Align/{sample}.log"
    input:
        index = "../index/genomes/bt2/mm10",
        fq = "../gz/gz_spleen/{sample}.fastq.gz"
    threads: 36
    shell:
        r"""
        bowtie2 -q --very-fast-local -p {threads} -x {input.index}
        -U {input.fq} -S {output} >& {log}
        """
rule samtoolsBT2:
    output: protected("bowtie2Align/{sample}.bam")
    # doesn't need to be sorted
    input: "bowtie2Align/{sample}.sam"
    threads: 4
    shell:
        r"""
        samtools view -b -@ {threads} {input} -o {output}
        """

rule RepEnrich2_setup:
    conda: "environment.yaml"
    output: directory("../index/repeats/RepEnrichSetup")
    input: 
        bed="../index/repeats/mm10_rmsk_filter_fix_rmblacklist.bed",
        fa="../index/genomes/mm10.fa"
    shell:
        r"""
        python2 scripts/RepEnrich2_setup.py {input.bed} {input.fa} {output} --is_bed TRUE
        """


rule RepEnrich2_subset:
    conda: "environment.yaml"
    output: 
        bam="RepEnrich2/{sample}_unique.bam",
        fq="RepEnrich2/{sample}_multimap.fastq"
    input: "bowtie2Align/{sample}.bam"
    params:
        MAPQ = 30,
        prefix = "{sample}"
    shell:
        r"""
        python2 scripts/RepEnrich2_subset.py {input} {params.MAPQ} {params.prefix}
        """
rule testRepEnrich2:
    input: expand("bowtie2Align/{geno}_{treatment}{rep}_unique.bam", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES)

rule samtoolsSortUnique:
    output: "RepEnrich2/{sample}_sorted.bam"
    input: "RepEnrich2/{sample}.bam"
    threads: 8
    shell:
        r"""
        samtools sort -@ {threads} {input} -o {output} 
        """
rule samtoolsIndex:
    output: "RepEnrich2/{sample}_unique_sorted.bai"
    input: "RepEnrich2/{sample}_unique_sorted.bam"
    shell:
        r"""
        samtools index -b {input}
        """
rule testIndex:
    input: expand("RepEnrich2/{geno}_{treatment}{rep}_unique_sorted.bai", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES)

rule RepEnrich2:
    conda: "environment.yaml"
    output: "RepEnrich2/table/{sample}_fraction_counts.txt"
    input:
       bed="../index/repeats/mm10_rmsk_filter_fix_rmblacklist.bed",
       folder="../index/repeats/RepEnrichSetup",
       bam = "RepEnrich2/{sample}_unique_sorted.bam",
       multifq = "RepEnrich2/{sample}_multimap.fastq"
    threads: 32
    params:
        prefix="RepEnrich2/table"
    shell:
        r"""
        python2 scripts/RepEnrich2.py {input.bed} {params.prefix} {wildcards.sample} {input.folder} {input.multifq} {input.bam} --cpus {threads} --is_bed TRUE
        """

rule testBT2:
    input: expand("bowtie2Align/{geno}_{treatment}{rep}.log", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES)

rule testRE2:
    input: expand("RepEnrich2/table/{geno}_{treatment}{rep}_fraction_counts.txt", geno=GENOTYPE, treatment=TREATMENT, rep=REPLICATES)