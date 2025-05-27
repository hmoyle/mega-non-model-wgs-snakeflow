


rule get_genome_length:
    input:
        "resources/genome.fasta.fai"
    output:
        "results/bqsr-round-{bqsr_round}/DS_control/genome_length.txt"
    shell:
        " awk 'NF>0 {{clen = $3}} END {{print clen}}' {input} > {output} "


rule get_ave_depths:
    input:
        gLen="results/bqsr-round-{bqsr_round}/DS_control/genome_length.txt",
        ss=expand("results/bqsr-round-{{bqsr_round}}/qc/samtools_stats/{sample}.txt", sample = sample_list)
    output:
        "results/bqsr-round-{bqsr_round}/DS_control/sample_info.tsv"
    shell:
        " ("
        " printf \"sample\\tave_depth\\n\";   "
        " for i in {input.ss}; do "
        " FN=$(basename $i);    "
        " FN=${{FN/.txt/}};       "
        " awk -F\"\t\" -v f=$FN -v NumBases=$(cat {input.gLen}) '  "
        "   BEGIN {{OFS=\"\t\";}} "
        "   $2==\"bases mapped (cigar):\" {{sub(/\\.stats/, \"\", f); print f,  $3/NumBases}} "
        " ' $i; done) > {output}  "




# in here the subsample seed is simply hardwired at 1.
# It didn't seem that we would want to do multiple reps of subsampling.
# If the file is already below the downsampling depth this will merely
# hard link the file to the new location.  I previously soft-linked relatively
# but ln -sr is not available on Mac.  So I just hard link it!
rule thin_bam:
    input:
        bam="results/bqsr-round-{bqsr_round}/overlap_clipped/{sample}.bam",
        bai="results/bqsr-round-{bqsr_round}/overlap_clipped/{sample}.bam.bai",
        dps="results/bqsr-round-{bqsr_round}/DS_control/sample_info.tsv"
    output:
        bam="results/bqsr-round-{bqsr_round}/downsample-{cov}X/overlap_clipped/{sample}.bam",
        bai="results/bqsr-round-{bqsr_round}/downsample-{cov}X/overlap_clipped/{sample}.bam.bai"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/thin_bams/{sample}.log"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/thin_bams/{sample}.bmk"
    conda:
        "../envs/samtools.yaml"
    shell:
        " ( "
        " OPT=$(awk 'NR>1 && $1==\"{wildcards.sample}\" {{ wc = \"{wildcards.cov}\"; if(wc == \"FD\") {{print \"NOSAMPLE\"; exit}} fract = wc / $NF; if(fract < 1) print fract; else print \"NOSAMPLE\"; }}' {input.dps});  "
        " if [ $OPT = \"NOSAMPLE\" ]; then "
        "     ln  {input.bam} {output.bam}; "
        "     ln  {input.bai} {output.bai}; " 
        " else "
        "     samtools view --subsample $OPT --subsample-seed 1  -b {input.bam} > {output.bam}; "
        "     samtools index {output.bam}; "
        " fi "
        " ) 2> {log} "

rule make_ds_gvcf_sections:
    input:
        bam="results/bqsr-round-{bqsr_round}/downsample-{cov}X/overlap_clipped/{sample}.bam",
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        fai="resources/genome.fasta.fai",
        interval_list="results/bqsr-round-{bqsr_round}/interval_lists/{sg_or_chrom}.list"
    output:
        gvcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz.tbi",
    conda:
        "../envs/gatk4.2.6.1.yaml"
    log:
        stderr="results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/haplotypecaller/{sample}/{sg_or_chrom}.stderr",
        stdout="results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/haplotypecaller/{sample}/{sg_or_chrom}.stdout",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/make_gvcfs/{sample}/{sg_or_chrom}.bmk"
    params:
        java_opts="-Xmx4g",
        conf_pars=config["params"]["gatk"]["HaplotypeCaller"]
    resources:
        time="1-00:00:00",
        mem_mb = 4600,
        cpus = 1
    threads: 1
    shell:
        "gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
        " -R {input.ref} "
        " -I {input.bam} "
        " -O {output.gvcf} "
        " -L {input.interval_list} "
        " --native-pair-hmm-threads {threads} "
        " {params.conf_pars} "
        " -ERC GVCF > {log.stdout} 2> {log.stderr} "

