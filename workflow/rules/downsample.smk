


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


# calling
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
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/make_ds_gvcfs/{sample}/{sg_or_chrom}.bmk"
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


rule concat_ds_gvcf_sections:
    input: 
        expand("results/bqsr-round-{{bqsr_round}}/downsample-{{cov}}X/gvcf_sections/{{sample}}/{sgc}.g.vcf.gz", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        gvcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/gvcf/{sample}.g.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/gvcf/{sample}.g.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/concat_ds_gvcf_sections/{sample}.txt"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/concat_ds_gvcf_sections/{sample}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools concat {params.opts} -O z {input} > {output.gvcf} 2>{log}; "
        " bcftools index -t {output.gvcf} "


rule genomics_db_import_chromosomes_ds:
    input:
        gvcfs=lambda wc: expand("results/bqsr-round-{{bqsr_round}}/downsample-{{cov}}X/gvcf_sections/{sample}/{{chromo}}.g.vcf.gz", sample=sample_list),
        gvcf_idxs=lambda wc: expand("results/bqsr-round-{{bqsr_round}}/downsample-{{cov}}X/gvcf_sections/{sample}/{{chromo}}.g.vcf.gz.tbi", sample=sample_list),
    output:
        db=directory("results/bqsr-round-{bqsr_round}/downsample-{cov}X/genomics_db/{chromo}")
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/genomicsdbimport/{chromo}.log"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/genomics_db_import_chromosomes_ds/{chromo}.bmk"
    params:
        my_opts=chromo_import_gdb_opts,
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $(echo {input.gvcfs} | awk '{{for(i=1;i<=NF;i++) printf(\" -V %s \", $i)}}') "
        " {params.my_opts} {output.db} > {log} 2>&1  "
        


rule genomics_db_import_scaffold_groups_ds:
    input:
        gvcfs=lambda wc: expand("results/bqsr-round-{{bqsr_round}}/downsample-{{cov}}X/gvcf_sections/{sample}/{{scaff_group}}.g.vcf.gz", sample=sample_list),
        gvcf_idxs=lambda wc: expand("results/bqsr-round-{{bqsr_round}}/downsample-{{cov}}X/gvcf_sections/{sample}/{{scaff_group}}.g.vcf.gz.tbi", sample=sample_list),
        scaff_groups = config["scaffold_groups"],
    output:
        interval_list="results/bqsr-round-{bqsr_round}/downsample-{cov}X/gdb_intervals/{scaff_group}.list",
        db=directory("results/bqsr-round-{bqsr_round}/downsample-{cov}X/genomics_db/{scaff_group}")
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/genomicsdbimport/{scaff_group}.log"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/genomics_db_import_scaffold_groups_ds/{scaff_group}.bmk"
    params:
        my_opts=scaff_group_import_gdb_opts,
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output.interval_list}; "
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $(echo {input.gvcfs} | awk '{{for(i=1;i<=NF;i++) printf(\" -V %s \", $i)}}') "
        " {params.my_opts} {output.db} >{log} 2>&1; "


  
rule genomics_db2vcf_scattered_ds:
    input:
        genome="resources/genome.fasta",
        scatters="results/bqsr-round-{bqsr_round}/scatter_interval_lists/{sg_or_chrom}/{scatter}.list",
        db="results/bqsr-round-{bqsr_round}/downsample-{cov}X/genomics_db/{sg_or_chrom}",
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz",
        tbi="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/genotypegvcfs/{sg_or_chrom}/{scatter}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/genomics_db2vcf_ds/{sg_or_chrom}/{scatter}.bmk",
    params:
        gendb="results/bqsr-round-{bqsr_round}/downsample-{cov}X/genomics_db/{sg_or_chrom}",
        java_opts="-Xmx8g",  # I might need to consider a temp directory, too in which case, put it in the config.yaml
        pextra=" --genomicsdb-shared-posixfs-optimizations --only-output-calls-starting-in-intervals "
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "1-00:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenotypeGVCFs "
        " {params.pextra} "
        " -L {input.scatters} "
        " -R {input.genome} "
        " -V gendb://{params.gendb} "
        " -O {output.vcf} > {log} 2> {log} "


rule gather_scattered_ds_vcfs:
    input:
        vcf=lambda wc: get_scattered_ds_vcfs(wc, ""),
        tbi=lambda wc: get_scattered_ds_vcfs(wc, ".tbi"),
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sections/{sg_or_chrom}.vcf.gz",
        tbi="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sections/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gather_scattered_ds_vcfs/{sg_or_chrom}.txt"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/gather_scattered_ds_vcfs/{sg_or_chrom}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Oz {input.vcf} > {output.vcf}; "
        " bcftools index -t {output.vcf})  2>{log}; "



rule mark_dp0_as_missing_ds:
    input:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sections/{sg_or_chrom}.vcf.gz"
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/mark_dp0_as_missing_ds/{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}Xbenchmarks/mark_dp0_as_missing_ds/{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools +setGT {input.vcf} -- -t q -n . -i 'FMT/DP=0 | (FMT/PL[:0]=0 & FMT/PL[:1]=0 & FMT/PL[:2]=0)' | "
        " bcftools +fill-tags - -- -t 'NMISS=N_MISSING' | "
        " bcftools view -Oz - > {output.vcf}; "
        " bcftools index -t {output.vcf}) 2> {log} "




rule ds_bcf_concat:
    input:
        expand("results/bqsr-round-{{bqsr_round}}/downsample-{{cov}}X/hard_filtering/both-filtered-{sgc}.bcf", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        bcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/bcf/all.bcf",
        tbi="results/bqsr-round-{bqsr_round}/downsample-{cov}X/bcf/all.bcf.csi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/ds_bcf_concat/bcf_concat_log.txt"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/ds_bcf_concat/bcf_concat.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2> {log}; "





rule ds_bcf_concat_mafs:
    input:
        expand("results/bqsr-round-{{bqsr_round}}/downsample-{{cov}}X/hard_filtering/both-filtered-{sgc}-maf-{{maf}}.bcf", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        bcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/bcf/pass-maf-{maf}.bcf",
        tbi="results/bqsr-round-{bqsr_round}/downsample-{cov}X/bcf/pass-maf-{maf}.bcf.csi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/ds_bcf_concat_mafs/maf-{maf}.txt"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/ds_bcf_concat_mafs/maf-{maf}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "


# hard filtering
rule make_snp_vcf_ds:
    input:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/selectvariants/select-snps-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/make_snp_vcf_ds/selectvariants-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type SNP -O {output.vcf} > {log} 2>&1 "



rule make_indel_vcf_ds:
    input:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/bqsr-round-{bqsr_round}/downsample-{cov}X/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/selectvariants/select-indels-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/make_indel_vcf_ds/selectvariants-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type INDEL -O {output.vcf} > {log} 2>&1 "




rule hard_filter_snps_ds:
    input:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/variantfiltration/snps-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/hard_filter_snps_ds/variantfiltration-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        "  -filter 'QD < 2.0' --filter-name 'QD2' "
        "  -filter 'QUAL < 30.0' --filter-name 'QUAL30' "
        "  -filter 'SOR > 3.0' --filter-name 'SOR3' "
        "  -filter 'FS > 60.0' --filter-name 'FS60' "
        "  -filter 'MQ < 40.0' --filter-name 'MQ40' "
        "  -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' "
        "  -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' "
        " -O {output.vcf} > {log} 2>&1 "




rule hard_filter_indels_ds:
    input:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/gatk/variantfiltration/indels-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/hard_filter_indels_ds/variantfiltration-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        "  -filter 'QD < 2.0' --filter-name 'QD2' "
        "  -filter 'QUAL < 30.0' --filter-name 'QUAL30' "
        "  -filter 'FS > 200.0' --filter-name 'FS200' "
        "  -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' "
        " -O {output.vcf} > {log} 2>&1 "




rule bung_filtered_ds_vcfs_back_together:
    input:
        snp="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz",
        indel="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        snp_idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz.tbi",
        indel_idx="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/both-filtered-{sg_or_chrom}.bcf",
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/bung_filtered_ds_vcfs_back_together/bung-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/bung_filtered_ds_vcfs_back_together/bcftools-{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools concat -a {input.snp} {input.indel} | "
        " bcftools view -Ob > {output.vcf}; ) 2> {log} "


rule maf_filter_ds:
    input:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/both-filtered-{sg_or_chrom}.bcf"
    output:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/hard_filtering/both-filtered-{sg_or_chrom}-maf-{maf}.bcf"
    log:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/logs/maf_filter_ds/{sg_or_chrom}-maf-{maf}.log",
    params:
        maf="{maf}"
    benchmark:
        "results/bqsr-round-{bqsr_round}/downsample-{cov}X/benchmarks/maf_filter/{sg_or_chrom}-maf-{maf}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ob -i 'FILTER=\"PASS\" & MAF > {params.maf} ' "
        " {input} > {output} 2>{log} "
