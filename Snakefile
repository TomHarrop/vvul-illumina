#!/usr/bin/env python3

import pathlib2

#############
# FUNCTIONS #
#############


def readset_wildard_resolver(wildcards):
    if wildcards.read_set == 'norm':
        return({'fq': 'output/030_norm/vvul.fq.gz'})
    elif wildcards.read_set == 'trim-decon':
        return({'fq': 'output/010_trim-decon/vvul.fq.gz'})
    else:
        raise ValueError('unknown read_set')


def write_config_file(fastq, threads, config_string, config_file):
    '''
    Accept fastq file, threads config string and output location and write
    config
    '''
    my_fastq = str(pathlib2.Path(fastq).resolve())
    my_conf = config_string.format(my_fastq, threads)
    with open(config_file, 'wt') as f:
        f.write(my_conf)
    return True


###########
# GLOBALS #
###########

r1_raw = ['data/C7N3UANXX-1782-01-04-01_L008_R1.fastq',
          'data/C7N3UANXX-1782-01-06-01_L008_R1.fastq']
r2_raw = ['data/C7N3UANXX-1782-01-04-01_L008_R2.fastq',
          'data/C7N3UANXX-1782-01-06-01_L008_R2.fastq']
bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'
meraculous_config_file = 'src/meraculous_config.txt'
meraculous_threads = 50

# containers
kraken_container = 'shub://TomHarrop/singularity-containers:kraken_2.0.7beta'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
mer_container = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.1'


########
# MAIN #
########

# read the meraculous config
with open(meraculous_config_file, 'rt') as f:
    meraculous_config_string = ''.join(f.readlines())


#########
# RULES #
#########

rule target:
    input:
        expand(('output/040_meraculous/{read_set}_k{k}_diplo{diplo}/'
                'meraculous_final_results/final.scaffolds.fa'),
               read_set=['norm', 'trim-decon'],
               k=['31', '71', '101'],
               diplo=['0', '1']),
        'output/030_norm/kmer_plot.pdf',
        'output/011_kraken/kraken_report.txt',
        'output/010_trim-decon/vvul.fq.gz'

# 04 meraculous
rule meraculous:
    input:
        unpack(readset_wildard_resolver),
        config = ('output/040_meraculous/'
                  '{read_set}_k{k}_diplo{diplo}/config.txt')
    output:
        ('output/040_meraculous/{read_set}_k{k}_diplo{diplo}/'
         'meraculous_final_results/final.scaffolds.fa')
    params:
        outdir = 'output/040_meraculous/{read_set}_k{k}_diplo{diplo}',
        dmin = '0'
    threads:
        meraculous_threads
    log:
        'output/logs/040_meraculous/{read_set}_k{k}_diplo{diplo}.log'
    singularity:
        mer_container
    shell:
        'run_meraculous.sh '
        '-dir {params.outdir} '
        '-config {input.config} '
        '-cleanup_level 2 '
        '&> {log}'

rule meraculous_config:
    input:
        unpack(readset_wildard_resolver)
    output:
        config = ('output/040_meraculous/'
                  '{read_set}_k{k}_diplo{diplo}/config.txt'),
    threads:
        1
    params:
        threads = meraculous_threads
    run:
        write_config_file(
            input.fq,
            params.threads,
            meraculous_config_string,
            output.config)


# 03 normalise input
rule plot_kmer_coverage:
    input:
        hist_before = 'output/030_norm/hist.txt',
        hist_after = 'output/030_norm/hist_out.txt',
        peaks = 'output/030_norm/peaks.txt'
    output:
        plot = 'output/030_norm/kmer_plot.pdf'
    threads:
        1
    log:
        log = 'output/logs/030_norm/plot_kmer_coverage.log'
    singularity:
        r_container
    script:
        'src/plot_kmer_coverage.R'

rule bbnorm:
    input:
        fq = 'output/010_trim-decon/vvul.fq.gz'
    output:
        fq_norm = 'output/030_norm/vvul.fq.gz',
        fq_toss = 'output/030_norm/toss.fq.gz',
        hist = 'output/030_norm/hist.txt',
        hist_out = 'output/030_norm/hist_out.txt',
        peaks = 'output/030_norm/peaks.txt'
    log:
        norm = 'output/logs/030_norm/bbnorm.log'
    params:
        target = 60
    threads:
        25
    singularity:
        bbduk_container
    shell:
        'bbnorm.sh '
        'in={input.fq} '
        'threads={threads} '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target={params.target} '
        'min=5 '
        'peaks={output.peaks} '
        '2> {log.norm} '  


# 02 attempt to merge overlapping reads
rule bbmerge:
    input:
        fq = 'output/010_trim-decon/vvul.fq.gz'
    output:
        merged = 'output/020_merge/merged.fq.gz',
        unmerged = 'output/020_merge/unmerged.fq.gz',
        ihist = 'output/020_merge/ihist.txt'
    log:
        merge = 'output/logs/020_merge.log'
    threads:
        25
    singularity:
        bbduk_container
    shell:
        'bbmerge.sh '
        'threads={threads} '
        'in={input.fq} '
        'verystrict=t '
        'out={output.merged} '
        'outu={output.unmerged} '
        'ihist={output.ihist} '
        '2> {log.merge} '

# 011 run kraken on decontaminated reads
rule kraken:
    input:
        r1 = 'output/011_kraken/r1.fq.gz',
        r2 = 'output/011_kraken/r2.fq.gz',
        db = directory('data/20180917-krakendb')
    output:
        out = 'output/011_kraken/kraken_out.txt',
        report = 'output/011_kraken/kraken_report.txt'
    log:
        'output/logs/011_kraken/kraken.log'
    threads:
        50
    priority:
        10
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--report-zero-counts '
        '--use-mpa-style '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

rule split:
    input:
        fq = 'output/010_trim-decon/vvul.fq.gz'
    output:
        r1 = temp('output/011_kraken/r1.fq.gz'),
        r2 = temp('output/011_kraken/r2.fq.gz')
    threads:
        25
    log:
        'output/logs/011_kraken/reformat.log'
    singularity:
        bbduk_container
    shell:
        'reformat.sh '
        'in={input.fq} '
        'int=t '
        'out={output.r1} '
        'out2={output.r2} '
        '2> {log}'

# 01 trim and decontaminate reads
rule trim_decon:
    input:
        r1 = 'output/010_trim-decon/r1.fq',
        r2 = 'output/010_trim-decon/r2.fq'
    output:
        fq = 'output/010_trim-decon/vvul.fq.gz',
        f_stats = 'output/010_trim-decon/vvul_filter-stats.txt',
        t_stats = 'output/010_trim-decon/vvul_trim-stats.txt'
    log:
        filter = 'output/logs/010_trim-decon/filter.log',
        trim = 'output/logs/010_trim-decon/trim.log',
        repair1 = 'output/logs/010_trim-decon/repair1.log',
        repair2 = 'output/logs/010_trim-decon/repair2.log'
    params:
        filter = bbduk_ref,
        trim = bbduk_adaptors
    threads:
        25
    singularity:
        bbduk_container
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '2> {log.repair1} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.f_stats} '
        '2> {log.filter} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'repair.sh '
        'in=stdin.fastq '
        'out={output.fq} '
        '2> {log.repair2} '

rule join_reads:
    input:
        r1 = r1_raw,
        r2 = r2_raw
    output:
        r1 = temp('output/010_trim-decon/r1.fq'),
        r2 = temp('output/010_trim-decon/r2.fq')
    singularity:
        bbduk_container
    shell:
        'cat {input.r1} > {output.r1} & '
        'cat {input.r2} > {output.r2} & '
        'wait'

