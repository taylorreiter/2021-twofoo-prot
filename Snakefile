rule all:
    input:
        "outputs/twofoo_k31_r1_multifasta_x/multifasta_x.cdbg_annot.csv",

rule build:
    input:
        reads="data/twofoo.fq.gz",
        conf="conf.yml"
    output: 
        "outputs/twofoo_k31_r1/catlas.csv"
    resources: 
        mem_mb = 16000,
    benchmark: "benchmarks/build.tsv"
    params: outdir = "outputs/"
    conda: "envs/spacegraphcats_prot_gather.yml"
    shell:
        "python -m spacegraphcats run {input.conf} build --nolock  --outdir={params.outdir} --rerun-incomplete" 

rule build_cdbg_list_by_record_x:
    input:
        catlas = "outputs/twofoo_k31_r1/catlas.csv",
        conf="conf.yml"
    output: 
        "outputs/twofoo_k31_r1_multifasta_x/multifasta_x.cdbg_annot.csv",
    resources: 
        mem_mb = 900000,
    threads: 1
    benchmark: "benchmarks/mutlifasta_x_gather_nbhd_all_pfam.tsv"
    params: outdir = "outputs/"
    conda: "envs/spacegraphcats_prot_gather.yml"
    shell:'''
    python -m spacegraphcats run {input.conf} build_cdbg_list_by_record_x  --nolock --outdir={params.outdir} --rerun-incomplete 
    '''

#
# akker-reads.abundtrim.gz is a collection of reads from podar data
# that maps to the Akkermansia muciniphila ATCC BAA-835 genome via bwa aln.
# "Real" data, with known answer.  There does not appear to be significant
# overlap with other genomes in the Podar data set; so, no significant strain
# variation.
#

rule download_akker:
    output:
        "data/akker-reads.abundtrim.gz"
    shell:
    	"curl -o {output} -L https://osf.io/dk7nb/download"

#
# shew-reads.abundtrim.gz is a collection of reads from podar data
# that maps to the Shewanella OS223 genome via bwa aln.  "Real" data,
# with known answer.  Note that there is significant overlap with the
# Shewanella OS185 genome; this is a data set with significant strain
# variation.
#

rule download_shew:
    output:
        "data/shew-reads.abundtrim.gz"
    shell:
    	"curl -o {output} -L 'https://osf.io/7az9p/?action=download'"

#
# twofoo, use a synthetic mixture of reads from podar data -
# the shew-reads.abundtrim.gz (mapping to Shewanella baltica OS223) and
# akker-reads.abundtrim.gz (mapping to Akkermansia muciniphila ATCC BAA-835).
# Many of the shew-reads also map to S. baltica OS185, while the akker-reads
# do not; so this is a good mixture for testing the effects of strain variation
# on catlas foo.

rule make_twofoo:
    input:
        "data/shew-reads.abundtrim.gz",
        "data/akker-reads.abundtrim.gz"
    output:
        "data/twofoo.fq.gz"
    shell:
        "gunzip -c {input} | gzip -9c > {output}"
