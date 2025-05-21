# RNAseq Analysis STAR + RSEM 

RNA-seq data processing
Raw FASTQ files were aligned to the human reference genome (GRCh38.p14) using STAR (v2.7.10a) with default parameters. Gene-level expression was quantified with RSEM (v1.3.3) using annotation files from NCBI.

Post-alignment quality control was performed with Qualimap (v2.3), generating summary statistics including mapping rate, exon/intron distribution, and coverage metrics. QC reports were aggregated using built-in scripts. All inputs, parameters, and intermediate steps are fully documented in the repository to ensure computational reproducibility.


## QC report 

**qualimap - single end**

```bash
BAM="$1"
GTF="/data/apps/star/GCF_000001405.40_GRCh38.p14_genomic.gtf"
qualimap="/opt/qualimap_v2.3/qualimap"
output=$(basename $1 .bam)


export JAVA_OPTS="-Djava.io.tmpdir=/scratch/tmp -Xmx10G"

$qualimap rnaseq \
-outdir qc_$output \
-a proportional \
-bam $BAM \
-gtf $GTF \
--java-mem-size=10G
```

**qualimap results**

```bash
#!/bin/bash

echo -e "id_amostra\tsample_name\ttotal_alignments\tread_pairs_aligned\treads_at_junctions\tnot_aligned\texonic\tintronic\tintergenic"

for d in qc_*/ ; do
    id_amostra="${d#qc_}"
    id_amostra="${id_amostra%/}"
    qc_file="$d/rnaseq_qc_results.txt"

    if [[ -f "$qc_file" ]]; then
        sample_name=$(grep "bam file" "$qc_file" | sed -E 's/.*\/(.*)\.bam/\1/')
        total_alignments=$(grep "total alignments" "$qc_file" | sed -E 's/.*= *//; s/,//g')
        read_pairs_aligned=$(grep "read pairs aligned" "$qc_file" | sed -E 's/.*= *//; s/,//g')
        reads_at_junctions=$(grep "reads at junctions" "$qc_file" | sed -E 's/.*= *//; s/,//g')
        not_aligned=$(grep "not aligned" "$qc_file" | sed -E 's/.*= *//; s/,//g')
        exonic=$(grep "exonic" "$qc_file" | sed -E 's/.*= *([0-9,]+).*/\1/; s/,//g')
        intronic=$(grep "intronic" "$qc_file" | sed -E 's/.*= *([0-9,]+).*/\1/; s/,//g')
        intergenic=$(grep "intergenic" "$qc_file" | sed -E 's/.*= *([0-9,]+).*/\1/; s/,//g')

        echo -e "${id_amostra}\t${sample_name}\t${total_alignments}\t${read_pairs_aligned}\t${reads_at_junctions}\t${not_aligned}\t${exonic}\t${intronic}\t${intergenic}"
    fi
done
```

## STAR + RSEM (single-end)

run the script
```bash
sh star-rsem.sh MY_SAMPLE

```

**code bash**
```bash
GENOMEDIR="/data/apps/star/ref/human_rsem_refseq"
RSEM="/opt/RSEM"
STAR="/opt/STAR/bin/Linux_x86_64/STAR"

mkdir -p RNASEQ_data

#ulimit
ulimit -n 10024

sample=$1

echo "[$0] zcat $sample FASTQ files"
zcat fastq/$sample\_*.fastq.gz > $sample.fastq

mkdir  "RNASEQ_data/$sample"
	
$STAR --genomeDir $GENOMEDIR \
        --readFilesIn $sample.fastq \
        --limitBAMsortRAM 80000000000 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --twopassMode Basic \
        --outFilterMultimapNmax  4 \
        --quantMode TranscriptomeSAM \
        --runThreadN 16 \
        --outFileNamePrefix "RNASEQ_data/$sample/";

# delete tmp fastq
rm -f $sample.fastq

mkdir  "RNASEQ_data/rsem.$sample";
 
$RSEM/rsem-calculate-expression \
	--bam \
	--no-bam-output \
	--alignments \
	-p 16 \
	--fragment-length-mean 250 \
	--fragment-length-sd 50 \
	RNASEQ_data/$sample/Aligned.toTranscriptome.out.bam \
	$GENOMEDIR \
	RNASEQ_data/rsem.$sample/rsem


```

## Expression Table (Expected_count, TPM and FPKM)

run the script
```bash
sh call-merge-files.sh
```

```bash
#!/bin/bash


cd RNASEQ_data || exit

# by gees
mkdir -p gene-level
cd gene-level || exit
ls -d1 ../rsem.* | gawk '{print("ln -s",$1"/rsem.genes.results",gensub("../rsem.","","g",$1))}' | sh
cd ..

# array
declare -A typeFiles
typeFiles[5]="expected_count"
typeFiles[6]="tpm"
typeFiles[7]="fpkm"

# 5 = expected_count, 6=tpm and 7=fpkm
for i in 5 6 7
do
    R --file=../merge.files.R --args gene-level $i gene-level-${typeFiles[$i]}
done
```

**code R (merge-files.R)**

```R
args <- commandArgs(TRUE)
nameDir <- args[1]
colExps <- as.numeric(args[2])
output  <- args[3]

setwd(nameDir)

nomeArqs <- dir()

umArq <- read.delim(nomeArqs[1])
vals <- sapply(nomeArqs, function(x) read.delim(x)[,colExps])
colnames(vals) <- nomeArqs
tb <- data.frame("id1"=umArq[,1],"id2"=umArq[,2],vals)

write.table(tb, paste("../merge-STAR-RSEM-",output,".tsv",sep=""), quote=FALSE, row.names=FALSE, sep='\t')
```

### Fake Pool Generator for RNA-Seq FASTQ Files

This Perl script automates the generation of in silico technical replicates ("fake pools") for RNA-Seq data using pre-processed FASTQ files. It is designed to simulate pooled sequencing samples by subsampling a fixed number of reads from multiple individual FASTQ files belonging to the same experimental group.

**Main Features:**
Group Configuration:
The %groups hash defines which samples (E0155, E3060, E3732, E4915) belong to each group (e.g., 0, 10, 20, 50).

**Parameters:**

- reads_per_sample: number of reads to sample from each input file (7.5 million).
- replicates: number of fake technical replicates to generate per group (5).
- seed_base: base seed for random subsampling (incremented per replicate to ensure reproducibility).

**Step 1 – FASTQ Lane Concatenation:**

- Concatenates the paired lane files (L001 and L002) for each sample using zcat.
- Merged files are stored in the concat/ directory.

**Step 2 – Fake Pool Generation:**

- For each group and each replicate, a new output directory is created under fake_pools/.
- For each sample in the group:
- Subsamples a fixed number of reads using seqtk sample with a reproducible random seed.
- Subsampled reads from all samples in the group are combined into a single gzip-compressed FASTQ file representing the fake pool replicate.

**Directories:**

- concat/: holds concatenated FASTQ files.
- fake_pools/: stores output fake pool files grouped by group and replicate.
- logs/: placeholder directory (not used in this version, but reserved for future logging).

**Dependencies:**
- Perl (with File::Path and File::Basename)
-  seqtk for random read subsampling
- zcat and gzip for file processing

**code**

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename;

# CONFIGURAÇÃO
my %groups = (
    '0'  => [qw(E0155 E3060 E3732 E4915)],
    '10' => [qw(E0155 E3060 E3732 E4915)],
    '20' => [qw(E0155 E3060 E3732 E4915)],
    '50' => [qw(E0155 E3060 E3732 E4915)],
);

my $reads_per_sample = 7500000;
my $replicates = 5;
my $seed_base = 42;

my $concat_dir = "concat";
my $output_dir = "fake_pools";
my $log_dir    = "logs";

make_path($concat_dir, $output_dir, $log_dir);

# 1. Concatenar L001 + L002 por amostra
print ">> Concatenando lanes...\n";
for my $group (keys %groups) {
    for my $sample (@{$groups{$group}}) {
        my $merged = "$concat_dir/${sample}-${group}.fastq.gz";
        next if -e $merged;

        print "  + $sample-$group\n";
        my @files = glob("${sample}-${group}_S*_L001_R1_001.fastq.gz ${sample}-${group}_S*_L002_R1_001.fastq.gz");

        unless (@files) {
            warn "  ! Arquivos não encontrados para $sample-$group\n";
            next;
        }

        system("zcat @files | gzip > $merged") == 0 or die "Erro ao concatenar $sample-$group com zcat\n";
	#system("zcat @files > $merged") == 0 or die "Erro ao concatenar $sample-$group\n";
    }
}

# 2. Criar os fake pools
print ">> Create fake pools...\n";
for my $group (keys %groups) {
    for my $rep (1..$replicates) {
        my $pool_dir = "$output_dir/group$group/rep$rep";
        make_path($pool_dir);
        my $pool_fastq = "$pool_dir/fake_pool_group${group}_rep${rep}.fastq";
        open my $fh_out, "|-", "gzip > $pool_fastq.gz" or die "Erro ao abrir $pool_fastq.gz: $!\n";

        print "  > Grupo $group - Réplica $rep\n";

        for my $sample (@{$groups{$group}}) {
            my $input_fastq = "$concat_dir/${sample}-${group}.fastq.gz";
            my $seed = $seed_base + $rep;

            print "    + Subamostrando $sample-$group com seed $seed\n";
            my $cmd = "seqtk sample -s$seed $input_fastq $reads_per_sample";
            open my $fh_in, "-|", $cmd or die "Erro ao rodar $cmd: $!\n";
            while (<$fh_in>) {
                print $fh_out $_;
            }
            close $fh_in;
        }

        close $fh_out;
    }
}

print "Fake pools done!\n";

```
