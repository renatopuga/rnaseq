# RNAseq Analysis STAR + RSEM 



### QC report 

**qualimap**

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
--paired \
-gtf $GTF \
--java-mem-size=10G
```

**qualimap results**

```bash
#!/bin/bash

# Cabeçalho
echo -e "id_amostra\tsample_name\ttotal_alignments\tread_pairs_aligned\treads_at_junctions\tnot_aligned\texonic\tintronic\tintergenic"

# Loop nos diretórios
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

**STAR + RSEM**

```bash
GENOMEDIR="/data/apps/star/ref/human_rsem_refseq"
RSEM="/opt/RSEM"
STAR="/opt/STAR/bin/Linux_x86_64/STAR"

# make dir RNASEQ_data
mkdir -p RNASEQ_data

#ulimit
ulimit -n 10024

sample=$1

	echo "[$0] zcat $sample FASTQ files"
	zcat fastq/$sample\_*.fastq.gz > $sample.fastq

        # make dir sample output dir
        mkdir  "RNASEQ_data/$sample"
	
        # run star
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

        # run rsem
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

rm -f $sample.fastq
```

**merge table**

```bash
# definindo diretório onde as amostras *abundance.tsv
args <- commandArgs(TRUE)

#aligned
nameDir <- args[1]

# definir a coluna onde começa o sinal de expressão ou contagem
colExps <- as.numeric(args[2])

# output file
output  <- args[3]

setwd(nameDir)

# listar os nomes dos arquivos
nomeArqs <- dir()

# abrir um arquivo e armazenar seu conteúdo
umArq <- read.delim(nomeArqs[1])
#umArq <- read.delim(nomeArqs[1], skip=4, header=FALSE)

# abrir todos os arquivos e armazenar apenas o sinal de expressão
vals <- sapply(nomeArqs, function(x) read.delim(x)[,colExps])

# adicionar nomes nas colunas da variável (vals)
colnames(vals) <- nomeArqs

# adicionar todas as colunas informativas
tb <- data.frame("id1"=umArq[,1],"id2"=umArq[,2],vals)

write.table(tb, paste("../merge-STAR-RSEM-",output,".tsv",sep=""), quote=FALSE, row.names=FALSE, sep='\t')
```
