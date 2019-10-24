#!/bin/bash

dossier_reads_brut=$1
dossier_sortie=$2

list_fastq=$(ls fastq/)
mkdir $dossier_sortie
mkdir $dossier_sortie/fastqc
mkdir $dossier_sortie/trimmed/
mkdir $dossier_sortie/merged/
mkdir $dossier_sortie/amplicon
mkdir $dossier_sortie/dereplication/
mkdir $dossier_sortie/chimera/
mkdir $dossier_sortie/otu/
touch $dossier_sortie/amplicon/amplicon.fasta
gunzip $dossier_reads_bruts/*.fastq.gz

for i in $list_fastq; do
    fastqc -o $dossier_sortie/fastqc --noextract fastq/$i
done

for file in $(ls ./$dossier_reads_brut/*_R1.fastq); do
    echo $file
    nameR1="$file"
    nameR2=$(echo $file|sed "s:R1:R2:g")
    java -jar soft/AlienTrimmer.jar -if ./$dossier_reads_brut/$nameR1 -ir ./$dossier_reads_brut/$nameR2 -c databases/contaminants.fasta -q 20 -of $dossier_sortie/trimmed/trim$nameR1 -or $dossier_sortie/trimmed/trim$nameR1
done

for file in $(ls ./$dossier_sortie/trimmed/*_R1.fastq); do
    vsearch --fastq_mergepairs $file --reverse ${file/_R1/_R2} --fastaout $dossier_sortie/merged/"$file".fasta --label_suffix ";sample=$file;"
done

sed "s: ::g" $dossier_sortie/merged/* >> $dossier_sortie/amplicon/amplicon.fasta

vsearch --derep_fulllength $dossier_sortie/amplicon/amplicon.fasta --output $dossier_sortie/dereplication/dereplication.fasta --minuniquesize 10
vsearch --uchime_denovo  $dossier_sortie/dereplication/dereplication.fasta --nonchimeras $dossier_sortie/chimera/chimera_free.fasta
vsearch --cluster_size $dossier_sortie/chimera/chimera_free.fasta --centroids $dossier_sortie/otu/otu.fasta --id 0.97 --clusterout_id --sizeout --relabel "OTU_"
vsearch --usearch_global $dossier_sortie/amplicon/amplicon.fasta --db $dossier_sortie/otu/otu.fasta --otutabout $dossier_sortie/otu/otu_table.txt --id 0.97 --sizeout
vsearch --usearch_global $dossier_sortie/otu/otu.fasta --db databases/mock_16S_18S.fasta --id 0.9 --top_hits_only --userfields query+target --userout $dossier_sortie/otu/otu_table_annoted.txt
