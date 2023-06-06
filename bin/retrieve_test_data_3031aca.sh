#!/usr/bin/env bash

ml seqkit/0.7.0

mkdir test_data

cp /output/genomic/fairGenomes/Fungus/Neonectria/ditissima/sex_na/1x/assembly_rs324p/v1/Nd324_canupilon_all.sorted.renamed.fasta \
./test_data/test_data_original.fasta

seqkit sample -p 0.8 -s 33 ./test_data/test_data_original.fasta | gzip -c > ./test_data/test_data1.fasta.gz
seqkit sample -p 0.8 -s 49 ./test_data/test_data_original.fasta > ./test_data/test_data2.fasta
seqkit sample -p 0.5 -s 22 ./test_data/test_data_original.fasta > ./test_data/test_data3.fasta
seqkit sample -p 0.5 -s 33 ./test_data/test_data_original.fasta | gzip -c > ./test_data/test_data4.fasta.gz

rm -f ./test_data/test_data_original.fasta

sequence_ids=$(zcat ./test_data/test_data1.fasta.gz | grep '^>' | awk '{print substr($1, 2)}' | sort -u)
cat /output/genomic/fairGenomes/Fungus/Neonectria/ditissima/sex_na/1x/assembly_rs324p/v1/augustus.hints.fixed.gff3 \
| grep -wFf <(echo "$sequence_ids") \
> ./test_data/test_data1.gff3

gzip ./test_data/test_data1.gff3

sequence_ids=$(cat ./test_data/test_data2.fasta | grep '^>' | awk '{print substr($1, 2)}' | sort -u)
cat /output/genomic/fairGenomes/Fungus/Neonectria/ditissima/sex_na/1x/assembly_rs324p/v1/augustus.hints.fixed.gff3 \
| grep -wFf <(echo "$sequence_ids") \
> ./test_data/test_data2.gff3

cat ./test_data/test_data1.fasta.gz | gzip -cdf | grep ">*chr" | head -3 | sed 's/>//g' | awk '{print $1, "h1_"NR}' OFS="\t" > ./test_data/test_data1.seq.list
cat ./test_data/test_data2.fasta | grep ">*chr" | tail -2 | sed 's/>//g' | awk '{print $1, "h2_"NR}' OFS="\t" >  ./test_data/test_data2.seq.list
cat ./test_data/test_data3.fasta | grep ">*chr" | head -5 | tail -2 | sed 's/>//g' | awk '{print $1, "GA_"NR}' OFS="\t" >  ./test_data/test_data3.seq.list
cat ./test_data/test_data4.fasta.gz | gzip -cdf | grep ">*chr" | tail -3 | sed 's/>//g' | awk '{print $1, "GB_"NR}' OFS="\t" >  ./test_data/test_data4.seq.list