
for file in data/fastq/*_1.fq.gz
do
mv "$file" "${file/_1.fq.gz/.sra_1.fastq.gz}"
done

for file in data/fastq/*_2.fq.gz
do
mv "$file" "${file/_2.fq.gz/.sra_2.fastq.gz}"
done

