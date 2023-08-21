#! bin/bash/

set -e
cd ~/variant_calling/results

genome=~/variant_calling/data/ref_genome/KJ660346.2.fasta

bwa index $genome

# makes directories (should already be made)
# mkdir -p sam bam bcf vcf

for fq1 in ~/variant_calling/data/trimmed_fastq/*_1.trim.fastq.gz
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.fastq.gz)
    echo "base name is $base"

    fq1=~/variant_calling/data/trimmed_fastq/${base}_1.trim.fastq.gz
    fq2=~/variant_calling/data/trimmed_fastq/${base}_2.trim.fastq.gz
    sam=~/variant_calling/results/sam/${base}.aligned.sam
    bam=~/variant_calling/results/bam/${base}.aligned.bam
    sorted_bam=~/variant_calling/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/variant_calling/results/bcf/${base}_raw.bcf
    variants=~/variant_calling/results/bcf/${base}_variants.vcf
    final_variants=~/variant_calling/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done
