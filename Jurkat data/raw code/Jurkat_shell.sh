#!/bin/bash
source /home/yihanlin_pkuhpc/lustre2/wuyan/.bashrc ;
cd /home/yihanlin_pkuhpc/lustre3/wuyan/scTFA/Jurkat/HW/analysis/analyzed20200622_1 ;

n_threads=24 ;

mkdir STAR_index ;
STAR  --runMode genomeGenerate --genomeDir ./STAR_index --runThreadN ${n_threads} \
--genomeFastaFiles /home/yihanlin_pkuhpc/lustre2/wuyan/resource/human_genome/gencode/GRCh38.p13.genome.fa \
--sjdbGTFfile /home/yihanlin_pkuhpc/lustre2/wuyan/resource/human_genome/gencode/gencode.v34.annotation.gtf --sjdbOverhang 149 ;

data_dir=/home/yihanlin_pkuhpc/lustre3/wuyan/scTFA/Jurkat/HW/data ;
mkdir map;

for data in `ls ${data_dir}| grep "_combined_R1.fastq.gz"`
do
  {
  sample_name=${data/_combined_R1.fastq.gz/} ;

  java -jar /home/yihanlin_pkuhpc/lustre2/wuyan/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads ${n_threads} -trimlog ./map/trim.log \
  ${data_dir}/${sample_name}_combined_R1.fastq.gz ${data_dir}/${sample_name}_combined_R2.fastq.gz \
  ./map/${sample_name}_output_forward_paired_R1.fq.gz ./map/${sample_name}_output_forward_unpaired_R1.fq.gz \
  ./map/${sample_name}_output_reverse_paired_R2.fq.gz ./map/${sample_name}_output_reverse_unpaired_R2.fq.gz \
  ILLUMINACLIP:/home/yihanlin_pkuhpc/lustre2/wuyan/app/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36 ;

  STAR \
      --runMode alignReads \
      --outFileNamePrefix ./map/${sample_name}_ \
      --genomeDir ./STAR_index \
      --runThreadN ${n_threads} \
      --readFilesIn ./map/${sample_name}_output_forward_paired_R1.fq.gz ./map/${sample_name}_output_reverse_paired_R2.fq.gz \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingThreadN ${n_threads} ;

  rm -f ./map/${sample_name}_output_forward_paired_R1.fq.gz ./map/${sample_name}_output_forward_unpaired_R1.fq.gz \
  ./map/${sample_name}_output_reverse_paired_R2.fq.gz ./map/${sample_name}_output_reverse_unpaired_R2.fq.gz ;

  TPMCalculator -g /home/yihanlin_pkuhpc/lustre2/wuyan/resource/human_genome/gencode/gencode.v34.annotation.gtf \
  -b ./map/${sample_name}_Aligned.sortedByCoord.out.bam -c 150 -p -k gene_name ;

} &
wait
done
