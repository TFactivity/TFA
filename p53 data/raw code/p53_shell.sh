#!/bin/bash
source /home/yihanlin_pkuhpc/lustre2/wuyan/.bashrc ;
cd /home/yihanlin_pkuhpc/lustre3/wuyan/scTFA/P53/17NSMB/analysis/analyzed20200618_1 ;

# # download data from ncbi through wget
# cd /home/yihanlin_pkuhpc/lustre3/wuyan/scTFA/P53/17NSMB/data ;
# for site in `cat URL_for_RNA_seq_of_WT_only.txt`
# do
#   {
#   wget ${site} -q;
#   } &
# wait
# done

n_threads=20;

mkdir STAR_index ;
STAR  --runMode genomeGenerate --genomeDir ./STAR_index --runThreadN ${n_threads} \
--genomeFastaFiles /home/yihanlin_pkuhpc/lustre2/wuyan/resource/human_genome/gencode/GRCh38.p13.genome.fa \
--sjdbGTFfile /home/yihanlin_pkuhpc/lustre2/wuyan/resource/human_genome/gencode/gencode.v34.annotation.gtf --sjdbOverhang 73 ;

data_dir=/home/yihanlin_pkuhpc/lustre3/wuyan/scTFA/P53/17NSMB/data ;
mkdir map;

for data in `ls ${data_dir}| grep "SRR"`
do
  {
  fastq-dump ${data_dir}/${data} -O ./map ;

  java -jar /home/yihanlin_pkuhpc/lustre2/wuyan/app/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads ${n_threads} -trimlog ./map/trim.log \
  ./map/${data}.fastq \
  ./map/${data}_output.fq.gz \
  ILLUMINACLIP:/home/yihanlin_pkuhpc/lustre2/wuyan/app/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10:8:true SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36 ;

  STAR \
      --runMode alignReads \
      --outFileNamePrefix ./map/${data}_ \
      --genomeDir ./STAR_index \
      --runThreadN ${n_threads} \
      --readFilesIn ./map/${data}_output.fq.gz \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingThreadN ${n_threads} ;

  rm -f ./map/${data}.fastq ./map/${data}_output.fq.gz ;

  TPMCalculator -g /home/yihanlin_pkuhpc/lustre2/wuyan/resource/human_genome/gencode/gencode.v34.annotation.gtf \
  -b ./map/${data}_Aligned.sortedByCoord.out.bam -c 74 -k gene_name ;

} &
wait
done
