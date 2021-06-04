#!/bin/bash
### READ ME ###
echo '
Note:
Requires BBmap: https://sourceforge.net/projects/bbmap/
Requires prodigal, which can be installed using:
conda install -c bioconda prodigal
Requires bowtie2, which can be installed:
conda install -c bioconda bowtie2
Requires pigz to be installed:
sudo apt-get install -y pigz
'
### PARAMETERS ###
#Full directory with contigs from assembly, no '/' at the end:
HOME_DIR=/media/julio/Storage1/CRC/Summer_2019_colon_cancer_data/functional

#Full directory with per sample fastq.gz files, no '/' at the end:
FASTQ_DIR=/media/julio/Storage1/CRC/Summer_2019_colon_cancer_data/QF_nonhuman_seqs

#How your fastq.gz files end. E.g. for sample123_1.fastq.gz it would be _1.fastq.gz
FILE_END_FORWARD=.1.fastq.gz
FILE_END_REVERSE=.2.fastq.gz

#File name of contigs:
CONTIGS=megahit_contigs.fa

#Full BBmap directory, no '/' at the end:
BBMAP=/media/julio/Storage1/Software/bbmap

#Number of cores/threads:
NUMCORE=32

#Minimum length of contigs:
MINLEN=2500

#Minimum read count of sample (paired end):
MINCOUNT=100000

### "Code" ###
cd $HOME_DIR

$BBMAP/reformat.sh in=$CONTIGS out=Reformated_contigs_$MINLEN.fa minlength=$MINLEN overwrite=true
prodigal -p meta -i Reformated_contigs_$MINLEN.fa -d prodigal_genes.fna -a prodigal_genes.faa
mkdir bt2_db pile_up sam_files
bowtie2-build --threads $NUMCORE -f prodigal_genes.fna bt2_db/prodigal_db

for i in $FASTQ_DIR/*$FILE_END_FORWARD; do basename ${i} $FILE_END_FORWARD; done > sample_list.txt

#change zcat to cat if using fastq files like a psycho or 
#comment everything out if you dont care about read counts.
while read i; do
echo -n "${i} "
echo $(zcat $FASTQ_DIR/${i}$FILE_END_FORWARD|wc -l)/4|bc;
done < sample_list.txt > read_counts.txt
awk -v var=$MINCOUNT '(NR $1) && ($2 > var)' read_counts.txt > sample_list_tmp.txt
awk '{$2=""; print $0}' sample_list_tmp.txt > sample_list.txt
rm sample_list_tmp.txt
#Comment out until here if you don't care about read counts (you should).

while read i; do
bowtie2 -x bt2_db/prodigal_db -1 $FASTQ_DIR/${i}$FILE_END_FORWARD -2 $FASTQ_DIR/${i}$FILE_END_REVERSE -p $NUMCORE -S sam_files/${i}.sam
$BBMAP/pileup.sh in=sam_files/${i}.sam out=pile_up/${i}_raw.txt rpkm=pile_up/${i}_norm.txt
pigz sam_files/${i}.sam;
done < sample_list.txt

forawk=$(ls pile_up/ | head -1)
awk '{print $1}' pile_up/$forawk > contig_table.txt

#Change $5 for raw reads (recommened) 
#or $8 for FKPM (which is only useful for within sample comparison, not recommended)
while read i; do
awk '{print $5}' pile_up/${i}_norm.txt > contig_table_tmp.txt
paste contig_table.txt contig_table_tmp.txt > contig_table_tmp2.txt && mv contig_table_tmp2.txt contig_table.txt
echo ${i} >> col_order.txt;
done < sample_list.txt

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' col_order.txt > col_order2.txt

sed -i 's/ /\t/g' col_order2.txt
sed -i 's/^/#NAME\t/' col_order2.txt
sed -i '/^[[:blank:]]*#/d;s/#.*//' contig_table.txt
cat col_order2.txt contig_table.txt > contig_table2.txt && mv contig_table2.txt contig_table.txt
rm contig_table_tmp.txt col_order2.txt col_order.txt
pigz pile_up/*.txt

