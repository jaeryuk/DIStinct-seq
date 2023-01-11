
##Required packages
#python3
#cutadapt
#seqkit
#bwa
#samtools
#picard

##############################

#Arg1 = thread
#Arg2 = sample name
#Arg3 = working directory
#Arg4 = illumina adapter sequence
#Arg5 = LTR 3prime terminal sequence
#Arg6 = reference genome #fusion of vector genome and human genome
#Arg7 = path of insertion_site_detection.py

##Step1: trim illumina adapter
fastp -i ${3}/${2}_1.fq \
  -I ${3}/${2}_2.fq \
  -o ${3}/${2}_1.trimmed.fq \
  -O ${3}/${2}_2.trimmed.fq \
  -h ${3}/${2}.fastp.html \
  -j ${3}/${2}.fastp.json \
  -w ${1}

##Step2: extract LTR containing reads
cat ${3}/${2}_1.trimmed.fq | seqkit grep -s -r -i -p ^${5} > ${3}/${2}_1.extracted.fq
seqkit seq ${3}/${2}_1.extracted.fq -n -i > ${3}/${2}_IDs.txt
cat ${3}/${2}_2.trimmed.fq | seqkit grep -f ${3}/${2}_IDs.txt > ${3}/${2}_2.extracted.fq

##Step3: trimming LTR sequences
cutadapt -j ${1} -g "^${5}" -o ${3}/${2}_1.extracted.trimmed.fq ${3}/${2}_1.extracted.fq
cutadapt -j ${1} -g "^${5}" -o ${3}/${2}_2.extracted.trimmed.fq ${3}/${2}_2.extracted.fq

##Step4: mapping to vector-human reference genome and sort
bwa mem -t ${1} -M ${6} ${3}/${2}_1.extracted.trimmed.fq ${3}/${2}_2.extracted.trimmed.fq \
> ${3}/${2}.mapped.bam

samtools sort -@ ${1} -O bam -o ${3}/${2}.mapped.sorted.bam ${3}/${2}.mapped.bam

##Step5: remove duplicated reads
picard MarkDuplicates I=${3}/${2}.mapped.sorted.bam O=${3}/${2}.dedup.bam M=${3}/${2}.dedup.metrics.txt REMOVE_DUPLICATES=true

##Step6: extract mapping quality >=20 & properly paired reads & primary alignment
samtools view -b -f 0x2 -F 0X100 ${3}/${2}.dedup.bam > ${3}/${2}.filtered.bam

##Step7: extract first read pair
samtools view -f 0x40 ${3}/${2}.filtered.bam > ${3}/${2}.first_pair.sam

##Step8: exclude reads mapped to lenti vector sequenceq
cat ${3}/${2}.first_pair.sam | awk '$3 !~ /^lenti/' > ${3}/${2}.lenti_excluded.sam

##Step9: exclude PCR recombination using fragment length and CIGAR string & exclude alt and decoy chromosome
cat ${3}/${2}_lenti_excluded.sam | awk '$7~/=/'| awk '($9 <= 2000) || (-$9 <= 2000)' \
| awk '($6~/^[0-9]*M$/) || ($6~/^[0-9]*M([0-9]*I|[0-9]*D)[0-9]*(M)$/)' | awk ' ($9 >0)' |  awk '($9 < 151)' \
> ${3}/${2}_temp1.sam

cat ${3}/${2}_lenti_excluded.sam | awk '$7~/=/'| awk '($9 <= 2000) || (-$9 <= 2000)' | \
awk '($6~/^[0-9]*M$/) || ($6~/^[0-9]*M([0-9]*I|[0-9]*D)[0-9]*(M)$/)' | awk ' ($9 <0)' |  awk '($9 > -151)' \
> ${3}/${2}_temp2.sam

cat ${3}/${2}_temp1.sam | awk '($12~/^(MC:Z:)(*)([0-9]*M)$/) || ($13~/^(MC:Z:)(*)([0-9]*M)$/)' \
> ${3}/${2}_temp3.sam

cat ${3}/${2}_temp2.sam | awk '($12~/^(MC:Z:[0-9]*M)/) || ($13~/^(MC:Z:[0-9]*M)/)' \
> ${3}/${2}_temp4.sam

cat ${3}/${2}_lenti_excluded.sam |  awk '$7~/=/'| awk '($9 <= 2000) || (-$9 <= 2000)' \
| awk '($6~/^([0-9]*M)/)' | awk ' ($9 >= 151) || ($9 <= -151)' \
| awk '($12~/^(MC:Z:)[0-9]*M$/) || ($12~/^(MC:Z:)[0-9]*M([0-9]*I|[0-9]*D)[0-9]*(M)$/) || ($13~/^(MC:Z:)[0-9]*M$/ ) || ($13~/^(MC:Z:)[0-9]*M([0-9]*I|[0-9]*D)[0-9]*(M)$/)' \
> ${3}/${2}_temp5.sam

cat ${3}/${2}_temp3.sam ${3}/${2}_temp4.sam ${3}/${2}_temp5.sam \
> ${3}/${2}_final.sam

##Step10: insertion site detection and convert to bed format
python3 "${7}" "${2}" "${3}"
