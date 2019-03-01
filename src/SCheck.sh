#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 02:00:00
#SBATCH --mem 60G

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
SAMPLE=TP-HDF-9-NS
echo $SAMPLE

echo "> Loading modules"
module purge
module load picard/2.18.14
module load gcc/6.4.0 samtools/1.9
module load gcccore/6.4.0 preseq/2.0.3
module load gcccore/6.4.0 bedtools/2.27.1


echo "> PICARD COLLECT ALIGNMENT METRICS"
java -jar $EBROOTPICARD/picard.jar \
	CollectAlignmentSummaryMetrics \
        R=${RESDIR}/${REF}.fa \
        I=${WORKDIR}/${SAMPLE}.bam \
        O=${WORKDIR}/${SAMPLE}.alignment_summary_metrics_all.txt \
        MAX_INSERT_SIZE=1000 \
        ADAPTER_SEQUENCE=null \
        ADAPTER_SEQUENCE=GCTGTCAGTTAA,TTAACTGACAGCAGGAATCCCACT,GTGAGTGATGGTTGAGGTAGTGTGGAG,CTCCACACTACCTCAACCATCACTCAC,TGTGTTGGGTGTGTTTGG,CCAAACACACCCAACACA,TGTTGTGGGTTGTGTTGG,CCAACACAACCCACAACA,TGTGTTGGGTGTGTTTGG,CCAAACACACCCAACACA \
	VALIDATION_STRINGENCY=SILENT

# Transform output format from alignment summary
awk '/^CATEGORY/ {split($0,header);n=1;next; } {if(n!=1) next; for(i=2;i<=NF;++i) printf("%s\t%s\t%s\n",$1,header[i],$i);}' ${WORKDIR}/${SAMPLE}.alignment_summary_metrics_all.txt | column -t | grep -e "^UNPAIRED" -e "^PAIR" > ${WORKDIR}/${SAMPLE}.alignment_summary_metrics.txt

# Potential contaminations
echo "> METAPHYLER"
samtools view -f 0x4 ${WORKDIR}/${SAMPLE}.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${WORKDIR}/${SAMPLE}.unmapped.fasta
~/apps/Metaphyler/MetaPhylerSRV0.115/metaphyler.pl 2 ${WORKDIR}/${SAMPLE}.unmapped.fasta ${WORKDIR}/${SAMPLE}


## NEW downsampling
raw_reads=$(grep "TOTAL_READS" ${WORKDIR}/${SAMPLE}.alignment_summary_metrics.txt | awk '{print $3}')
genome_length=$(cat ${RESDIR}/${REF}.fai | cut -f2 | awk '{sum+=$1}END{print sum}')
# Not sure following line works with IonTorrent data
read_length=$(grep "MEAN_READ_LENGTH" ${WORKDIR}/${SAMPLE}.alignment_summary_metrics.txt | awk '{print $3}')
#read_length=150
sequencing_depth=$(awk -v rawreads=$raw_reads -v genomelength=$genome_length -v readlength=$read_length 'BEGIN{print readlength*rawreads/genomelength}')
subsampling_depth=0.1
raw_bases=$(awk -v readlength=$read_length 'BEGIN{printf "%.0f\n", readlength}' | awk -v rawreads=$raw_reads '{print $0*rawreads}')

# Previous downsampling (before 18/09/18)
## Downsample file for lorenz and preseq gc_extrap
#echo "> DOWNSAMPLE"
## Careful!!!Reads removedwith cutadapt or filtered out will not be taken into account
#raw_reads=$(grep "TOTAL_READS" ${WORKDIR}/${SAMPLE}.alignment_summary_metrics.txt | awk '{print $3}')
## Get the number of bases which are aligned, soft-clipped or unmapped 
#aligned_soft_bases=$(samtools view $WORKDIR/$SAMPLE.dedup.bam | cut -f10 | awk '{total+=length}END{print total}')
## Get the number of bases which are hard-clipped
#hard_bases=$(samtools view $WORKDIR/$SAMPLE.dedup.bam | cut -f6 | grep H | sed 's/\([A-Z]\)[0-9]/\1\n/g' | grep H | sed 's/H//' | awk '{sum+=$1}END{print sum}')
#if [ -z "$hard_bases" ];then
#	hard_bases=0
#fi
#raw_bases=$(($aligned_soft_bases + $hard_bases))
#genome_length=$(cat ${RESDIR}/${REF}.fai | cut -f2 | awk '{sum+=$1}END{print sum}')
#sequencing_depth=$(bc -l <<< "scale=4; $raw_bases / $genome_length")

## Avoiding subsampling when coverage is lower than 0.1
#if [ $(echo "$sequencing_depth < $subsampling_depth" | bc) -eq 1 ]; then
if [ ! -z $(awk -v seq=$sequencing_depth -v down=$subsampling_depth 'BEGIN{if ((seq/down) < 1) {print "1"}}') ]; then
        ln -s ${WORKDIR}/${SAMPLE}.bam ${WORKDIR}/${SAMPLE}.ps.bam
        echo "Sequencing depth is lower than "${subsampling_depth}": "${sequencing_depth}". No need to downsample"
	preseq_inference=${sequencing_depth}
else
	## Calculating downsampling probability
	probability=`bc -l <<< "scale=10; $subsampling_depth / $sequencing_depth"`
	echo "Downsampling probability: "$probability
	preseq_inference=${subsampling_depth}
	## Selecting a downsampling strategy
	## Default downsampling strategy: Constant memory. Usually results in approximately PROBABILITY * input reads being retained but for very small PROBABILITIES this may not be the case
	## Chained: Attempts to provide a compromise strategy that offers some of the advantages of both the ConstantMemory and HighAccuracy strategies. Uses a ConstantMemory strategy to downsample the incoming stream to approximately the desired proportion, and then a HighAccuracy strategy to finish. Works in a single pass, and will provide accuracy close to (but often not as good as) HighAccuracy while requiring memory proportional to the set of reads emitted from the ConstantMemory strategy to the HighAccuracy strategy. Works well when downsampling large inputs to small proportions (e.g. downsampling hundreds of millions of reads and retaining only 2%. Should be accurate 99.9% of the time when the input contains >= 50,000 templates (read names). For smaller inputs, HighAccuracy is recommended instead
	if [ "$raw_reads" -lt "50000" ]
	then
	        strategy="HighAccuracy"
	        echo "Downsampling following "${strategy}" strategy. From "${sequencing_depth}"X to "${subsampling_depth}"X."
	# In case sequencing depth is higher than 2X (probability will be higher than 0.05) too much memory is required and we should avoid wasting time (although we will lose accuracy)
	#elif [ $(echo " $probability > 0.05 " | bc) -eq 1 ];then
	#        strategy="ConstantMemory"
	#        echo "Downsampling following "${strategy}" strategy. From "${sequencing_depth}"X to "${subsampling_depth}"X."
	else
	        strategy="Chained"
	        echo "Downsampling following "${strategy}" strategy. From "${sequencing_depth}"X to "${subsampling_depth}"X."
	fi
	java -jar $EBROOTPICARD/picard.jar \
		DownsampleSam \
	        INPUT=${WORKDIR}/${SAMPLE}.bam \
	        OUTPUT=${WORKDIR}/${SAMPLE}.ps.bam \
	        RANDOM_SEED=1 \
	        PROBABILITY=${probability} \
	        STRATEGY=$strategy
fi

# Predict breadth at higher sequencing depths
echo "> PRESEQ"
bam2mr -o ${WORKDIR}/${SAMPLE}.mr ${WORKDIR}/${SAMPLE}.ps.bam
preseq gc_extrap -o ${WORKDIR}/${SAMPLE}_inferredcov.txt ${WORKDIR}/${SAMPLE}.mr

# Obtain coverage stats at original depth
echo "> BEDTOOLS GENOMECOV"
bedtools genomecov \
        -ibam ${WORKDIR}/${SAMPLE}.bam | grep "^genome" > ${WORKDIR}/${SAMPLE}_genome.bed

echo "> FILTER"
# 1284: remove not primary alignments, unmapped and duplicates
samtools view -F 1284 -b ${WORKDIR}/${SAMPLE}.ps.bam > ${WORKDIR}/${SAMPLE}.flt.bam
samtools index ${WORKDIR}/${SAMPLE}.flt.bam


echo "> LORENZ CURVES PER BASE"
# new lorenz are obtained from the bed file
bedtools genomecov \
     -ibam ${WORKDIR}/${SAMPLE}.ps.bam | grep "^genome" > ${WORKDIR}/${SAMPLE}_genome.ps.bed

echo "> LORENZ CURVES PER WINDOW"
# 1284: remove not primary alignments, unmapped and duplicates
samtools view -F 1284 -b ${WORKDIR}/${SAMPLE}.ps.bam > ${WORKDIR}/${SAMPLE}.flt.bam
samtools index ${WORKDIR}/${SAMPLE}.flt.bam

## Get lorenz curves by windows (we don't need them anymore)
# Python 2.7.5
python /home/uvi/be/tpf/apps/pysamstats/scripts/pysamstats \
	-t coverage_ext_binned \
	${WORKDIR}/${SAMPLE}.flt.bam \
	--fasta ${RESDIR}/${REF}.fa \
	--max-depth=1000000000 \
	--window-size=1000000 \
	--omit-header | cut -f4 | sort -n | \
	awk -v SAMPLE=${SAMPLE} -v dir=${WORKDIR} 'BEGIN{sum=0}{sum+=$1; print sum}END{print sum > dir"TotalCumCov."SAMPLE }' > ${WORKDIR}/${SAMPLE}_pysamstats.txt

totalcumcov=$(head -1 ${WORKDIR}/TotalCumCov.${SAMPLE})
xmax=$(awk '$1 > 0 {print $0}' ${WORKDIR}/${SAMPLE}_pysamstats.txt | wc -l)
awk -v SAMPLE=${SAMPLE} -v TOTALCUMCOV=$totalcumcov -v XMAX=$xmax 'BEGIN{num=1}{ if ($1 > 0) {print num/XMAX"\t"$1/TOTALCUMCOV"\t"SAMPLE; num+=1}}' ${WORKDIR}/${SAMPLE}_pysamstats.txt > ${WORKDIR}/${SAMPLE}_LorenzCurve.txt


echo "> CREATE QC FILE"
rm ${WORKDIR}/${SAMPLE}_QC.txt

# THIS IS NOT OK!!
# Count uniquely mapped reads using  grep -v -e 'XA:Z:' (removes reads with multiple mappings) -e 'SA:Z:' (remove reads with supplementary alignments). I can consider supplementary alignments as uniquely mapped reads too, depends on me

alignments=$(samtools view -c ${WORKDIR}/${SAMPLE}.bam)
unmapped=$(samtools view -cf 4 ${WORKDIR}/${SAMPLE}.bam)
# 2308 are reads unmapped,secondary and supplementary; -F is equivalent to grep -v
primary=$(samtools view -cF 2308 ${WORKDIR}/${SAMPLE}.bam)
# 320:first read in pair, secondary alignment
primarynonunique_1=$(samtools view -f 320 ${WORKDIR}/${SAMPLE}.bam | awk '{print $1}' | sort -h | uniq -c | wc -l)
# 384:second in pair, not primary alignment
primarynonunique_2=$(samtools view -f 384 ${WORKDIR}/${SAMPLE}.bam | awk '{print $1}' | sort -h | uniq -c | wc -l)
# unique are those reads which only map to one region in the genome
unique=$((${primary}-${primarynonunique_1}-${primarynonunique_2}))
# 1024: read is PCR or optical duplicate
duplicate=$(samtools view -cf 1024 ${WORKDIR}/${SAMPLE}.bam)
MT=$(samtools idxstats ${WORKDIR}/${SAMPLE}.bam | grep -e "^MT" -e "^chrM" | awk '{print $3}')
genome=$(cut -f 5 ${WORKDIR}/${SAMPLE}_genome.bed | head -1)
breadth=$(head -1 ${WORKDIR}/${SAMPLE}_genome.bed | awk '{print (1-$5)*100}')
#genome_length=$(head -1 ${WORKDIR}/${SAMPLE}_genome.bed | awk '{print $4}')
adapter=$(grep "PCT_ADAPTER" ${WORKDIR}/${SAMPLE}.alignment_summary_metrics.txt | awk '{print $3*100}') # Although the variable name is PCT_ADAPTER, it is proportion, so I multiply by 100

# Count supplementary alignments
# without -M, a split read is flagged as 2048. With -M it will not be counted twice but as I use uniq I think it 
# should not be counted anyway twice but check!!!
# SA:Z allow us to distinguish suppl. reads even when they are not correctly flagged (2048)
chimeric_flag=$(samtools view -cf 2048 ${WORKDIR}/${SAMPLE}.bam)
bwa_mem=$(samtools view -H ${WORKDIR}/${SAMPLE}.bam | grep -c "bwa mem")
if [ $bwa_mem -gt 0 ]; then
# Select first in pair reads (64)
chimeric_1=$(samtools view -F 64 ${WORKDIR}/${SAMPLE}.bam | grep "SA:Z" | awk '{print $1}' | sort -k1,1 | uniq | wc -l)
# Select second in pair reads (128)
chimeric_2=$(samtools view -F 128 ${WORKDIR}/${SAMPLE}.bam | grep "SA:Z" | awk '{print $1}' | sort -k1,1 | uniq | wc -l)
chimeric=$(echo "$chimeric_1+$chimeric_2" | bc -l)
elif [ $chimeric_flag -eq 0 ]
then
chimeric="NA"
else
chimeric=$(samtools view -cf 2048 ${WORKDIR}/${SAMPLE}.bam)
fi


value=$( grep -c "UNPAIRED" ${WORKDIR}/${SAMPLE}.alignment_summary_metrics.txt )
if [ $value -gt 0 ]
then
	chimerapairs="NA"
else
	chimerapairs=$(grep "PCT_CHIMERAS" ${WORKDIR}/${SAMPLE}.alignment_summary_metrics.txt | awk '{print $3*100}') # Although the variable name is PCT_ADAPTER, it is proportion, so I multiply by 100
fi
class=$(awk '{if ($1 !~ "{") print $0}' ${WORKDIR}/${SAMPLE}.genus.tab | grep -v "^@" | awk '{if ($5 >= 98) print $1}' | tr -s '\n' ',' | sed 's/,$/\n/')

echo "Variables set up"

awk -v sample=${SAMPLE} \
	-v treads=${raw_reads} \
	-v tbases=${raw_bases} \
	-v seqdepth=${sequencing_depth} \
	-v talign=${alignments} \
	-v unique=${unique} \
	-v unmapped=${unmapped} \
	-v dup=${duplicate} \
	-v mt=${MT} \
	-v bact=${class} \
	-v breadth=${breadth} \
	-v genlen=${genome_length} \
	-v adapt=${adapter} \
	-v suppl=${chimeric} \
	-v chimerapairs=${chimerapairs} \
	-v preseq_inf=${preseq_inference} \
	-F $'\t' \
	'BEGIN{OFS=FS; print sample,treads,tbases/1000000000,seqdepth,unique/treads*100,unmapped/treads*100,adapt,dup/(treads-unmapped)*100,suppl/treads*100,chimerapairs,mt/(treads-unmapped)*100,bact,breadth,breadth*tbases/genlen,preseq_inf}' \
	> ${WORKDIR}/${SAMPLE}_QC.txt

echo "> FINISHED"
