#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 00:10:00
#SBATCH --mem 5G

source ReadConfig.sh $1

config=$(basename $1 | sed 's/.txt//' | sed 's/Config.//')

rm ${WORKDIR}/PreSeqCoverageDetailed.${config}.2.txt
#rm ${WORKDIR}/PreSeqCoverageDetailed.txt


while read sample
do 
echo $sample
grep -v "^TOTAL" ${WORKDIR}/${sample}_inferredcov.txt | awk -v sample=$sample '{print sample"\t"$0}' >> ${WORKDIR}/PreSeqCoverageDetailed.${config}.2.txt 
done < ${ORIDIR}/${SAMPLELIST}


#samples_fq=$(awk -v dir=$WORKDIR '{print dir$0"_reportfastqc.txt"}' ${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')
samples_qc=$(awk -v dir=$WORKDIR '{print dir"/"$0"_QC.txt"}' ${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')
#samples_lc=$(awk -v dir=$WORKDIR '{print dir$0".ForLorenzCurves.txt"}' ${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')
samples_lc2=$(awk -v dir=$WORKDIR '{print dir"/"$0"_LorenzCurve.txt"}' ${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')


SAMPLE=$(head -1 ${ORIDIR}/${SAMPLELIST})
genome_length=$(head -1 ${WORKDIR}/${SAMPLE}_genome.bed | awk '{print $4}')
#cut -f1 ${ORIDIR}/${ADAPTERS} | awk 'BEGIN { ORS = "\t" } { print }' | sed 's/\t$/\n/' | awk '{print "Sample\tADAPT\tRead Position\t"$0}' | cat - $samples_fq > ${WORKDIR}/Adapters.${config}.txt
cat $samples_qc | awk -v genome=$genome_length -F $'\t' 'BEGIN{OFS=FS; print genome"\nSample","Total reads","Yield (Gb)","Sequencing depth (X)","% Unique reads","% Unmapped reads","% Unmapped reads with WGA adapters","% Duplicated reads","% Chimeric alignments","% Chimeric paired-end reads","% Mitochondrial reads","Potential bacterial contamination","% Observed breadth","% Observed breadth normalized by sequencing depth","Sequencing depth used in Preseq inferences"}{print $0}' > ${WORKDIR}/QC.${config}.txt
#cat $samples_lc | awk 'BEGIN{print "Sample\tPositions\tCum"}{print $0}' > ${WORKDIR}/LorenzCurves.${config}.txt
awk 'BEGIN{print "Sample\tTotalBases\tExpectedCoveredBases\tLower95%\tUpper95%"}{print $0}' ${WORKDIR}/PreSeqCoverageDetailed.${config}.2.txt > ${WORKDIR}/PreSeqCoverageDetailed.${config}.txt
cat $samples_lc2 | awk 'BEGIN{print "Positions\tCum\tSample"}{print $0}' > ${WORKDIR}/LorenzCurvesByWindow.${config}.txt
rm ${WORKDIR}/PreSeqCoverageDetailed.${config}.2.txt
rm *.tab
rm *flt.bam
rm *flt.bai
rm *.bed
rm *unmapped.fasta
rm *ps.bam
rm *.mr
rm *.map
#rm *_QC.txt
rm *_all.txt
rm *metrics.txt
rm *_inferredcov.txt
