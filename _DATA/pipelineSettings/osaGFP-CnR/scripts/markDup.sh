#!/bin/bash
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH --mem 8g

module purge
module load picard/2.2.4
picardPath="/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar"

for f in Bam-sub/yw-3LW-wing-aGFP-CnR-sup-rep1_dm6-sacCer3_trim_q5.bam ; do
	echo $f
	fp="${f%.*}"
	fn=`basename $fp`
	java -Xmx8g -jar $picardPath SortSam INPUT= $f OUTPUT= Bam-sub/${fn}_sorted.bam SORT_ORDER=coordinate &&
	java -Xmx8g -jar $picardPath  MarkDuplicates INPUT= Bam-sub/${fn}_sorted.bam OUTPUT= Bam-sub/${fn}_sorted_markDups.bam METRICS_FILE= Bam-sub/${fn}_pcrDups REMOVE_DUPLICATES= false ASSUME_SORTED= true TAGGING_POLICY= All
done

