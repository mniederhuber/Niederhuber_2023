#!/bin/bash
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem=8g

(echo {
for f in Bam/*dupsMarked.bam; do
	fp="${f%.*}"
	fn=`basename $fp`
	jstring='"%s":{\n\t'
	printf "$jstring" "$fn" 
	for sub in subsample/${fn}_*.bam; do
		fp="${sub%.*}"
		fn=`basename $fp`
		frac=`echo $fn | awk -F'[_.]' '{print $6}'`
		j2='"%s":'
		printf "$j2" "$frac"
		flag=$(samtools flagstat $sub -O json)
		echo "$flag"
	done
	printf "}\n"
done
echo }) > subsample/out.json
