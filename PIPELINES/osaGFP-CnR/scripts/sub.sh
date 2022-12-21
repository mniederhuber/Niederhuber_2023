#!/bin/bash
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH --mem=10g

if [ ! -d subsample/ ];
then
				mkdir subsample
fi

module purge
module load samtools

for f in Bam/yw*rep1*q5.bam
do
	fp="${f%.*}"
	fn=`basename $fp`
	echo $fn
	echo "subsample 10%..."
	samtools view -s 0.1 -b $f > "subsample/${fn}_10.bam" 
	echo "subsample 20%..."
	samtools view -s 0.2 -b $f > "subsample/${fn}_20.bam" 
	echo "subsample 30%..."
	samtools view -s 0.3 -b $f > "subsample/${fn}_30.bam" 
	echo "subsample 40%..."
	samtools view -s 0.4 -b $f > "subsample/${fn}_40.bam" 
	echo "subsample 50%..."
	samtools view -s 0.5 -b $f > "subsample/${fn}_50.bam" 
	echo "subsample 60%..."
	samtools view -s 0.6 -b $f > "subsample/${fn}_60.bam" 
	echo "subsample 70%..."
	samtools view -s 0.7 -b $f > "subsample/${fn}_70.bam" 
	echo "subsample 80%..."
	samtools view -s 0.8 -b $f > "subsample/${fn}_80.bam" 
	echo "subsample 90%..."
	samtools view -s 0.9 -b $f > "subsample/${fn}_90.bam" 
done
