# Genomics Processing Details
### General steps...

For all pipelines:

load python if necessary (on UNC's Longleaf cluster):
```
module load python/3.5.1
```

and deploy pipeline:
```
sh slurmSubmission.sh
```

each branch should have the correct config settings

**NOTE** current paths in sampleInfo.tsv files are relative to the UNC HPC.


### osaGFP 3LW Wing CUT&RUN

```
git clone git@github.com:mniederhuber/cutNrun-pipeline.git
cd cutNrun-pipeline
git checkout osaGFP-CnR
```


#### CUT&RUN Summary:
*see pipeline repo for software versions and specific settings*

1. reads aligned to a combined SaCer and Drosophila `dm6` genome using `--very-sensitive-local` (see Snakefile for other optional settings)
2. Duplicates marked and removed
3. Bam files filtered on a negative-control CUT&RUN peak list - excluding regions of reproducible high signal in both IgG and no-epitope CUT&RUN experiments ~80 regions
   - see https://github.com/mckaylabunc/CUT-RUN_ExcludeLists for additional details
4. Peaks called with `macs2` with no control and `--nolambda` 
   *--nolambda turns off local lambda calculation for peak calling and instead sets a genome wide lambda which seems to reduce the rate of noise generating false-positive peak calls*
5. Coverage tracks generated with `bedtools`

----

### WT Wing FAIRE-seq Timecourse

```
git clone git@github.com:mniederhuber/faire-pipeline.git
cd faire-pipeline
git checkout wt-FAIRE 
```

#### FAIRE Summary:

1. Single-end reads aligned to `dm6` with `--very-sensitive`
2. Bam files filtered on published `dm6` blacklist 
    - blacklist from [Boyle Lab](https://github.com/Boyle-Lab)
    - [dm6-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz)
    - Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
3. Peaks called with `macs2` and sheared genomic control with default settings
4. Coverage tracks generated with `bedtools`


### osaGFP-deGrad CUT&RUN Pupal Wing FAIRE-seq

Clone the FAIRE pipeline and checkout the osa-deGrad branch...
```
git clone git@github.com:mniederhuber/faire-pipeline.git
cd faire-pipeline
git checkout osa-deGrad
```



