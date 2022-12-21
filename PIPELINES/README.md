# Genomics Processing Details

### osaGFP 3LW Wing CUT&RUN

```
git clone git@github.com:mniederhuber/cutNrun-pipeline.git
cd cutNrun-pipeline
git checkout osaGFP-CnR
```

load python if necessary (on UNC's Longleaf cluster):
```
module load python/3.5.1
```

and deploy pipeline:
```
sh slurmSubmission.sh
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

### deets

| sample      | .fastq.gz | trim.bam | trim_q5.bam | trim_q5_sorted_dupsRemoved.bam |
|-------------|-----------|----------|-------------|--------------------------------|
| wt_3LW-Rep1 | 11269118  | 10626209 | 3288337     | 3102225                        |
| wt_3LW-Rep2 | 11208318  | 10806584 | 1826897     | 1756857                        |
| wt_3LW-Rep3 | 12405263  | 10832145 | 6128307     | 5530457                        |
| wt_6h-Rep1  | 11339030  | 11105912 | 1620914     | 1477209                        |
| wt_6h-Rep2  | 11569572  | 11464756 | 717300      | 682811                         |
| wt_6h-Rep3  | 9362182   | 8112346  | 6752204     | 5094970                        |
| wt_24h-Rep1 | 8444795   | 7526179  | 6373406     | 5732635                        |
| wt_24h-Rep2 | 11802657  | 10845048 | 9364286     | 8209666                        |
| wt_24h-Rep3 | 9261271   | 8345949  | 7164934     | 6112607                        |
| wt_44h-Rep1 | 2947053   | 2707749  | 1574226     | 1500091                        |
| wt_44h-Rep2 | 3373275   | 2994911  | 2575318     | 2424838                        |
| wt_44h-Rep3 | 3211134   | 2849991  | 2403161     | 2269829                        |

---

### osaGFP-deGrad CUT&RUN Pupal Wing FAIRE-seq
