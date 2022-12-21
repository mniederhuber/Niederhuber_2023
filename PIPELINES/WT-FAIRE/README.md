# WT FAIRE Timecourse Details
### Pipeline

Based on McKay lab FAIRE pipeline...
```
git clone git@github.com:mniederhuber/faire-pipeline.git
cd faire-pipeline
git checkout wt-FAIRE 
```

With modifications:

- dm6 blacklist from [Boyle Lab](https://github.com/Boyle-Lab)
  - [dm6-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz)
  - Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

- `bowtie2` run with `--very-sensitve` 

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
