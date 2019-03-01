# SCheck

--------

#### SCheck is a single-cell whole-genome amplification quality control tool  
  
## Requirements

<!-- I still do not know the version requirements -->

* [picard](http://broadinstitute.github.io/picard/) (>=2.xx)
* [samtools](http://samtools.sourceforge.net/)
* [bedtools](http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html)
* [preseq](http://smithlabresearch.org/software/preseq/)
* [metaphyler](http://metaphyler.cbcb.umd.edu/ short illumina reads version) (not installed in cesga)
* [R](https://www.r-project.org/)
* [pysamstats]()

## Usage

### 1. Make sure src is in your path.

`export PATH=path-to-SCheck/src:$PATH`

### 2. Prepare your input data:
  - Bam files with index.  
  - Genome reference with index.

### 3. Create the sample list

Text file including a bam file name per row. Do not include the *.bam* suffix. Place the sample list in the original directory.

```
AM-HDF-10-NS
TP-HDF-9-NS
```

### 4. Create the configuration file

Type `ReadConfig.sh --help` to get more info. 
Here is an example:
  
```
sample_list=Samples.txt
original_directory=/mnt/netapp2/posadalab2/uvibetpf/SCheck/test
working_directory=/mnt/netapp2/posadalab2/uvibetpf/SCheck/test
resources_directory=/mnt/netapp2/posadalab2/uvibetpf/SCheck/test
reference_name=hs37d5
```

### 5. Run the main script


```bash
main=$(sbatch --array=1-2 SCheck.sh Config.test.txt | awk '{print $4}')
```

### 6. Generate the input for the shiny app

```bash
sbatch --dependency=afterok:$main Summarize.sh Config.test.txt
```
This command will generate the three txt files you have to feed into the shiny app.

### 7. Run the shiny app

Enter R and run the following lines

```r
library(shiny)
runApp("SCheck-shiny-app")
```

Please, note that the shiny app has some R packages dependencies

<!-- still have to check whether they have been use after having modified the code. -->

* shinydashboard
* shinycssloaders
* DT
* ggplot2
* reshape2
* data.table
* htmltools
* grid
* plotly
