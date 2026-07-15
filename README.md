![polyT-AS-Ter](images/polyT-AS-Ter-logo.png)
# polyT-AS-Ter: Validation of splice detection tools with manually configured simulated data or long reads.

Cryptic splicing describes the occurrence of unannotated, misspliced isoforms. It is associated with numerous neurodegenerative diseases. There are numerous tools available for studying differential splicing, however they often don't recover unannotated events or don't distinguish them from annotated alternative splicing. Recently, the [SpliCeAT](https://github.com/GTK-lab/SpliCeAT) pipeline was developed to discover and quantify cryptic splicing on the event- and isoform-level.
It includes consensus voting of the differential splicing detection tools [Whippet](https://github.com/timbitz/Whippet.jl), [LeafCutter](https://github.com/davidaknowles/leafcutter) and [MAJIQ](https://majiq.biociphers.org/) and .

This pipeline, polyT-AS-Ter, was developed to benchmark SpliCeAT comprehensively, but it can be used for other tools as well.
## polyT-AS-Ter offers two validation tracks
1)  generate a simulated RNAseq experiment from a ground truth gtf containing novel transcripts, run their tool of interest and validate its output against the ground truth or 
2) validate their tool of interest's output against condition-specific gtfs generated with long reads.


### How to run the pipeline for both tracks
1) Simulated Data Workflow

00_generate_augmented_transcriptome

|

|

01_get_junctions

|

|

02_classify_junctions

|

|

03_simulate_experiment

|

|

[run tool of interest on simulated data]

|

|

04_score_junctions

The Run_simulation Snakefile runs transcriptome generation and simulation, while the Run_scoring Snakefile runs the steps required for scoring.

Alternatively, each step can be run individually with the Snakefile provided in the respective directory.

2) Long Read Experiment
Executed with a single Snakefile in the Y_isoseq_validation directory. Note that it has a seperate config file in the same directory.


### Requirements before starting
-- MOST importantly: make sure that the tool you want to test produces output that can be converted into a 1-junction/1-contig-per row format (ideally .tsv) and provides a directional, quantitative measure of junction-/contig-usage between conditions (like ΔPSI)

-- Reference transcriptome gtf

-- reference genome fasta

### Parameters to adjust in the config files
1) Simulated data
- reference fasta and gtf files
- Number of novel transcripts to be generated per event type
- Expression level of original and novel transcripts relative to baseline, per condition (allowing exact control of ΔPSI)
- Number of samples per condition
- several polyester-related parameters

2) Iso-Seq
- whether or not to keep unclassifiable events (events outside of reference annotated transcripts)


### Execution
Execute a Snakemake dry run with

```bash
snakemake -np
```

to check the parameters of the run. Once ready to run, execute

```bash
snakemake --use-conda --cores 16
```


### Standalone Scripts

Besides the Snakemake workflow, a number of standalone scripts have been used for the benchmarking of SpliCeAT. These include:
- two standalone scripts for bidirectional filtering with individual tools to generate output of hypothetical truncated SpliCeAT pipelines skipping the consensus voting
- the iCLIP target enrichment analysis workflow