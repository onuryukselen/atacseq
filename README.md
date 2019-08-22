[![Travis-ci tests:](https://travis-ci.org/dolphinnext/atacseq.svg?branch=master)](https://travis-ci.org/dolphinnext/atacseq) [![DOI:10.1101/689539](https://zenodo.org/badge/DOI/10.1101/689539.svg)](https://doi.org/10.1101/689539)

This pipeline maps reads to selected genome (by using Bowtie2), identifies enriched regions (by using MACS2) and finally creates mean count table for analysis. 
 
#### Steps:
  1. For Quality Control, we use FastQC to create qc outputs. There are optional read quality filtering (trimmomatic), read quality trimming (trimmomatic), adapter removal (cutadapt) processes available.
  2. In the sequential mapping step, Bowtie2 is used to count or filter out common reads (eg. ercc, rmsk). 
  3. Bowtie2 is used to align reads to a selected genome, and duplicates removed with Picard or Samtools,
  4. For ATAC-Seq specifically the pipeline calls accessible chromatin regions by estimating the the Tn5 transposase cut site by first positioning on the 9-th base upstream of the 5’ read end then extending by 29-bases downstream. This extension process has been shown to more accurate reflect the exact position that was accessible to transposase (Donnard et al. 2018; Buenrostro et al. 2013). Once each read has been shortened, Peaks are called using MACS2 (Zhang et al. 2008).
  5. When processing several samples together, the ATAC pipeline provide consensus peak calls by merging all peaks individually called in each samples using Bedtools (Quinlan and Hall 2010). The number of reads in each peak location are then quantified using Bedtools (Quinlan and Hall 2010) coverage function.
  6. Optionally, genome-wide Bam analysis is done by RseQC, and Picard’s CollectRNASeqMetrics program.
  7. Optionally, you can create Integrative Genomics Viewer (IGV)  and Genome Browser Files (TDF and Bigwig, respectively)
  8. As a result, ATAC-Seq pipeline generates a matrix that has the count values for each peak region and samples. This matrix can be uploaded directly to the embedded version of DEBrowser (Kucukural et al. 2019) to perform differential analysis or downloaded to perform other analysis.

#### Inputs:

  - Reads
  - ATAC-prep section

There are three fields need to be entered: output-prefix, sample-prefix, and input-prefix. Please use sample names to fill this form.For instance, to enter following files control-rep1.fastq.gz, exper-rep1.fastq.gz,  as the following.

    | output-prefix | sample-prefix | input-prefix |
    |---------------|---------------|--------------|
    | exper         |   exper-rep1  | control-rep1 |
    | control       |  control-rep1 |              |


#### Program Versions:
  - Macs2 v2.1.2
  - Bowtie2 v2.3.5
  - Bowtie v1.2.2
  - FastQC v0.11.8
  - Star v2.6.1
  - Picard v2.18.27
  - Rseqc v2.6.2
  - Samtools v1.3
  - Multiqc v1.7
  - Trimmomatic v0.39
  - Igvtools v2.5.3
  - Bedtools v2.27.1
  - Fastx_toolkit v0.0.14
  - Ucsc-wigToBigWig v366
  - Pdfbox-App v2.0.0

#### Run through DolphinNext User Interface:

To start using the dolphinnext/chipseq pipeline please go to [*DolphinNext Web page*](https://dolphinnext.umassmed.edu/index.php?np=1&id=485) and click run button.

#### Run through Command Line:

To install and start using the dolphinnext/atacseq pipeline by using command line, please follow these steps: [*Installation*](https://github.com/dolphinnext/atacseq/blob/master/docs/local.md).