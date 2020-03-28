# rnaQUAST 2.0.0 manual

1\. [About rnaQUAST](#sec1)  
2\. [Installation & requirements](#sec2)  
    2.1\. [General requirements](#sec2.1)  
    2.2\. [Software for _de novo_ quality assessments](#sec2.2)  
    2.3\. [Read alignment software](#sec2.3)  
3\. [Options](#sec3)  
    3.1\. [Input data options](#sec3.1)  
    3.2\. [Basic options](#sec3.2)  
    3.3\. [Advanced options](#sec3.3)  
4\. [Understanding rnaQUAST output](#sec4)  
    4.1\. [Reports](#sec4.1)  
    4.2\. [Detailed output](#sec4.2)  
    4.3\. [Plots](#sec4.3)  
5\. [Citation](#sec5)  
6\. [Feedback and bug reports](#sec6)  

<a name="sec1"></a>

## 1 About rnaQUAST

rnaQUAST is a tool for evaluating RNA-Seq assemblies using reference genome and gene database. In addition, rnaQUAST is also capable of estimating gene database coverage by raw reads and _de novo_ quality assessment using third-party software.

rnaQUAST version 2.0.0 was released under GPLv2 on February 8th, 2020 and can be downloaded from [http://cab.spbu.ru/software/rnaquast/](http://cab.spbu.ru/software/rnaquast/) or [https://github.com/ablab/rnaquast/releases](https://github.com/ablab/rnaquast/releases).

**For impatient people:**  

*   You will need Python, [gffutils](https://pythonhosted.org/gffutils/installation.html), [matplotlib](http://matplotlib.org/) and [joblib](https://joblib.readthedocs.io/en/latest/). Also you will need [GMAP](http://research-pub.gene.com/gmap/) (or [BLAT](http://hgwdev.cse.ucsc.edu/~kent/exe/)) and [BLASTN](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) installed on your machine and added to the `$PATH` variable.

*   You may also install rnaQUAST via conda

         conda install -c bioconda rnaquast

*   To verify your installation run

         python rnaQUAST.py --test 

*   To run rnaQUAST on your data use the following command

         python rnaQUAST.py \
        --transcripts /PATH/TO/transcripts1.fasta /PATH/TO/ANOTHER/transcripts2.fasta [...] \
        --reference /PATH/TO/reference.fasta --gtf /PATH/TO/gene_coordinates.gtf

<a name="sec2"></a>

## 2 Installation & requirements

<a name="sec2.1"></a>

### 2.1 General requirements

rnaQUAST can be installed via conda:

        conda install -c bioconda rnaquast

If you wish to run rnaQUAST from [the release archive](https://github.com/ablab/rnaquast/releases) you need:  

*   Python3 or Python2 (2.5+)
*   [matplotlib](http://matplotlib.org/) python package
*   [joblib](https://joblib.readthedocs.io/en/latest/) python package
*   [gffutils](https://pythonhosted.org/gffutils/installation.html) python package (needs [biopython](http://biopython.org))
*   [NCBI BLAST+ (blastn)](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
*   [GMAP](http://research-pub.gene.com/gmap/) (or [BLAT](http://hgwdev.cse.ucsc.edu/~kent/exe/)) aligner

rnaQUAST still works under Python2 (2.5+), but since Python2 is outdated, its support is not maintained since version 2.0.

Note, that due to the limitations of `BLAT`, in order to work with reference genomes of size more than 4 Gb a `pslSort` is also required.

Paths to `blastn` and `GMAP` (or `BLA`T) should be added to the `$PATH` environmental variable. To check that everything is installed correctly we recommend to run:  

    python rnaQUAST.py --test

Note that `gffutils` is used to complete gene coordinates in case of missing transcripts / genes records. For more information, see [advanced options](#sec3.3).<a name="sec2.2"></a>

### 2.2 Software for _de novo_ quality assessment

When reference genome and gene database are unavailable, we recommend to run [BUSCO](http://busco.ezlab.org/) and [GeneMarkS-T](http://topaz.gatech.edu/GeneMark/) in rnaQUAST pipeline.

**BUSCO requirements**

BUSCO allows to detect core genes in the assembled transcripts. To use it you should install [BUSCO v3](http://busco.ezlab.org/), [tblastn](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), [HMMER](http://hmmer.janelia.org/) and [transeq](ftp://emboss.open-bio.org/pub/EMBOSS/) and add these tools to the `$PATH` variable.

Depending on the species you wish to assess, you should download the appropriate lineage-specific profile libraries: Metazoa, Eukaryota, Arthropoda, Vertebrata, Fungi, or Bacteria from [http://busco.ezlab.org](http://busco.ezlab.org) and provide it to rnaQUAST with `--busco_lineage` option.

**GeneMarkS-T requirements**

[GeneMarkS-T](http://topaz.gatech.edu/GeneMark/) allows to predict genes in the assembled transcripts without reference genome. If you wish to use it in rnaQUAST pipeline, GeneMarkS-T should be properly installed and added to the `$PATH` variable.

<a name="sec2.3"></a>

### 2.3 Read alignment software

rnaQUAST is also capable of calculating various statistics using raw reads (e.g. database coverage by reads). To obtain them you need to install [STAR](https://github.com/alexdobin/STAR) aligner and add it to the `$PATH` variable. To learn more see [input options](#readopts).

<a name="sec3"></a>

## 3 Options

<a name="sec3.1"></a>

### 3.1 Input data options

To run rnaQUAST you need to provide either FASTA files with transcripts (recommended), or align transcripts to the reference genome manually and provide the resulting PSL files.

`-r <REFERENCE>, --reference <REFERENCE>`  
    Single file with reference genome containing all chromosomes/scaffolds in FASTA format (preferably with `*.fasta, *.fa, *.fna, *.ffn or *.frn` extension) OR  
    **`*.txt`** file containing the one-per-line list of FASTA files with reference sequences.

`--gtf <GENE_COORDINATES>`  
    File with gene coordinates in GTF/GFF format (needs information about parent relations). We recommend to use files downloaded from [GENCODE](http://www.gencodegenes.org/) or [Ensembl](ftp://ftp.ensembl.org/pub/).

`--gene_db <GENE_DB>`  
    Path to the gene database generated by gffutils. The database is created during the first run. This option is not compatible with `--gtf` option. We recommend to use this option once the database is created in order to speed up the run.

`--gmap_index <INDEX FOLDER>,`  
    Folder containing pre-built GMAP index for the reference genome. Using previously constructed index decreases running time. Note, that you still need to provide the reference genome that was used for index construction when this option is used.

`-c <TRANSCRIPTS ...>, --transcripts <TRANSCRIPTS, ...>`  
     File(s) with transcripts in FASTA format separated by space.

`-psl <TRANSCRIPTS_ALIGNMENT ...>, --alignment <TRANSCRIPTS_ALIGNMENT, ...>`  
     File(s) with transcript alignments to the reference genome in PSL format separated by space.

<a name="readopts"></a>

`-sam <READS_ALIGNMENT>, --reads_alignment <READS_ALIGNMENT>`  
     File with read alignments to the reference genome in SAM format.

`-1 <LEFT_READS>, --left_reads <LEFT_READS>`  
     File with forward paired-end reads in FASTQ or gzip-compressed fastq format.

`-2 <RIGHT_READS>, --right_reads <RIGHT_READS>`  
     File with reverse paired-end reads in FASTQ or gzip-compressed fastq format.

`-s <SINGLE_READS>, --single_reads <SINGLE_READS>`  
     File with single reads in FASTQ or gzip-compressed fastq format.

<a name="sec3.2"></a>

### 3.2 Basic options

`-o <OUTPUT_DIR>, --output_dir <OUTPUT_DIR>`  
     Directory to store all results. Default is `rnaQUAST_results/results_<datetime>`.

`--test`  
     Run rnaQUAST on the test data from the `test_data` folder, output directory is `rnaOUAST_test_output`.

`-d, --debug`  
     Report detailed information, typically used only for detecting problems.

`-h, --help`  
     Show help message and exit.

<a name="sec3.3"></a>

### 3.3 Advanced options

`-t <INT>, --threads <INT>`  
     Maximum number of threads. Default is min(number of CPUs / 2, 16).

`-l <LABELS ...>, --labels <LABELS ...>`  
     Name(s) of assemblies that will be used in the reports separated by space and given in the same order as files with transcripts / alignments.

`--prokaryote`  
     Use this option if the genome is prokaryotic.

`-ss, --strand_specific`  
     Set if transcripts were assembled using strand-specific RNA-Seq data in order to benefit from knowing whether the transcript originated from the + or - strand.

`--min_alignment <MIN_ALIGNMENT>`  
     Minimal alignment length to be used, default value is 50.

`--no_plots`  
     Do not draw plots (makes rnaQUAST run a bit faster).

`--blat`  
     Run with [BLAT alignment tool](http://hgwdev.cse.ucsc.edu/~kent/exe/) instead of [GMAP](http://research-pub.gene.com/gmap/).

`--busco_lineage`  
     Run [BUSCO tool](http://busco.ezlab.org/), which detects core genes in the assembly (see [Installation & requirements](#sec2)). Use this option to provide path to the BUSCO lineage data (Eukaryota, Metazoa, Arthropoda, Vertebrata or Fungi).

`--gene_mark`  
     Run with [GeneMarkS-T](http://topaz.gatech.edu/GeneMark/) gene prediction tool. Use `--prokaryote` option if the genome is prokaryotic.

`--disable_infer_genes`  
     Use this option if your GTF file already contains genes records, otherwise gffutils will fix it. Note that gffutils may work for quite a long time.

`--disable_infer_transcripts`  
     Is option if your GTF file already contains transcripts records, otherwise gffutils will fix it. Note that gffutils may work for quite a long time.

`--lower_threshold`  
     Lower threshold for x-assembled/covered/matched metrics, default: 50%.

`--upper_threshold`  
     Upper threshold for x-assembled/covered/matched metrics, default: 95%.

<a name="sec4"></a>

## 4 Understanding rnaQUAST output

In this section we describe metrics, statistics and plots generated by rnaQUAST. Metrics highlighted with **_bold italic_** are considered as the most important and are included in the short summary report (`short_report.txt`).



<a name="sec4"></a>

For the simplicity of explanation, _transcript_ is further referred to as a sequence generated by the assembler and _isoform_ denotes a sequence from the gene database. Figure below demonstrates how rnaQUAST classifies transcript and isoform sequences using alignment information. ![](fig1.png)

 <a name="sec4.1"></a>

### <a name="sec4.1">4.1 Reports</a>

The following text files with reports are contained in `comparison_output` directory and include results for all input assemblies. In addition, these reports are contained in `<assembly_label>_output` directories for each assembly separately.

**`database_metrics.txt`**  
Gene database metrics.

*   **_Genes_** / Protein coding genes – number of genes / protein coding genes
*   Isoforms / Protein coding isoforms – number of isoforms / protein coding isoforms
*   Exons / Introns – total number of exons / introns
*   Total / Average length of all isoforms, bp
*   Average exon length, bp
*   Average intron length, bp
*   **_Average_** / Maximum number of exons per isoform

</a>

<a name="sec4.1"></a><a name="readcov"></a>Coverage by reads. The following metrics are calculated only when `--left_reads`, `--right_reads`, `--single_reads` or `--sam` options are used (see [options](#readopts) for details).

*   Database coverage – the total number of bases covered by reads (in all isoforms) divided by the total length of all isoforms.
*   x%-covered genes / isoforms / exons – number of genes / isoforms / exons from the database that have at least x% of bases covered by all reads, where x is specified with `--lower_threshold / --upper_threshold` options (50% / 95% by default).

**`basic_mertics.txt`**  
Basic transcripts metrics are calculated without reference genome and gene database.

*   **_Transcripts_** – total number of assembled transcripts.
*   **_Transcripts > 500 bp_**
*   Transcripts > 1000 bp
*   Average length of assembled transcripts
*   Longest transcript
*   Total length
*   Transcript N50 – a maximal number N, such that the total length of all transcripts longer than N bp is at least 50% of the total length of all transcripts.

**`alignment_metrics.txt`**  
Alignment metrics are calculated with reference genome but without using gene database. To calculate the following metrics rnaQUAST filters all short partial alignments (see [`--min_alignment` option](#sec3.3)) and attempts to select the best hits for each transcript.

*   **_Transcripts_** – total number of assembled transcripts.
*   **_Aligned_** – the number of transcripts having at least 1 significant alignment.
*   **_Uniquely aligned_** – the number of transcripts having a single significant alignment.
*   Multiply aligned – the number of transcripts having 2 or more significant alignments. Multiply aligned transcripts are stored in `<assembly_label>.paralogs.fasta` file.
*   Misassembly candidates reported by GMAP (or BLAT) – transcripts that have discordant best-scored alignment (partial alignments that are either mapped to different strands / different chromosomes / in reverse order / too far away).
*   **_Unaligned_** – the number of transcripts without any significant alignments. Unaligned transcripts are stored in `<assembly_label>.unaligned.fasta` file.

Number of assembled transcripts = Unaligned + Aligned = Unaligned + (Uniquely aligned + Multiply aligned + Misassembly candidates reported by GMAP (or BLAT)).

Alignment metrics for non-misassembled transcripts

*   **_Average aligned fraction._** Aligned fraction for a single transcript is defined as total number of aligned bases in the transcript divided by the total transcript length.
*   Average alignment length. Aligned length for a single transcript is defined as total number of aligned bases in the transcript.
*   Average blocks per alignment. A block is defined as a continuous alignment fragment without indels.
*   Average block length (see above).
*   **_Average mismatches per transcript_** – average number of single nucleotide differences with reference genome per transcript.
*   NA50 – N50 for alignments.

<a name="misassemblies"></a>**`misassemblies.txt`**  

*   **_Transcripts_** – total number of assembled transcripts.
*   Misassembly candidates reported by GMAP (or BLAT) – transcripts that have discordant best-scored alignment (partial alignments that are either mapped to different strands / different chromosomes / in reverse order / too far away).
*   Misassembly candidates reported by BLASTN – transcripts are aligned to the isoform sequences extracted from the genome using gene database with BLASTN and then transcripts that have partial alignments to multiple isoforms are selected.
*   **_Misassemblies_** – misassembly candidates confirmed by both methods described above. Using both methods simultaneously allows to avoid considering misalignments that can be caused, for example, by paralogous genes or genomic repeats. Misassembled transcripts are stored in `<assembly_label>.misassembled.fasta` file.

**`sensitivity.txt`**  
Assembly completeness (sensitivity). For the following metrics (calculated with reference genome and gene database) rnaQUAST attempts to select best-matching database isoforms for every transcript. Note that a single transcript can contribute to multiple isoforms in the case of, for example, paralogous genes or genomic repeats. At the same time, an isoform can be covered by multiple transcripts in the case of fragmented assembly or duplicated transcripts in the assembly.

*   **_Database coverage_** – the total number of bases covered by transcripts (in all isoforms) divided by the total length of all isoforms.
*   Duplication ratio – total number of aligned bases in assembled transcripts divided by the total number of isoform covered bases. This metric does not count neither paralogous genes nor shared exons, only real overlaps of the assembled sequences that are mapped to the same isoform.
*   Average number of transcripts mapped to one isoform.
*   **_x%-assembled genes / isoforms_**/ exons – number of genes / isoforms / exons from the database that have at least x% captured by a single assembled transcript, where x is specified with `--lower_threshold / --upper_threshold` options (50% / 95% by default). 95%-assembled isoforms are stored in `<assembly_label>.95%assembled.fasta` file.
*   x%-covered genes / isoforms– number of genes / isoforms from the database that have at least x% of bases covered by all alignments, where x is specified with `--lower_threshold / --upper_threshold` options (50% / 95% by default).
*   **_Mean isoform assembly_** – assembled fraction of a single isoform is calculated as the largest number of its bases captured by a single assembled transcript divided by its length; average value is computed for isoforms with > 0 bases covered.
*   Mean isoform coverage – coverage of a single isoform is calculated as the number of its bases covered by all assembled transcripts divided by its length; average value is computed for isoforms with > 0 bases covered.
*   Mean exon coverage – coverage of a single exon is calculated as the number of its bases covered by all assembled transcripts divided by its length; average value is computed for exons with > 0 bases covered.
*   Average percentage of isoform x%-covered exons, where x is specified with `--lower_threshold / --upper_threshold` options (50% / 95% by default). For each isoform rnaQUAST calculates the number of x%-covered exons divided by the total number of exons. Afterwards it computes average value for all covered isoforms.

[BUSCO](http://busco.ezlab.org/) metrics. The following metrics are calculated only when `--busco_lineage` option is used (see [options](#sec3.3) for details).

*   **_Complete_** – percentage of completely recovered genes.
*   **_Partial_** – percentage of partially recovered genes.

[GeneMarkS-T](http://topaz.gatech.edu/GeneMark/) metrics. The following metrics are calculated when reference and gene database are not provided or `--gene_mark` option is used (see [options](#sec3.3) for details).

*   **_Genes_** – number of predicted genes in transcripts.

**`specificity.txt`**  
Assembly specificity. To compute the following metrics we use only transcripts that have at least one significant alignment and are not misassembled.

*   **_Unannotated_** – total number of transcripts that do not cover any isoform from the database. Unannotated transcripts are stored in `<assembly_label>.unannotated.fasta` file.
*   **_x%-matched_** – total number of transcripts that have at least x% covering an isoform from the database, where x is specified with `--lower_threshold / --upper_threshold` options (50% / 95% by default).
*   **_Mean fraction of transcript matched_** – matched fraction of a single transcript is calculated as the number of its bases covering an isoform divided by the transcript length; average value is computed for transcripts with > 0 bases matched.
*   Mean fraction of block matched – matched fraction of a single block is calculated as the number of its bases covering an isoform divided by the block length; average value is computed for blocks with > 0 bases matched.
*   x%-matched blocks – percentage of blocks that have at least x% covering an isoform from the database, where x is specified with `--lower_threshold / --upper_threshold` options (50% / 95% by default).
*   Matched length – total number of transcript bases covering isoforms from the database.
*   Unmatched length – total alignment length - Matched length.

**`relative_database_coverage.txt`**  
Relative database coverage metrics are calculated only when raw reads (or read alignments) are provided. rnaQUAST uses read alignments to estimate the upper bound of the database coverage and the number of x-covered genes / isoforms / exons (see [read coverage](#readcov)) and computes the following metrics:

*   **_Relative database coverage_** – ratio between transcripts database coverage and reads database coverage.
*   Relative x%-assembled genes / isoforms / exons – ratio between transcripts x%-assembled and reads x%-covered genes / isoforms / exons.
*   Relative x%-covered genes / isoforms / exons – ratio between transcripts x%-covered and reads x%-covered genes / isoforms / exons.

<a name="sec4.2"></a>

### 4.2 Detailed output

These files are contained in `<assembly_label>_output` directories for each assembly separately.  

*   `<assembly_label>.unaligned.fasta` – transcripts without any significant alignments.
*   `<assembly_label>.paralogs.fasta` – transcripts having 2 or more significant alignments.</a>
*   `<assembly_label>.misassembled.fasta` – misassembly candidates detected by methods described above. See</a> [`misassemblies.txt`](#misassemblies) description for details.
*   `<assembly_label>.correct.fasta` – transcripts with exactly 1 significant alignment that do not contain misassemblies.
*   `<assembly_label>.x%-assembled.list` – IDs of the isoforms from the database that have at least x% captured by a single assembled transcript, where x is specified by the user with an option `--upper_threshold` (95% by default).
*   `<assembly_label>.unannotated.fasta` – transcripts that do not cover any isoform from the database.

The following text file is contained in `comparison_output` directory and `<assembly_label>_output` directories for each assembly separately.

*   `reads.x%-covered.list` – IDs of the isoforms from the database that have at least x% bases covered by all reads, where x is specified with `--lower_threshold / --upper_threshold` options (50% / 95% by default).

<a name="sec4.3"></a>

### 4.3 Plots

The following plots are similarly contained in both `comparison_output` directory and `<assembly_label>_output` directories. Please note, that most of the plots represent cumulative distributions and some plots are given in logarithmic scale.

**Basic**

*   **_`transcript_length.png`_** – assembled transcripts length distribution (+ database isoforms length distribution).
*   `block_length.png` – alignment blocks length distribution (+ database exons length distribution).
*   `x-aligned.png` – transcript aligned fraction distribution.
*   `blocks_per_alignment.png` – distribution of number of blocks per alignment (+ distribution of number of database exons per isoform).
*   `alignment_multiplicity.png` – distribution for the number of significant alignment for each multiply-aligned transcript.
*   **_`mismatch_rate.png`_** – substitution errors per alignment distribution.
*   `Nx.png` – Nx plot for transcripts. Nx is a maximal number N, such that the total length of all transcripts longer than N bp is at least x% of the total length of all transcripts.
*   `NAx.png` – Nx plot for alignments.

**Sensitivity**

*   **_`x-assembled.png`_** – a histogram in which each bar represents the number of isoforms from the database that have at least x% captured by a single assembled transcript.
*   `x-covered.png` – a histogram in which each bar represents the number of isoforms from the database that have at least x% of bases covered by all alignments.
*   `x-assembled_exons.png` – a histogram in which each bar represents the number of exons from the database that have at least x% captured by a single assembled transcript.
*   `x-covered_exons.png` – a histogram in which each bar represents the number of exons from the database that have at least x% of bases covered by all alignments.
*   `alignments_per_isoform.png` – plot showing number of transcript alignments per isoform

**Specificity**

*   `x-matched.png` – a histogram in which each bar represents the number of transcripts that have at least x% matched to an isoform from the database.
*   `x-matched_blocks.png` – a histogram in which each bar represents the number of all blocks from all transcript alignments that have at least x% matched to an isoform from the database.

<a name="sec5"></a>

## 5 Citation

[Bushmanova, E., Antipov, D., Lapidus, A., Suvorov, V. and Prjibelski, A.D., 2016. rnaQUAST: a quality assessment tool for de novo transcriptome assemblies. Bioinformatics, 32(14), pp.2210-2212.](https://academic.oup.com/bioinformatics/article/32/14/2210/1743439)

<a name="sec6"></a>

## 6 Feedback and bug reports

Your comments, bug reports, and suggestions are very welcomed. They will help us to further improve rnaQUAST. If you have any troubles running rnaQUAST, please send us `logs/rnaQUAST.log` from the output directory. Address for communications: [rnaquast_support@ablab.spbau.ru](mailto:rnaquast_support@ablab.spbau.ru).
You may also submit your issue to our [GitHub repository](https://github.com/ablab/rnaquast/issues).

