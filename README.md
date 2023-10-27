# CCycDB: an integrative knowledgebase to fingerprint microbial carbon cycling processes
<!-- button -->

**Please see https://ccycdb.github.io/ for more details**

CCycDB is a knowledge-based functional gene database for accurate metagenomic profiling of carbon cycling microbial communities. CCycDB contains 4,676 gene families within 6 categories. These gene families are further categorized into 45 sub-categories within sub-category I and 188 sub-categories within sub-category II with a total of 10,991,724 targeted sequences. A series of validations demonstrated that CCycDB outperformed large public orthology databases in terms of coverage, specificity, and accuracy, and can be used to accurately profile carbon cycling microbial communities in real metagenomic datasets.
Please see "CCycDB: an integrative knowledgebase to fingerprint microbial carbon cycling processes" for more details.

## Contents

- [[Usage and Examples]](#usage-and-examples)
- [Options](#options)
- [Examples](#examples)
- [Contact](#contact)

### Usage and Examples

#### Step 0. Prepare your files 

Input files supports 3 types

**1.Read-based**

Raw reads (clean reads) or Merged Reads (Forward and reverse reads were merged into longer sequences by the program e.g. PEAR).

**2.Assembly-based**

Contigs generated through metagenome assmebly by (e.g., MEGAHIT, MetaSPAdes, SPAdes).

**3.Tabular Files**

BLAST table (delimited with "\t") generated through sequences similarity searching tools (e.g., BLAST, USEARCH, DIAMOND).

---
#### Step 1. Download `#install`
```
$ wget -c https://zenodo.org/record/8324818/files/data.tar.gz?download=1

Alternatively, you can download database directly through the website or third-party download tools.

$ git clone https://github.com/ccycdb/CCycDB.PL
```
---
#### Step 2. Annotation

```
perl GetFun_CCycdb.pl [-situation read-based|assembly-based|tabular] [-wd work_directory] [-m diamond|usearch|blast] [-f filetype] [-s seqtype] [-id] [-e] [-tpm] [-norm xx] [-rs xx] [-thread xx] [-od xx]
```

### Options
```
-situation  : The situation for input files (read-based|assembly-based|tabular).

-wd  : Work directory. Ensure that the files downloaded in Step 1 and your input files be included in this directory.

-od  : Output file. This directory may or may not exist.

-m   : Database searching program you plan to use (diamond|usearch|blast).

-f   : Specify the extensions of your sequence files (E.g. fastq, fastq.gz, fasta, fasta.gz, fq, fq.gz, fa, fa.gz) or (faa, fna) or (diamond|usearch|blast).
     
      When using "-situation tablular", -f supports "diamond|usearch|blast".Ensure that filetype is support for the tool selected by -m option.
      
      (E.g., if -m usearch, the supported file types for -f are "fastq|fasta," and for "-m blast," they are "fasta|fa").
   
-s   : (nucl|prot) Sequence type.

-tpm : (0|1)  "1" need $sample.tpm exist in the work directory (default: 0)."-situation assembly-based" is a prerequisite for this option.

-id  : Minimum identity to report an alignment (default: 30).

-e   : Maximum e-value to report alignments (default: 1e-5).

-norm : (0|1) 0: don`t need random sampling; 1: need random sampling.

-rs   : The number of sequences for random subsampling. (default: the lowest number of sequences).Note: "-norm 1" is a prerequisite for this parameter

-thread : Number of threads (default: 2
```

---
### Examples

**1.read-based**
```
$ perl GetFun_CCycdb.pl -situation read-based -wd ./ -m diamond -f fasta -s nucl -norm 0 -thread 10 -od ./output

$ perl GetFun_CCycdb.pl -situation read-based -wd ./ -m diamond -f fasta -s nucl -norm 1 -rs 10000000 -thread 10 -od ./output

```
<p> Output:</p>

<li>FunProfile_read-based_$method_random.txt  OR  FunProfile_read-based_$method_norandom.txt:</li>

```
Gene    Mean identity   SampleA    SampleB
geneA         70           5         20
geneB         80           10        12
```

<li>SEQ2GENE/$sample.SEQ2G.txt :</li>

```
Query sequence                   Gene
k141_433371_length_91162_1      geneA
k141_455489_length_11328_1      geneB
```

**2.Assembly-based**

```
$ perl GetFun_CCycdb.pl -situation assembly-based -wd ./ -m diamond -f fatsa -s nucl -norm 0 -thread 10 -od ./output

$ perl GetFun_CCycdb.pl -situation assembly-based -wd ./ -m diamond -f fatsa -s nucl -tpm 1 -norm 0 -thread 10 -od ./output
```

<p>Output:</p>

- FunProfile_read-based_$method_random.txt OR FunProfile_read-based_$method_norandom.txt

- ORF2GENE/$sample.ORF2GENE.txt

- ORF2GENE.tpm (If "-tpm =1" and exist "$sample.tpm")


**3.Tabular Files**

```
$ perl GetFun_CCycdb.pl -situation tabular -wd ./ -m diamond -f diamond  -norm 0 -thread 10 -od ./output

$ perl GetFun_CCycdb.pl -situation tabular -wd ./ -m diamond -f diamond -norm 1 -thread 10 -od ./output
```

<p>Output:</b>

- FunProfile_read-based_$method_random.txt OR FunProfile_read-based_$method_norandom.txt

- SEQ2GENE/$sample.SEQ2GENE.txt
<br>

 
<p><b>Depending on the tools used, you may want to cite also:</b></p>

DIAMOND: Buchfink B, Xie C, Huson D H. Fast and sensitive protein alignment using DIAMOND[J]. Nature methods, 2015, 12(1): 59-60.

BLASTX: Boratyn G M, Camacho C, Cooper P S, et al. BLAST: a more efficient report with usability improvements[J]. Nucleic acids research, 2013, 41(W1): W29-W33.

USEARCH: Edgar R C. Search and clustering orders of magnitude faster than BLAST[J]. Bioinformatics, 2010, 26(19): 2460-2461.

CSVTK: Csvtk—CSV/TSV Toolkit. Available online: https://bioinf.shenwei.me/csvtk/

SEQKIT: Shen W, Le S, Li Y, et al. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation[J]. PloS one, 2016, 11(10): e0163962.
