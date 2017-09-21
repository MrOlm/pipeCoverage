# pipeCoverage
Quickly calculate coverage and quasi-breadth of sequences from a mapping. Usefull to avoid creating large .sam files when only coverage inforamtion is desired.

## Description
`pipe_coverage.awk` is used to convert .sam formatted data into .pcov formatted data. .pcov files contain coverage information for each scaffold, as well as the coverage of windows accross the scaffolds (to get an idea of breadth). 

`pcov_tools.py` is used to process .pcov files into more useful formats. These includes coverage tables needed for automated-binned, and the files needed to make ESOMs.

## `pipe_coverage.awk`

This awk script can be used to .pcov files in a variety of ways. For example:

**Directly from mapping algorithm:**
```
Bowtie2:
$ bowtie2  -x bt2/foo.fa -1 r1.fq.gz -2 r2.fq.gz 2> foo-vs-sample01_mapped.log  | shrinksam | awk -f pipe_coverage.awk > foo.fa-vs-sample01.pcov

Snap:
$ snap-aligner paired foo.fa r1.fq.gz r2.fq.gz -t 6 -F a -o -sam - 2> foo.fa-vs-sample01.mapped.log | awk -f pipe_coverage.awk > foo.fa-vs-sample01.pcov
```

**From .sam file:**
```
$ cat foo.fasta-vs-sample01.sam | awk -f pipe_coverage.awk > foo.fasta-vs-sample01.pcov
```

**From .bam file:**
```
samtools view foo.fasta-vs-sample01.bam | awk -f pipe_coverage.awk > foo.fasta-vs-sample01.pcov
```

## `pcov_tools.py`

This python script is used to process .pcov files into more useful outputs.


### .pcov file format

.pcov files provide the number of reads and coverage of each scaffold *as well as* the reads and coverage of windows of a specified length accross each scaffold. For example:

```
$ head foo.fa-vs-sample01.pcov
# window length: 3000
# read length: 96.8757 (based on first 10000 reads)
# scaffold       coverage        length
b003-d044_scaffold_429_read_length_100_read_count_624   0       7542
>b003-d044_scaffold_429_read_length_100_read_count_624_1        0       3000
>b003-d044_scaffold_429_read_length_100_read_count_624_2        0       3000
>b003-d044_scaffold_429_read_length_100_read_count_624_3        0       1542
b003-d044_scaffold_1588_read_length_100_read_count_96   1.71969 1859
>b003-d044_scaffold_1588_read_length_100_read_count_96_1        1.71969 1859
```

The `window length` specifies how large each window should be. By default this is 3000

The `read length` is the average read length, calculated on-the-fly based on the first X reads in the file. By default this is 10000

Following this header is the coverage of all scaffolds and windows. Windows are differentiated from scaffolds in that they start with `>`, and end with `_x`.

Following the above example, we can see the first scaffold in this sample is nammed `b003-d044_scaffold_429_read_length_100_read_count_624`, it's coverage is 0, and it's 7542bp long.  In this example the window length is 3000, so the next 3 lines will specify the coverage of the 3 windows that make up this 7542bp sequence- from 1-3000, 3001-6000, and 6001-7542. The coverage and length of these windows are shown on the next 3 lines (in this case they also have 0 coverage).

### `pipe_coverage.awk` variables

There are a two variables that can be passed into `pipe_coverage.awk`:

**BL** specifies the window length (default = 3000). Smaller will lead to better estimations of breadth, but larger file sizes.

**RC** specifices the number of reads to use to calculate the average read-length (which is used to calculate coverage) (default = 10000). Higher will lead to more accurate RL/coverage calculations, but longer run-times.

The variables are changed by directly passing them into the awk script. For example:

```
$ cat foo.fasta-vs-sample01.sam | awk -f pipe_coverage.awk -v BL=1000 > foo.fasta-vs-sample01.pcov
```
