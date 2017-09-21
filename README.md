# pipeCoverage
Quickly calculate coverage and breadth of sequences from a mapping. Usefull to avoid creating large .sam files when only coverage inforamtion is desired.

## Quick start
`pipe_coverage.awk` is used to convert .sam formatted data into .pcov formatted data. .pcov files contain coverage information for each scaffold, as well as the coverage of windows accross the scaffolds (to get an idea of breadth). 

`pcov_tools.py` is used to process .pcov files into more useful formats. These includes coverage tables needed for automated-binned, and the files needed to make ESOMs.

## pipe_coverage.awk

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

## pcov_tools.py

This python script is used to process .pcov files into more useful outputs. For example:

**Generate coverage table from pcovs mapping to the same original sequnce:**
```
$ pcov_tools.py -p foo.fa-vs-sample01.pcov foo.fa-vs-sample02.pcov foo.fa-vs-sample03.pcov -o foo.fa.coverage -c --header

$ head foo.fa.coverage.cov
# scaffold      foo.fa-vs-sample01.pcov       foo.fa-vs-sample02.pcov       foo.fa-vs-sample03.pcov
scaffold_1516_read_length_100_read_count_18682        0.245003        0.0614925       0.0
scaffold_1573_read_length_100_read_count_174  0.0     0.0     0.0
scaffold_1432_read_length_100_read_count_331  0.0847442       0.0     2.75576
scaffold_2632_read_length_100_read_count_76   5.40076 60.4315 221.066
scaffold_1923_read_length_100_read_count_119  0.0     0.0     0.0
scaffold_592_read_length_100_read_count_427   0.0     0.0977909       0.0
scaffold_2407_read_length_100_read_count_54   0.0     0.0     0.0
scaffold_2522_read_length_100_read_count_1478 0.0     0.0     0.0
scaffold_1377_read_length_100_read_count_195  0.0     3.85477 0.24264
```

**See help for full list of commands:**
```$ pcov_tools.py -h
usage: pcov_tools.py [-h] -p [BCOVS [BCOVS ...]] [-o OUT] [-a] [-c] [-n] [-l]
                     [-auto] [--min_window MIN_WINDOW] [--header]
                     [--fix_names]

INPUT/OUTPUT:
  -h                    show this help message and exit
  -p [BCOVS [BCOVS ...]], --bcovs [BCOVS [BCOVS ...]]
                        .bcov file(s) (from pipe_coverage.awk) (default: None)
  -o OUT, --out OUT     output basename (default: None)

OPPERATIONS:
  -a, --all             generate all possible outputs (default: False)
  -c, --coverage        generate coverage file of complete scaffolds (CONCOCT
                        format) (default: False)
  -n, --names           generate esom.names files (default: False)
  -l, --learn           generate single UN-NORMALIZED esom.lrn file (default:
                        False)
  -auto, --auto         group .pcov filenames based on the first 12 characters
                        (default: False)

OTHER:
  --min_window MIN_WINDOW
                        minimum window size to allow for .esom files (default:
                        3000)
  --header              write a header in coverage table (will not work
                        natively with CONCOCT) (default: False)
  --fix_names           attempt to fix scaffold names messed up by SNAP
                        (default: False)
```

**API for other methods:**

There are some methods in pcov_tools.py that are useful for use in Jupyter notebooks. I need to document these better, but for example:

```
gen_genome_coverage_table(pcovs, stb, min_c = 1)
"""
Calcuate the coverage and breadth of genomes from pcov files, return a pandas dataframe with all this info

pcovs = list of .pcov files
stb = file containing scaffold to bin information
min_c = minimum coverage a window much have to be considered "present" during breadth calculattion
"""
```

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
