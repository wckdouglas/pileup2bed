# pileup2bed #

This is a program to parse *samtools mpileup* results to a tab-deliminated file storing base count at each position, this also allows one to choose a quality cut off and coverage cut off to report

to install


```
pip install git+https://github.com/wckdouglas/pileup2bed.git
```

usage:

```
$/pileup_to_bed.py -h
```

```
usage: pileup_to_bed.py [-h] -i INPUT [-c COV] [-q QUAL]

This program tries to convert samtools mipleup format to base count table with
coverage threshold and quality threshold

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        mpileup format input filename (default: <->)
  -c COV, --cov COV     coverage threshold (default: 10)
  -q QUAL, --qual QUAL  quality threshold (default: 33)
```
