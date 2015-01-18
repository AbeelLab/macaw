# Macaw
Macaw is a tool that lineage-types MTB samples and identifies strain mixes in WGS data sets.

Quick-run recipe:

  `java -jar macaw.jar -o <output.txt> <input.bam>`
  
  `python macaw-utilities.py interpret -m <output.txt> -o <result>`

# Downloads and documentation

Versioned production releases: https://github.com/AbeelLab/macaw/releases

Nightly development builds: http://www.abeellab.org/jenkins/macaw/ 

Source-code: https://github.com/AbeelLab/macaw.git

Submit bugs and features requests: https://github.com/AbeelLab/macaw/issues

# Documentation

## Running marker detection

Usage: `java -jar macaw.jar [options] -o <output> <file>...`

Required arguments:
```
  -o <output> | --output <output>
        File where you want the output to be written
  <file>
        Input files (BAM)
```

Optional arguments:
```
  --marker <value>
        File containing marker sequences. This file has to be a multi-fasta file with the headers indicating the name of the markers.
  -t <value> | --threshold <value>
        Threshold to determine absence or presence of a marker (default=5)
```

## Interpreting detected markers

Usage: `python macaw-utilities.py interpret -m <marker detection file> -o <output>`

Required arguments:
```
  -m <marker detection file> | --macaw <marker detection file>
        File with absence/presence as determined by macaw.jar
  -o <output> | --output <output>
        Output file prefix (.macaw.interpreted automatically appended)
```
 
## Condensing results

Usage: `python macaw-utilities.py condense -d <input file directory> -o <output>`

Required arguments:
```
  -d <input file directory> | --directory <input file directory>
        Directory with all .interpreted.macaw files in it
  -o <output> | --output <output>
        Output file
```


# Requirements
- Python 2.7 + NumPy+ SCiPy 
- Java 1.7+
