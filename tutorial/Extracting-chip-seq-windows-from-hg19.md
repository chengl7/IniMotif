# **Extracting-chip-seq-windows-from-hg19**
chipWinExtract.py extracts chip-seq windows from reference chromosomes and builds a fasta file for each chromosome.

## **Usage**
```
usage: chipWinExtract.py [-h] -t INPUTTSV -r INPUTREF -o OUTPUTFASTA -i
COLNAME -s COLSTART -e COLEND

optional arguments:
    -h, --help            show this help message and exit
    -t INPUTTSV, --input-tsv INPUTTSV
                          Path to the input chip-seq tsv
    -r INPUTREF, --input-ref INPUTREF
                          Path to the directory where the reference chromosomes are
    -o OUTPUTFASTA, --output-dir OUTPUTFASTA
                          Output directory path for the generated fasta files
    -i COLNAME, --col-name COLNAME
                          Column number for the chromosome name
    -s COLSTART, --col-start COLSTART
                          Column number for the start coordinate of the window
    -e COLEND, --col-end COLEND
                          Column number for the end coordinate of the window
    -c, --concat-fasta    Concatenate the final fasta files into one file

```
chipWinExtract.py requires an input chip-seq tsv file, with columns for the name of each chromosome and the start and end window coordinates. You will also need the corresponding reference chromosomes to extract the sequence windows from.

It splits the input tsv by chromosome, then creates a fasta file for each chromosome's extracted windows. The name of the chromosome in the tsv file should be the same as the reference chromosome, which is assumed to have the file extension .fa* (e.g. .fa or .fasta).

## **Example for hg19**
An example chip-seq file VR_HB13_hg19.tsv is included in the directory chipseqwindows.

To create a fasta file from VR_HB13_hg19.tsv we will need the reference chromosomes for hg19. They can be downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/). Once the chromosomes have been downloaded, place them in a directory called chromosomes.

So to run chipWinExtract.py for this data:
```
python chipWinExtract.py -t chipseqwindows/VR_HB13_hg19.tsv -r chromosomes -o chipseqwindows -i 1 -s 2 -e 3 --concat-fasta
```
this will create a fasta file for each chromosome and a final concatenated fasta which contains all the chromosomes in the chipseqwindows directory.
