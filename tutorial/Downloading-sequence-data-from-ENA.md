# Downloading-sequence-data-from-ENA
enaFastqFetch.py can be used to query the ENA for different types of data and bulk download the associated fastq files.

CAUTION: Please be aware you may be downloading very large datasets. Before downloading the fastqs, the program will print to the terminal the total size of the files to be downloaded.


## **Usage**
```
usage: enaFastqFetch.py [-h] -s SEARCH -d DATATYPE [-n NUMRUNS]

optional arguments:
-h, --help            show this help message and exit
-s SEARCH, --search-term SEARCH
                      term you wish to search for, e.g. Mycobacterium, 1763,
                      SRR5188398, SRX2504319, PRJNA360902, SELEX
-d DATATYPE, --data-type DATATYPE
                      datatype you wish to search for, e.g. run, study,
                      experiment, sample
-n NUMDOWN, --num-downloads NUMDOWN
                      number of runs/studies/experiments/samples you wish to
                      download
-r, --report-file     generate a report file
```
## **How enaFastqFetch.py works**
enaFastqFetch works by querying the ENA's API for the corresponding XML file, then text mining the XML for the fastq download links. By providing enaFastqFetch with a search term and datatype it will automatically download the fastqs associated with your request.

### **Search term**
The search term is specified through the -s flag. You can search for data by using its accession ID, taxon ID, or by using free text searches.

### **Datatype**
The datatype is specified through the -d flag. The following values are valid: run, study, experiment, sample. The ENA uses a hierarchical system to define the datatype:

* STUDY:
A study defines an overarching investigation. In most cases it's a dataset associated with a publication.
* SAMPLE:
A sample is a biological sample which is used in a study.
* EXPERIMENT:
An experiment is conducted on a sample. This defines things like the instrument used for sequencing and the library preparation.
* RUN:
A run is the actual sequencing reads which are associated with a sample and experiment. I.e. these are the fastq files.

So for example the run [SRR3206414](https://www.ebi.ac.uk/ena/data/view/SRR3206414) is associated with the sample [SRS1318643](https://www.ebi.ac.uk/ena/data/view/SRS1318643) and the experiment [SRX1615315](https://www.ebi.ac.uk/ena/data/view/SRX1615315), and overall it belongs to the study [SRP071047/PRJNA313382](https://www.ebi.ac.uk/ena/data/view/PRJNA313382).

### **Number of runs/studies/experiments/samples to download**
The number of runs/studies/experiments/samples to download can be specified using the -n flag. This flag is optional. If not specified, all the data which matches the requested search term will be downloaded.

### **Report file**
A report file including information on the downloaded data will be generated if you pass the --report-file flag.

## **Examples of using the accession ID to download**
E.g. to download the fastq associated with the run SRR5188398:
```
python enaFastqFetch.py -s SRR5188398 -d run
```
E.g. to download all of the fastqs associated with the study PRJNA360902:
```
python enaFastqFetch.py -s PRJNA360902 -d study
```
E.g. to download all of the fastqs associated with the experiment SRX2504319:
```
python enaFastqFetch.py -s SRX2504319 -d experiment
```
E.g. to download all of the fastqs associated with the sample SAMN06240265:
```
python enaFastqFetch.py -s SAMN06240265 -d sample
```

## **Examples of using the taxon ID to download**
E.g. to download all the runs found for the taxon 47839:
```
python enaFastqFetch.py -s 47839 -d run
```
 E.g. to download 10 runs for the taxon 1773 and generate a report file:
```
python enaFastqFetch.py -s 1773 -d run -n 10 --report-file
```
## **Examples of using free text search to download**
E.g. to download 5 selex studies:
```
python enaFastqFetch.py -s "SELEX" -d study -n 5
```
