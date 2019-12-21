![](./GUIgraphics/logo.png)

## Installation
The recommended way to setup your computer for running **Inimotif** is to install the official distribution of Python 3.7. You can download the official Python 3.7 distribution [here](https://www.python.org/downloads/release/python-375). 

Once you've installed Python, download the zip file of **IniMotif** from this repository (click **Clone or download** on top right corner). Decompress this zip file and open the terminal/prompt and change to this directory. Enter the following command:

```bash
pip install -r requirements.txt
```

This will install the required python packages that Inimotif needs to run. Note that Inimotif may work for versions of Python 3 prior to 3.7, but this is untested. For Mac you may need to update Tcl/Tk to use the GUI for versions of Python prior to 3.7.2. See [here](https://www.python.org/download/mac/tcltk/) for more information.

Alternatively, you can also run **Inimotif** using [Anaconda](https://www.anaconda.com/distribution/). To setup a conda environment and install the requirements:

```bash
conda create -y --name inimotif python==3.7
conda install -fyq --name inimotif -c conda-forge --file requirements.txt
conda activate inimotif
```

## Run IniMotif
Open the GUI by opening the directory containing the python files in a terminal, and run the inimotif_gui.py using python.

```bash
$ python inimotif_gui.py
```

This will start the GUI of **IniMotif**. Select from one of the three options: ChP-seq, SELEX-seq, or Masker. Check the following tutorials for detailed function explanations.

## Tutorials
We will go through all functionalities of **IniMtoif** in the following list of tutorials.

* [Terms explanation](./tutorial/term_def.md)
	* kmer
	* Reverse complement
	* Palindrome
	* Motif
	* Motif logo
* [ChIP-seq Analysis](./tutorial/chipseq_tutorial.md) 
	* Motif discovery from ChIP-seq data
	* Query kmers
	* Scan motifs
* [SELEX-seq Analysis](./tutorial/selexseq_tutorial.md)
	* Motif discovery of SELEX-seq data
	* Query kmers
	* Scan motifs
* [Sequence masking](./tutorial/masker_tutorial.md)
	* Mask repetitive patterns
	* Mask motifs
* [Downloading datasets from European Nucleotide Archive](./tutorial/Downloading-sequence-data-from-ENA.md)
	* Downloading all relevant datasets by "accession number" using scripts
* [Extracting ChIP-seq windows from reference genome](./tutorial/Extracting-chip-seq-windows-from-hg19.md)
	* Extract sequences from given ChIP-seq windows on reference genome
