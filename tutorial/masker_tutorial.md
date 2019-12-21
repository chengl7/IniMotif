### Masker
The masker section will allow repeats or motifs within a file to be converted to N bases, and therefore not be read by IniMotif.

Masker will ask for an input file, and the output file name with location. Clicking the "Add Pattern" button will produce a new input row where the masks can be specified. Inputs for a pattern require: Type (repeat or motif), Reverse compliment (mask revcom of sequence), Sequence (to be masked), and depending on the type Number of minimum repeats, or Number of maximum mutations.

![Masker form entry](https://github.com/kearseya/IniMotif-py/blob/master/tutorial/screenshots/MaskerexampleGUI.png "Masker")

In this screenshot, repeats of 4 or more A's or T's and the barcode (TCTA/TAGA) are being masked from the NF1-1.fa file.
