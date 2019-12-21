# Term Definition

This tutorial explains what is motif defined by **IniMotif**.

## Kmer

kmer is a DNA sequence of length **k**. Considering all combination possibilities, there are in total 4<sup>k</sup> different kmers, given kmer length **k**.

## Reverse complement
Since there are two strands on the DNA, there are two DNA sequences for a DNA fragment. One on the forward strand and the other on the reverse strand. The DNA sequence on the reverse trand is called **reverse complement** sequence. 

For example, the reverse complement of kmer "AACTG" is "CAGTT".

## Palindrome
If the reverse complement of a kmer is the same as the kmer, this kmer is called a palindrome. For example, 5'-GGATCC-3' is a palindrome.

## Dimer
If the start and end parts of a kmer are reverse complement of each other, then this kmer is called a dimer. For example, 5'-TTGGCXXXXXGCCAA-3' is a dimer, where "X" can be any DNA base.

## Consensus sequence

In a sequence file such as ChIP-seq or SELEX-seq, we assume there exist a kmer that is bound by transcription factor (TF) with the highest affinity. As a consequence, this kmer appears with the highest frequency in the input sequence file.
We call the kmer with the highest count in the input file as **consensus sequence**.

## Motif
The transcription factor binds to the consensus sequence with the highest affinity. It will also bind to other kmers that are a few mutations to the consensus sequence.

Based on this observation, we define the motif as a set of kmers that are within a few mutations to the consensus sequence.

For example, if the consensus sequence is "**AACTG**" and we allow maximumly **1** mutation to the consensus sequence for the kmer to be bind by TF. Then the motif is given by the following list of kmers.

```python
AACTG (consensus) 
CACTG
GACTG
TACTG
ACCTG
.....
AACTT
```
The reverse complements of the kmers in the above list are also included in the motif by default. 

## Motif logo
Based on the counts of the kmers scanned from the input DNA sequence file, we could derive a position frequence matrix (PFM) or position weight matrix (PWM). Based on these matrixes, we could generate a logo, where the total height at each position is the information content and the heights of the letters are proportional to their frequencies.