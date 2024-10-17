# KMP_parallel

 In string computation, the exact pattern matching problem is the problem of 
finding all the occurrences of a pattern (string) P, in a text (string) S, where usually P is much 
shorter than S.
The Knuth-Morris-Pratt algorithm uses this approach: it first of all builds an index on P and then 
uses it to scan S, applying simple rules to the index to decide how to shift the pattern.

It takes as input:
1) FASTA file containing the genome
2) FASTA file containing a set of short sequences

Uses KMP to identify in parallel all matches of the short sequences 
within the longer genome sequence, and stores the matches, one per line, in a single file in 
a .txt format

Default parameter:
1) num_cores = 4
2) output file name = patterns_matches.txt
Here is an example of usage: python .\KMP.py --DNA .\genome.fasta --patterns .\100_random.fasta


For further information: python .\KMP.py --help 
