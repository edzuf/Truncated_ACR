# Truncated_ACR
This repository gives access to the code used to perform data analysis on CRISPR arrays deep-sequenced using Illumina MiSeq.

All sequence alignments were performed using command-line BLAST. An example of a command used to perform an alignment between a CRISPR locus sequence and the repeat sequence is shown below:
```bash
blastn -query locus_sequence.fasta -subject Repeat.fasta -task blastn-short -ungapped -num_threads 32 -outfmt 6 -out locus_sequence.blastn
```
Similar BLAST commands were used to perform the alignments between the newly acquired spacers and the reference sequences.

The Python file `acquisition.py` counts of the number of new spacers acquired and then extracts them and sends them to a FASTA file containing all new spacers. The file `spacer_align.py` helps find the origin of each newly acquired spacer. 
`PAM_extraction.py` was then used to extract PAMs and see if they matched the expected motif.
