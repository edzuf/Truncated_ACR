# Truncated_ACR
This repository gives access to the code used to perform data analysis on CRISPR arrays previously deep-sequenced using Illumina MiSeq.

All sequence alignments were performed using command-line BLAST. An example of the command used to perform an alignment between a CRISPR locus sequence and the repeat sequence is shown below:
```bash
blastn -query locus_sequence.fasta -subject References.fasta -task blastn-short -ungapped -num_threads 32 -outfmt 6 -out locus_sequence.blastn
```

The Python file `acquisition.py` counts of the number of new spacers acquired, while the file `spacer_align.py` helps find the origin of each new spacer.
