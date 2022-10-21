+ the sample names is extracted from the BAM metadata
+ reads having a MAPQ>=30 and surounding the interval `chr7:27_135_310-27_135_339` are considered (not partial overlap).
+ for each read
    + the DNA Sequence of the read is extracted
    + the cigar string describing matches/deletions/insertions is extracted
    + the cigar string is scanned and , in the interval `chr7:27_135_310-27_135_339` we build the interval sequence.
    + the DNA is reverse complemented (minus strand) and translated
    + the peptide is scanned an its 'symbol' is extracted
+ each DNA sequence is counted
+ DNA sequences with occurence <=2 are discarded.
+ The two most frequent Sequences are displayed.