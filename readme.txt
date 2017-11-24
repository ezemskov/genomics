This is a toolkit for bioinformatical calculations with peptides (short proteins/amino acid sequences) on Apache Spark.

Features : 
 - Estimation of binding affinities between peptides and human major histocompatibility complex (MHC) ??genes??
 - Generation of all possible peptides with given length
 - Inter-peptide distance calculation using PAM/BLOSUM substitution matrices 
 - Clusterization of peptides using k-medoids algorithm
 - Clusterization of (generated) peptide set around given cluster centers set (of real-world peptides)

Dependencies :
 - spark-kmedoids : Apache Spark extension implementing k-medoids algorithm. 
Available on https://github.com/tdebatty/spark-kmedoids.git , version used is 0.1.2 of 2017.09.24.    
 
- PSSMHCpan-1.0  : a toolkit for estimation of peptide binding affinities. 
Available on https://github.com/BGI2016/PSSMHCpan
It includes an amount of pre-calculated weight matrices (one for each HLA allele to peptide length pair) and a Perl script for binding affinity estimation. 
Original Perl code from this toolkit was rewritten in Java for Spark. One sample weight matrix (for HLA-A0201 allele and 9-meer peptides) is included into this package, the rest should be copied into the working tree from the original PSSMHCpan package as needed.
