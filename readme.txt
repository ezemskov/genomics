This is a toolkit for bioinformatical calculations with peptides (short proteins/amino acid sequences) on Apache Spark.

Features : 
 - Estimation of binding affinities between peptides and human major histocompatibility complex (MHC) genes
 - Generation of all possible peptides with given length
 - Inter-peptide distance calculation using PAM/BLOSUM substitution matrices 
 - Clustering of peptides using k-medoids algorithm
 - Clustering of (generated) peptide set around given cluster centers set (of real-world peptides)

Dependencies :
 - Apache Spark : version spark-2.2.0-bin-hadoop2.7 was used.
    Available on https://spark.apache.org/downloads.html
 - JUnit : version 4.12 was used.
    Available on https://github.com/junit-team/junit4/
 - Apache Maven builsystem : bundled with NetBeans IDE (?)
    Available on http://maven.apache.org/
 - spark-kmedoids : Apache Spark extension implementing k-medoids algorithm.
    Available on https://github.com/tdebatty/spark-kmedoids.git , version 0.1.2 was used.
 
3rd-party code included :
- PSSMHCpan-1.0  : a toolkit for estimation of peptide binding affinities. 
    Available on https://github.com/BGI2016/PSSMHCpan
    It includes an amount of pre-calculated weight matrices (one for each HLA allele to peptide length pair) and a Perl script for binding affinity estimation. 
    Original Perl code from this toolkit was rewritten in Java for Spark. 
    One sample weight matrix (for HLA-A0201 allele and 9-meer peptides) is included into this package, the rest should be copied into the working tree from the original PSSMHCpan package as needed.
 - NW-align : Java implementation of Needleman-Wunsch global alignment, included (slightly modified) for benchmarking and comparison purposes only.
    Available on http://zhanglab.ccmb.med.umich.edu/NW-align/

Environment :
    Toolkit was built with Netbeans 8.2 IDE, using Java 1.8.0 and tested under Windows 7 x64 and various Linux distributions (Ubuntu 14, CentOs 6 and AltLinux 7.0.5)

Building : 
 - Checkout spark-kmedoids and peptide-clustering repositories
 - Buld spark-kmedoids
 - Add spark-mkedoids dependency ?
 - Import POM project into NetBeans and build it. Maven will download all the required dependencies and produce two JAR packages : \PSSMHCpan\target\PSSMHCpan-1.0.jar and \peptide-clustering\target\peptide-clustering-1.0.jar 

Running : 
 - Start Apache Spark master and worker(s) processes
 - Submit task : 
    >spark-submit --class <org.package.MainClassName> --master spark://master_hostname:7077 <package_filename.jar> [Config.xml]
   Default config filename is PeptideCfg.xml in current directory. Config file format is explained in xml comments to sample config file \PSSMHCpan\src\main\resources\PeptideCfg.xml

For peptide generation, class name is org.PSSMHC.PSSMHCSpark and package filename is PSSMHCpan-1.0.jar 
For k-medoids clustering, class name is org.PeptideClustering.PeptideClusteringMain and package filename is peptide-clustering-1.0.jar
For clustering around given centers, class name org.PeptideClustering.AssignBindersToClusters and package filename is peptide-clustering-1.0.jar

So, example command line is : 
    >spark-submit --class org.PSSMHC.PSSMHCSpark --master spark://master_hostname:7077 PSSMHCpan-1.0.jar

Inputs : 

Outputs : 