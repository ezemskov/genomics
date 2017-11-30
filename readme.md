# Peptide toolkit
This is a toolkit for bioinformatical calculations with peptides (short proteins/amino acid sequences) on Apache Spark.

## Features
 - Estimation of binding affinities (IC50 value) between peptides and human major histocompatibility complex (MHC) genes
 - Generation of all possible peptides with given length
 - Inter-peptide distance calculation using PAM/BLOSUM substitution matrices 
 - Clustering of peptides using k-medoids algorithm
 - Clustering of (generated) peptide set around given cluster centers set (of real-world peptides)

## Dependencies
 - **Apache Spark 2.2.0**  
    Available on <https://spark.apache.org/downloads.html>
 - **JUnit 4.12**  
    Available on <https://github.com/junit-team/junit4/>
 - **Apache Maven** buildsystem  
    Available on <http://maven.apache.org/>
 - **spark-kmedoids 0.1.2** : Apache Spark extension implementing k-medoids algorithm.  
    Available on <https://github.com/tdebatty/spark-kmedoids>

## 3rd-party code included or used
 - **PSSMHCpan-1.0**  : a toolkit for estimation of peptide binding affinities.  
It includes an amount of pre-calculated weight matrices (one for each HLA allele to peptide length pair) and a Perl script for binding affinity estimation.  
Original Perl code from this toolkit was rewritten in Java for Spark.  
One sample weight matrix (for HLA-A0201 allele and 9-meer peptides) is included into this package, the rest should be copied into the working tree from the original PSSMHCpan package as needed.  
Available on <https://github.com/BGI2016/PSSMHCpan>  
 - **NW-align** : Java implementation of Needleman-Wunsch global alignment, included (slightly modified) for benchmarking and comparison purposes only.  
Available on <http://zhanglab.ccmb.med.umich.edu/NW-align/>
 - **BLOSUM** and **PAM** amino acid substitution matrices :  reference matrices from NCBI are hardcoded.  
Available on <ftp://ftp.ncbi.nih.gov/blast/matrices/>

## Algorithms
In clustering applications, similarity `Sab` between peptides `A` and `B`, using substitution matrix `M`, is calculated as following :  
```
Sab = 2*SCab/(SCaa + SCbb), where  
SCxy = sum(M(x[i],y[i])), where  
M(k,n) is a substitution matrix value for amino acids k and n, where  
x[i] is the amino acid in peptide x in position i  
```

## Environment
 - Netbeans 8.2 IDE
 - Java 1.8.0
 - executed under Windows 7 x64 and various Linux distributions (Ubuntu 14, CentOs 6 and AltLinux 7.0.5)

## Running
 - Start Apache Spark master and worker(s) processes
 - Submit Spark task : 
>spark-submit --class <org.package.MainClassName> --master spark://master_hostname:7077 <package_filename.jar> [Config.xml]  

  * For peptide generation, class name is `org.PSSMHC.PSSMHCSpark` and package filename is `PSSMHCpan-1.0.jar`
  * For clustering around given centers, class name is `org.PeptideClustering.AssignBindersToClusters` and package filename is `peptide-clustering-1.0.jar`
  * For k-medoids clustering, class name is `org.PeptideClustering.PeptideClusteringMain` and package filename is `peptide-clustering-1.0.jar`

Default config filename is `PeptideCfg.xml` in current directory.  
Configuration is explained in xml comments to the sample config file `\PSSMHCpan\src\main\resources\PeptideCfg.xml`

So, example command line is : 
>spark-submit --class org.PSSMHC.PSSMHCSpark --master spark://master_hostname:7077 PSSMHCpan-1.0.jar

## Output : 
Running each of Spark applications will produce 
 - some logging in `stdout`
 - depending on xml configuration, some datasets (Spark RDDs saved as text) in `output-*` folders in Spark work directory.  
In Spark standalone mode, `stdout` is written the command line windows where you've executed `spark-submit`, and Spark working directory is your current directory.  
    
    * Peptide generation produces the folder `output-pssmhc` containing a gzip-ed set of binder peptides, formatted as :  
    `<peptide>, <binding affinity>`
    * Clustering around given centers produces folders 
    - `output-centers` containing a set of binder cluster centers (format described above) 
    - `output-clusters` containing the cluster elements, formatted as :  
    `(<center peptide>,[<element peptide 1>,<binding affinity 1>, ..., <element peptide N>,<binding affinity N>])`
    * K-medoids clustering only logs the solution into stdout in following format : 
```
Medoid : <center peptide> totalSim <sum of similarities between center and elements> avgSim <average similarity between center and elements>  
       Element : <element peptide 1> <binding affinity 1>  
       ...  
       Element : <element peptide N> <binding affinity N>  
```
