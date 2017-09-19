@rmdir /s /q OutputPSSMHC

spark-submit --class org.PSSMHC.PSSMHCSpark --master spark://192.168.56.1:7077 C:\BigData\workspace\local_repo\genomics\peptide_universe\PSSMHCpan_maven\target\PSSMHCpan-1.0.jar dummy.fa 9 HLA-A0201 database\PSSM\pssm_file.list %*
