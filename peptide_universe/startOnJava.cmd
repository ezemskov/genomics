@echo %time%

set PSSMHC_JAR=C:\BigData\workspace\local_repo\genomics\peptide_universe\PSSMHCpan_maven\target\PSSMHCpan-1.0.jar
::@java org.PSSMHC.PSSMHCpanJava %1 9 HLA-A0201 database\PSSM\pssm_file.list > %2
::java -cp %PSSMHC_JAR% %1 9 HLA-A0201 database\PSSM\pssm_file.list > %2
java -cp %PSSMHC_JAR% org.PSSMHC.PeptideGen %1 %2

@echo %time%
