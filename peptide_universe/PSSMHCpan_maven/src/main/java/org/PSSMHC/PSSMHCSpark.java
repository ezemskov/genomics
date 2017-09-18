package org.PSSMHC;

import java.io.Serializable;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
//import org.apache.spark.sql.SparkSession;
//import org.apache.spark.sql.SQLContext;
//import org.apache.spark.sql.Dataset;
//import org.apache.spark.sql.Row;
import org.apache.spark.sql.*;

import java.util.ArrayList;
import java.util.List;

class PSSMHCpanSparkFunc extends PSSMHCpan
                        implements Serializable,
                        Function<String,ScoredPeptide>
{
    public ScoredPeptide call(String peptide)
    {
        return new ScoredPeptide(peptide, ScoreOnePeptide(peptide));
    }
}

class PeptideGenFunc implements Serializable, 
                     Function<Long, String>
{
    public String call(Long idx)
    {
        return PeptideGen.Generate(idx.longValue());
    }
}

final public class PSSMHCSpark
{   
    public static void main(String[] args) throws Exception 
    {
        try
        {
            SparkSession spark = SparkSession
                .builder()
                .appName("PSSMHCSpark")
                .getOrCreate();

            JavaSparkContext jsc = new JavaSparkContext(spark.sparkContext());
            SQLContext sqlc = new SQLContext(jsc);
            PSSMHCpanSparkFunc pssmhc = new PSSMHCpanSparkFunc();
            pssmhc.InitFromCmdline(args);

            int partitions = 2;
            PeptideGenFunc gen = new PeptideGenFunc();
            //sqlc.range((long)1010);
            //sqlc.range(1000, 1010);
            //sqlc.range(1000, 1010).map(gen);
            
            
            //Dataset<String> peptideDataset = peptideDS.map(gen);
            //Dataset<ScoredPeptide> scPepDataset = peptideDataset.map(pssmhc);

            //peptideDataset.toJavaRDD().saveAsTextFile("OutputPSSMHC_Dataset");
            
            //JavaRDD<String> peptideRDD = jsc.parallelize(pssmhc.peptides, partitions);
            //JavaRDD<ScoredPeptide> peptideRDD2 = peptideRDD.map(pssmhc);
            //peptideRDD2.saveAsTextFile("OutputPSSMHC_RDD");

            spark.stop();
        }
        catch (Exception ex)
        {
            System.err.print(ex + "\n");
            ex.printStackTrace();
        }
    }
}
