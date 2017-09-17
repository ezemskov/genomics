package org.PSSMHC;

import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SparkSession;

import java.util.ArrayList;
import java.util.List;

class PSSMHCpanSpark extends PSSMHCpan
                     implements Function<String,ScoredPeptide>
{
    public ScoredPeptide call(String peptide)
    {
        return new ScoredPeptide(peptide, ScoreOnePeptide(peptide));
    }
}

public final class PSSMHCSparkTest
{   
    public static void main(String[] args) throws Exception 
    {
        try
        {
            SparkSession spark = SparkSession
                .builder()
                .appName("PSSMHCSparkTest")
                .getOrCreate();

            JavaSparkContext jsc = new JavaSparkContext(spark.sparkContext());
            PSSMHCpanSpark app = new PSSMHCpanSpark();
            app.InitFromCmdline(args);

            int partitions = 2;
            JavaRDD<String> peptideRDD = jsc.parallelize(app.peptides2, partitions);
            JavaRDD<ScoredPeptide> peptideRDD2 = peptideRDD.map(app);

            peptideRDD2.saveAsTextFile("OutputPSSMHC");

            spark.stop();
        }
        catch (Exception ex)
        {
            System.err.print(ex.getMessage());
        }
    }
}
