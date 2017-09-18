package org.PSSMHC;

import java.io.Serializable;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.Encoders;
import org.apache.spark.api.java.function.MapFunction;

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
                     MapFunction<Row, String>
{
    public String call(Row idx)
    {
        return PeptideGen.Generate(idx.getLong(0));
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

            int partitions = 3;
            PeptideGenFunc gen = new PeptideGenFunc();
            sqlc.range(1000, 1050, 1, partitions)
                    .map(gen, Encoders.STRING())
                    .toJavaRDD()
                    .map(pssmhc)
                    .saveAsTextFile("OutputPSSMHC_Dataset");
            
            jsc.parallelize(pssmhc.peptides, partitions)
                    .map(pssmhc)
                    .saveAsTextFile("OutputPSSMHC");

            spark.stop();
        }
        catch (Exception ex)
        {
            System.err.print(ex + "\n");
            ex.printStackTrace();
        }
    }
}
