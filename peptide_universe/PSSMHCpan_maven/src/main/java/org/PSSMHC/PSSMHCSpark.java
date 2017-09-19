package org.PSSMHC;

import java.io.Serializable;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.Encoders;

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
    public static class Range
    {
        public long start;
        public long end;
    }
    
    public static Range ParseCmdline(String[] args, int firstArgIdx)
    {
        if (args.length < firstArgIdx+2)
        {
            throw new RuntimeException("Usage : java org.PSSMHC.PSSMHCpanJava peptides_list.fa <peptide_length> <allele name> database/PSSM/pssm_file.list [peptide idx start] [peptide qnty[\n");
        }

        Range res = new Range();
        res.start =           Long.parseUnsignedLong(args[firstArgIdx]);
        res.end = res.start + Long.parseUnsignedLong(args[firstArgIdx+1]);            
        return res;
    }
    
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
            int nextArgIdx = pssmhc.InitFromCmdline(args);

            Range idx = ParseCmdline(args, nextArgIdx);
            
            PeptideGenFunc gen = new PeptideGenFunc();
            sqlc.range(idx.start, idx.end)
                    .map(gen, Encoders.STRING())
                    .toJavaRDD()
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
