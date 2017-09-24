package org.PSSMHC;

import java.io.Serializable;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.Encoders;
import org.apache.spark.storage.StorageLevel;
import org.apache.hadoop.io.compress.GzipCodec;

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

class PeptideBloomFilterFunc extends PeptideBloomFilter
                        implements Serializable,
                        VoidFunction<ScoredPeptide>
{
    public void call(ScoredPeptide scPep)
    {
        Add(scPep.peptide);
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
        public int partitions;
    }
    
    public static Range ParseCmdline(String[] args, int firstArgIdx)
    {
        if (args.length < firstArgIdx+3)
        {
            throw new RuntimeException(PSSMHCpan.CmdlineHelpStr);
        }

        Range res = new Range();
        res.start =           Long.parseUnsignedLong(args[firstArgIdx]);
        res.end = res.start + Long.parseUnsignedLong(args[firstArgIdx+1]);    
        res.partitions =          Integer.parseUnsignedInt(args[firstArgIdx+2]);
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
            JavaRDD<ScoredPeptide> scPepts = sqlc.range(idx.start, idx.end, 1, idx.partitions)
                        .map(gen, Encoders.STRING())
                        .toJavaRDD()
                        .map(pssmhc);
            //scPepts.persist(StorageLevel.MEMORY_AND_DISK());
            
            JavaRDD<ScoredPeptide> binderPepts = scPepts.filter(scPep -> (scPep.ic50 < 1500.0));
            binderPepts.persist(StorageLevel.MEMORY_AND_DISK());
            binderPepts.saveAsTextFile("OutputPSSMHC", GzipCodec.class);
            
            System.out.format("Found %d binder peptides in total of %d\n", binderPepts.count(), (idx.end-idx.start));
            
            spark.stop();
        }
        catch (Exception ex)
        {
            System.err.print(ex.getMessage() + "\n");
        }
    }
}
;