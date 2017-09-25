package org.PSSMHC;

import java.io.Serializable;
import java.util.HashMap;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.Encoders;
import org.apache.spark.storage.StorageLevel;
import org.apache.hadoop.io.compress.GzipCodec;

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

class PSSMHCSparkCmdlineConfig
{
    public long start;
    public long end;
    public int partitions;
        
    public PSSMHCSparkCmdlineConfig(String[] args, int firstArgIdx)
    {
        if (args.length < firstArgIdx+3)
        {
            throw new RuntimeException(PSSMHCpan.CmdlineHelpStr);
        }

        start =       ParseLongWithSuffix(args[firstArgIdx]);
        end = start + ParseLongWithSuffix(args[firstArgIdx+1]);    
        partitions =  Integer.parseUnsignedInt(args[firstArgIdx+2]);
    }

    private static long ParseLongWithSuffix(String val)
    {
        int suffixPos = val.length()-1;
        Character suffix = val.toUpperCase().charAt(suffixPos);
        Long multiplier = suffixes.get(suffix);
        if ((val.length() < 2) || (multiplier == null))
        {
            return Long.parseUnsignedLong(val);
        }
        
        return multiplier * Long.parseUnsignedLong(val.substring(0, suffixPos));
    }

    private static final HashMap<Character, Long> suffixes;
    static
    {
        suffixes = new HashMap<>();
        suffixes.put('K', (long)1E3);
        suffixes.put('M', (long)1E6);
        suffixes.put('G', (long)1E9);
        suffixes.put('T', (long)1E12);
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
            int nextArgIdx = pssmhc.InitFromCmdline(args);

            PSSMHCSparkCmdlineConfig cfg = new PSSMHCSparkCmdlineConfig(args, nextArgIdx);
            
            PeptideGenFunc gen = new PeptideGenFunc();
            JavaRDD<ScoredPeptide> scPepts = sqlc.range(cfg.start, cfg.end, 1, cfg.partitions)
                        .map(gen, Encoders.STRING())
                        .toJavaRDD()
                        .map(pssmhc);
            //scPepts.persist(StorageLevel.MEMORY_AND_DISK());
            
            JavaRDD<ScoredPeptide> binderPepts = scPepts.filter(scPep -> (scPep.ic50 < 1500.0));
            binderPepts.persist(StorageLevel.MEMORY_AND_DISK());
            binderPepts.saveAsTextFile("OutputPSSMHC", GzipCodec.class);
            
            System.out.format("Found %d binder peptides in total of %d\n", binderPepts.count(), (cfg.end-cfg.start));
            
            spark.stop();
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
}
;