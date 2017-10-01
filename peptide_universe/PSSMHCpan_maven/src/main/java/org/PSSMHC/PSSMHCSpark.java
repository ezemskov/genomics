package org.PSSMHC;

import java.io.Serializable;
import java.util.HashMap;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.Encoders;
import org.apache.spark.storage.StorageLevel;
import org.apache.spark.SparkConf; 
import org.apache.hadoop.io.compress.GzipCodec;

import org.apache.spark.serializer.KryoRegistrator; 
import com.esotericsoftware.kryo.Kryo; 
import com.esotericsoftware.kryo.serializers.FieldSerializer; 


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

class PeptideGenFunc extends PeptideGen
                     implements Serializable, 
                     MapFunction<Row, String>
{
    public String call(Row idx)
    {
        return Generate(idx.getLong(0));
    }
}

class PeptideGenAndScoreFunc extends PSSMHCpan
                             implements Serializable, 
                             MapFunction<Row, ScoredPeptide>
{
    private PeptideGen pepGen = new PeptideGen();
    
    public ScoredPeptide call(Row idx)
    {
        String peptide = pepGen.Generate(idx.getLong(0));
        double ic50 = ScoreOnePeptide(peptide);
        return (ic50 < 1500.0) ? new ScoredPeptide(peptide, ic50) : null;
    }
}

class PSSMHCSparkCmdlineConfig
{
    public long start;
    public long end;
    public int partitions;
    public boolean doGenAndScore, doBinderPersist, doBinderStore, doBinderCount;
        
    public PSSMHCSparkCmdlineConfig(String[] args, int firstArgIdx)
    {
        if (args.length < firstArgIdx+7)
        {
            throw new RuntimeException(PSSMHCpan.CmdlineHelpStr);
        }

        start =       ParseLongWithSuffix(args[firstArgIdx]);
        end = start + ParseLongWithSuffix(args[firstArgIdx+1]);    
        partitions =  Integer.parseUnsignedInt(args[firstArgIdx+2]);

        doGenAndScore     = (Integer.parseInt(args[firstArgIdx+3]) == 1);
        doBinderPersist   = (Integer.parseInt(args[firstArgIdx+4]) == 1);
        doBinderStore     = (Integer.parseInt(args[firstArgIdx+5]) == 1);
        doBinderCount     = (Integer.parseInt(args[firstArgIdx+6]) == 1);
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
    public static class ScPepRegistrator implements KryoRegistrator 
    { 
        public void registerClasses(Kryo kryo) 
        { 
          kryo.register(ScoredPeptide.class, new FieldSerializer(kryo, ScoredPeptide.class)); 
        } 
    } 
    
    
    public static void main(String[] args) throws Exception 
    {
        try
        {
            SparkConf spconf = new SparkConf()
                .setAppName("PSSMHCSpark"); 
            
            spconf.set("spark.kryo.registrator", ScPepRegistrator.class.getName()); 
            JavaSparkContext jsc = new JavaSparkContext(spconf);
            SQLContext sqlc = new SQLContext(jsc);
            
            PeptideGenAndScoreFunc genPssmhc = new PeptideGenAndScoreFunc();
            int nextArgIdx = genPssmhc.InitFromCmdline(args);
            PSSMHCSparkCmdlineConfig cfg = new PSSMHCSparkCmdlineConfig(args, nextArgIdx);


            PeptideGenFunc gen = new PeptideGenFunc();
            PSSMHCpanSparkFunc pssmhc = new PSSMHCpanSparkFunc();
            pssmhc.InitFromCmdline(args);

            
            JavaRDD<ScoredPeptide> binderPepts;
            if (cfg.doGenAndScore)
            { 
                binderPepts = sqlc.range(cfg.start, cfg.end, 1, cfg.partitions)
                        .map(genPssmhc, Encoders.kryo(ScoredPeptide.class))
                        .toJavaRDD()
                        .filter(scPep -> (scPep != null));
            }
            else
            {
                binderPepts = sqlc.range(cfg.start, cfg.end, 1, cfg.partitions)
                        .map(gen, Encoders.STRING())
                        .toJavaRDD()
                        .map(pssmhc)
                        .filter(scPep -> (scPep.ic50 < 1500.0));
            }
                                    
            if (cfg.doBinderPersist)
            {
                binderPepts.persist(StorageLevel.MEMORY_AND_DISK_SER());
            }
            
            if (cfg.doBinderStore)
            {            
                binderPepts.saveAsTextFile("OutputPSSMHC", GzipCodec.class);
            }
            
            long binderCount = cfg.doBinderCount ? binderPepts.count() : -1;
            System.out.format("Found %d binder peptides in total of %d\n", binderCount, (cfg.end-cfg.start));
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
}
;