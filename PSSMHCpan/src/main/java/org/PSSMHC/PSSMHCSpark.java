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

class PeptideIc50FilterFunc implements Serializable,
                        Function<ScoredPeptide,Boolean>
{
    public PeptideIc50FilterFunc(double ic50Threshold_)
    {
        ic50Threshold = ic50Threshold_;
    }
        
    public Boolean call(ScoredPeptide scPep)
    {
        return (scPep.ic50 < ic50Threshold);
    }

    double ic50Threshold;
}

class PSSMHCSparkCmdlineConfig
{
    public long start;
    public long end;
    public int partitions;
    public boolean doScore, doBinderPersist, doBinderStore, doBinderCount;
    public int ic50Threshold;
        
    public PSSMHCSparkCmdlineConfig(String[] args, int firstArgIdx)
    {
        if (args.length < firstArgIdx+7)
        {
            throw new RuntimeException(PSSMHCpan.CmdlineHelpStr);
        }

        start =       ParseLongWithSuffix(args[firstArgIdx]);
        end = start + ParseLongWithSuffix(args[firstArgIdx+1]);    
        partitions =  Integer.parseUnsignedInt(args[firstArgIdx+2]);

        doScore     = (Integer.parseInt(args[firstArgIdx+3]) == 1);
        doBinderPersist   = (Integer.parseInt(args[firstArgIdx+4]) == 1);
        doBinderStore     = (Integer.parseInt(args[firstArgIdx+5]) == 1);
        doBinderCount     = (Integer.parseInt(args[firstArgIdx+6]) == 1);
        if (args.length >= 8)
        {
            ic50Threshold = Integer.parseInt(args[firstArgIdx+7]);
        }
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

            PeptideGenFunc gen = new PeptideGenFunc();
            PSSMHCpanSparkFunc pssmhc = new PSSMHCpanSparkFunc();
            int nextArgIdx = pssmhc.InitFromCmdline(args);
            PSSMHCSparkCmdlineConfig cfg = new PSSMHCSparkCmdlineConfig(args, nextArgIdx);
            PeptideIc50FilterFunc filterFunc = new PeptideIc50FilterFunc(cfg.ic50Threshold);

            JavaRDD<String> pepts = sqlc.range(cfg.start, cfg.end, 1, cfg.partitions)
                    .map(gen, Encoders.STRING())
                    .toJavaRDD();
            
            JavaRDD<ScoredPeptide> binderPepts;
            if (!cfg.doScore)
            {
                System.out.format("Generated %d peptides\n", pepts.count());
                return;
            }
            
            binderPepts = pepts.map(pssmhc)
                               .filter(filterFunc);
                                    
            if (cfg.doBinderPersist)
            {
                binderPepts.persist(StorageLevel.MEMORY_AND_DISK_SER());
            }
            
            if (cfg.doBinderStore)
            {            
                binderPepts.saveAsTextFile("OutputPSSMHC", GzipCodec.class);
            }
            
            long binderCount = cfg.doBinderCount ? binderPepts.count() : -1;
            System.out.format("Found %d binder peptides (IC50 < %d) in total of %d\n", 
                binderCount, cfg.ic50Threshold, (cfg.end-cfg.start));
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
}
;