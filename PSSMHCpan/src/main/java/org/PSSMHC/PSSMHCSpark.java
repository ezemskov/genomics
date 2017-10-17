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

            Impl.PeptideGenSparkFunc gen = new Impl.PeptideGenSparkFunc();
            Impl.PSSMHCpanSparkFunc pssmhc = new Impl.PSSMHCpanSparkFunc();
            int nextArgIdx = pssmhc.InitFromCmdline(args);
            Impl.CmdlineCfg cfg = new Impl.CmdlineCfg(args, nextArgIdx);

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
                               .filter(new Impl.Ic50FilterFunc(cfg.ic50Threshold));
                                    
            if (cfg.doBinderPersist)
            {
                binderPepts.persist(StorageLevel.MEMORY_AND_DISK_SER());
            }
            
            if (cfg.doBinderStore)
            {            
                binderPepts.saveAsTextFile("output-pssmhc", GzipCodec.class);
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