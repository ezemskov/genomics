package org.PSSMHC;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SQLContext;
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

            String xmlFilename = args[0];
            Impl.XmlCfg cfg = new Impl.XmlCfg(xmlFilename);

            JavaRDD<String> pepts = sqlc.range(cfg.start, cfg.end, 1, cfg.partitions)
                    .map(new Impl.PeptideGenSparkFunc(), Encoders.STRING())
                    .toJavaRDD();
            
            JavaRDD<ScoredPeptide> binderPepts;
            if (!cfg.doScore)
            {
                System.out.format("Generated %d peptides\n", pepts.count());
                return;
            }
            
            binderPepts = pepts.map(new Impl.PSSMHCpanSparkFunc(xmlFilename))
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