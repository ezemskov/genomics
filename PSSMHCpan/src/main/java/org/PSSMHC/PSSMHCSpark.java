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
        @Override
        public void registerClasses(Kryo kryo) 
        { 
          kryo.register(Impl.ScoredPeptide.class, 
                new FieldSerializer(kryo, Impl.ScoredPeptide.class)); 
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

            Xml.Cfg cfg = new Xml.Cfg(Xml.Utils.firstOrDef(args));
            
            JavaRDD<String> pepts;
            long peptQnty = -1;
            if (!cfg.peptidesFilename.isEmpty())
            {
                pepts = jsc.textFile(cfg.peptidesFilename, cfg.partitions);
                peptQnty = pepts.count();
                if (!pepts.isEmpty() && (pepts.first().length() != cfg.peptideLength))
                {
                    throw new Exception(String.format(
                        "Invalid peptide length : %d in file, %d in xml config", 
                        pepts.first().length(), cfg.peptideLength
                    ));
                }
            }
            else
            {
                pepts = sqlc.range(cfg.start, cfg.end, 1, cfg.partitions)
                        .map(new Impl.PeptideGenSparkFunc(), Encoders.STRING())
                        .toJavaRDD();
                peptQnty = cfg.end - cfg.start;
            }
                        
            JavaRDD<Impl.ScoredPeptide> binderPepts;
            if (!cfg.doScore)
            {
                System.out.format("Generated %d peptides\n", peptQnty);
                return;
            }
            
            binderPepts = pepts.map(new Impl.PSSMHCpanSparkFunc(Xml.Utils.firstOrDef(args)))
                               .filter(new Impl.ScoreFilterSparkFunc(cfg.ic50Threshold));
                                    
            if (cfg.doBinderPersist)
            {
                binderPepts.persist(StorageLevel.MEMORY_AND_DISK_SER());
            }
            
            if (cfg.doBinderStore)
            {
                final String resDirName = String.format("output-%s-%d", cfg.alleleName, cfg.peptideLength);
                binderPepts.saveAsTextFile(resDirName, GzipCodec.class);
            }
            
            long binderCount = cfg.doBinderCount ? binderPepts.count() : -1;
            System.out.format("Found %d binder peptides (IC50 < %d) in total of %d\n", 
                binderCount, cfg.ic50Threshold, peptQnty);
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
}
;