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
            int peptideLengthInFile = -1;
            if (!cfg.peptideFiles.isEmpty())
            {
                Xml.StringIntPair peptFile = cfg.peptideFiles.get(0);
                peptideLengthInFile = peptFile.second;
                pepts = jsc.textFile(peptFile.first, cfg.partitions);
                peptQnty = pepts.count();
                if (!pepts.isEmpty() && (pepts.first().length() != peptideLengthInFile))
                {
                    throw new Exception(String.format(
                        "Invalid peptide length : %d in file, %d in xml config", 
                        pepts.first().length(), peptideLengthInFile
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
                        
            if (!cfg.doScore)
            {
                System.out.format("Generated %d peptides\n", peptQnty);
                return;
            }
            
            Impl.PSSMHCpanSparkFunc   PSSMHCpanFunc = null;
            Impl.ScoreFilterSparkFunc filterFunc = new Impl.ScoreFilterSparkFunc(cfg.ic50Threshold);
            for (Xml.PSSMCfg pssmConfig : cfg.pssmConfigs)
            {
                pssmConfig.peptideLength = peptideLengthInFile; //or in generator
                
                try
                {
                    PSSMHCpanFunc = new Impl.PSSMHCpanSparkFunc(pssmConfig);
                }
                catch(Exception ex)
                {
                    System.err.println("Skipped : " + ex.getMessage());
                    continue;
                }
                    
                JavaRDD<Impl.ScoredPeptide> binders = pepts.map(PSSMHCpanFunc)
                                                           .filter(filterFunc);
                if (cfg.doBinderPersist)
                {
                    binders.persist(StorageLevel.MEMORY_AND_DISK_SER());
                }

                if (cfg.doBinderStore)
                {
                    final String resDirName = String.format("output-%s-%d", 
                        pssmConfig.allele, pssmConfig.peptideLength);
                    binders.saveAsTextFile(resDirName, GzipCodec.class);
                }

                long binderCount = cfg.doBinderCount ? binders.count() : -1;
                System.out.format("%s/%d : found %d binder peptides (IC50 < %d) in total of %d\n", 
                    pssmConfig.allele, pssmConfig.peptideLength, 
                    binderCount, cfg.ic50Threshold, peptQnty);
            }
        }
        catch (Exception ex)
        {
            System.err.println("Error : " + ex.getMessage());
        }
    }
}
;