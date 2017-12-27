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

final class PSSMHCSpark
{
    Xml.Cfg _cfg;
    JavaRDD<String> _pepts;
    long _peptQnty;
    int _peptideLength;
    
    public static class ScPepRegistrator implements KryoRegistrator 
    { 
        @Override
        public void registerClasses(Kryo kryo) 
        { 
          kryo.register(Impl.ScoredPeptide.class, 
                new FieldSerializer(kryo, Impl.ScoredPeptide.class)); 
        } 
    } 

    private PSSMHCSpark(Xml.Cfg cfg)
    {
        _cfg = cfg;
        _pepts = null;
        _peptQnty = -1;
        _peptideLength = -1;
    }
    
    private void processPeptideSet()
    {
        Impl.PSSMHCpanSparkFunc   PSSMHCpanFunc = null;
        Impl.ScoreFilterSparkFunc filterFunc = new Impl.ScoreFilterSparkFunc(_cfg.ic50Threshold);
        for (Xml.PSSMCfg pssmConfig : _cfg.pssmConfigs)
        {
            pssmConfig.peptideLength = _peptideLength;

            try
            {
                PSSMHCpanFunc = new Impl.PSSMHCpanSparkFunc(pssmConfig);
            }
            catch(Exception ex)
            {
                System.err.println("Skipped : " + ex.getMessage());
                continue;
            }

            JavaRDD<Impl.ScoredPeptide> binders = _pepts.map(PSSMHCpanFunc)
                                                       .filter(filterFunc);
            if (_cfg.doBinderPersist)
            {
                binders.persist(StorageLevel.MEMORY_AND_DISK_SER());
            }

            if (_cfg.doBinderStore)
            {
                final String resDirName = String.format("output-%s-%d", 
                    pssmConfig.allele, pssmConfig.peptideLength);
                binders.saveAsTextFile(resDirName, GzipCodec.class);
            }

            long binderCount = _cfg.doBinderCount ? binders.count() : -1;
            System.out.format("%s/%d : found %d binder peptides (IC50 < %d) in total of %d\n", 
                pssmConfig.allele, pssmConfig.peptideLength, 
                binderCount, _cfg.ic50Threshold, _peptQnty);
        }
    }

    public void run() throws Exception 
    {
        SparkConf spconf = new SparkConf()
            .setAppName("PSSMHCSpark"); 

        spconf.set("spark.kryo.registrator", ScPepRegistrator.class.getName()); 
        JavaSparkContext jsc = new JavaSparkContext(spconf);
        SQLContext sqlc = new SQLContext(jsc);

        if (!_cfg.peptideFiles.isEmpty())
        {
            Xml.StringIntPair peptFile = _cfg.peptideFiles.get(0);
            _peptideLength = peptFile.second;
            _pepts = jsc.textFile(peptFile.first, _cfg.partitions);
            _peptQnty = _pepts.count();
            if (!_pepts.isEmpty() && (_pepts.first().length() != _peptideLength))
            {
                throw new Exception(String.format(
                    "Invalid peptide length : %d in file, %d in xml config", 
                    _pepts.first().length(), _peptideLength
                ));
            }
        }
        else
        {
            _pepts = sqlc.range(_cfg.start, _cfg.end, 1, _cfg.partitions)
                    .map(new Impl.PeptideGenSparkFunc(), Encoders.STRING())
                    .toJavaRDD();
            _peptQnty = _cfg.end - _cfg.start;
            _peptideLength = Impl.PeptideGenSparkFunc.pepLen;
        }

        if (!_cfg.doScore)   //todo : remove if
        {
            System.out.format("Generated %d peptides\n", _peptQnty);
            return;
        }

        processPeptideSet();
    }
    
    public static void main(String[] args) throws Exception 
    {
        try
        {
            Xml.Cfg cfg = new Xml.Cfg(Xml.Utils.firstOrDef(args));
            new PSSMHCSpark(cfg).run();
        }
        catch (Exception ex)
        {
            System.err.println("Error : " + ex.getMessage());
        }
    }
}
;