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
    Xml.Cfg _cfg;
    JavaRDD<String> _pepts;
    long _peptQnty;
    int _peptideLength;
    String _peptFilename;
    JavaSparkContext _jsc;
    
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
        _peptFilename = "";

        SparkConf spconf = new SparkConf().setAppName("PSSMHCSpark"); 
        spconf.set("spark.kryo.registrator", ScPepRegistrator.class.getName()); 
        _jsc = new JavaSparkContext(spconf);
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
                final String resDirName = String.format("%s-%s-%d", 
                    _peptFilename, pssmConfig.allele, pssmConfig.peptideLength);
                binders.saveAsTextFile(resDirName, GzipCodec.class);
            }

            long binderCount = _cfg.doBinderCount ? binders.count() : -1;
            System.out.format("%s/%d : found %d binder peptides (IC50 < %d) in total of %d\n", 
                pssmConfig.allele, pssmConfig.peptideLength, 
                binderCount, _cfg.ic50Threshold, _peptQnty);
        }
    }

    public void generatePeptideSet()
    {
        _peptFilename = "gen";

        Xml.PeptideGenCfg genCfg = _cfg.genCfg;
        SQLContext sqlc = new SQLContext(_jsc);
        _peptideLength = genCfg.peptideLength;
        _peptQnty = genCfg.end - genCfg.start;
        _pepts = sqlc.range(genCfg.start, genCfg.end, 1, _cfg.partitions)
                .map(new Impl.PeptideGenSparkFunc(_peptideLength), Encoders.STRING())
                .toJavaRDD();
    }

    public boolean readPeptideSetFromFile(Xml.StringIntPair peptFileCfg)
    {
        _peptideLength = peptFileCfg.second;
        _peptFilename = peptFileCfg.first;
        _pepts = _jsc.textFile(_peptFilename, _cfg.partitions);
        _peptQnty = _pepts.count();
        
        final boolean res = _pepts.isEmpty() || 
                           (_pepts.first().length() == _peptideLength);
        if (!res)
        {
            System.err.format("Ignore file %s with invalid peptide length %d, %d in xml config\n", 
                _peptFilename, _pepts.first().length(), _peptideLength);
        }
        return res;
    }
    
    public void run() throws Exception 
    {
        if (_cfg.peptideFiles.isEmpty())
        {
            generatePeptideSet();
            processPeptideSet();
            return;
        }

        for (Xml.StringIntPair peptFile : _cfg.peptideFiles)
        {
            if (readPeptideSetFromFile(peptFile))
            {
                processPeptideSet();
            }
        }
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