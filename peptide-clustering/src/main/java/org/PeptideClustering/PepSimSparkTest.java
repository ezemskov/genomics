package org.PeptideClustering;

import org.PSSMHC.Impl;
import org.PSSMHC.Xml;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.Encoders;
import org.apache.spark.sql.SQLContext;

public class PepSimSparkTest
{
    public static void main(String[] args) throws Exception 
    {
        try
        {
            SparkConf conf = new SparkConf()
                    .setAppName("PeptideClustringSparkTest");
            JavaSparkContext jsc = new JavaSparkContext(conf);
            SQLContext sqlc = new SQLContext(jsc);

            Xml.Cfg pssmhcCfg = new Xml.Cfg(Xml.Utils.firstOrDef(args));
            XmlCfg appCfg = new XmlCfg(Xml.Utils.firstOrDef(args));
            
            JavaRDD<String> pepts = 
                sqlc.range(pssmhcCfg.start, pssmhcCfg.end, 1, appCfg.partitions)
                .map(new Impl.PeptideGenSparkFunc(), Encoders.STRING())
                .toJavaRDD();
            
            pepts.cache();
            System.out.format("Generated %d peptides\n", pepts.count());
            
            JavaPairRDD<String, String> pairs = pepts.cartesian(pepts);
            pairs.cache();
            System.out.format("Generated %d pairs\n", pairs.count());
            PeptideSimilarity sim = new PeptideSimilarity();
            sim.SetMatrix(new SubstMatrices.Blosum62());
            pairs.foreach(tuple -> { sim.posDiff(tuple._1, tuple._2); });
            System.out.format("Calculated %d diffs\n", pairs.count());
            
            jsc.close();
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
    
}
