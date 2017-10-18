package org.PeptideClustering;

import java.io.Serializable;
import java.util.Map;
import java.util.List;
import org.PSSMHC.Impl;
import org.PSSMHC.Xml;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.Encoders;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.storage.StorageLevel;
import org.w3c.dom.Element;
import scala.Tuple2;
import scala.Tuple3;

public class AssignBindersToClusters
{
    static class PepSimSparkFunc
                    extends PeptideSimilarity
                    implements Serializable,
                    Function<Tuple2<String, String>, Tuple3<String, String, Double>>
    {
        @Override
        public Tuple3<String, String,Double> call(Tuple2<String, String> pepPair)
        {
            return new Tuple3(pepPair._1, pepPair._2,
                similarity(pepPair._1, pepPair._2)
            );
        }
    }

    static class XmlCfg
    {
        String peptidesFilename = "";
        int partitions = 0;

        public XmlCfg(String xmlFilename) throws Exception
        {
            Element root = Xml.Utils.parseXml(xmlFilename);

            partitions       = Integer.parseInt(Xml.Utils.getChildAttr(root, "spark", "partitions"));
            peptidesFilename = Xml.Utils.getChildAttr(root, "clustering", "peptidesFilePath");            
        }
    }
    
    static <T> void FormatMap(Map<String, T> map)
    {
        for (Map.Entry<String, T> elem : map.entrySet())
        {
            System.out.format("%s %s\n", elem.getKey(), elem.getValue().toString());
        }
    }

    static <T> void FormatList(List<Tuple2<String, T>> list)
    {
        for (Tuple2<String, T> elem : list)
        {
            System.out.format("%s %s\n", elem._1, elem._2.toString());
        }
    }

    static void FormatList3(List<Tuple3<String, String, Double>> list)
    {
        for (Tuple3<String, String, Double> elem : list)
        {
            System.out.format("%s %s %.2f\n", elem._1(), elem._2(), elem._3());
        }
    }
    
    public static void main(String[] args) throws Exception 
    {
        try
        {
            SparkConf conf = new SparkConf()
                    .setAppName("PeptideClustring2");
            JavaSparkContext jsc = new JavaSparkContext(conf);
            SQLContext sqlc = new SQLContext(jsc);

            Xml.Cfg pssmhcCfg = new Xml.Cfg(Xml.Utils.firstOrDef(args));
            XmlCfg appCfg = new XmlCfg(Xml.Utils.firstOrDef(args));

            JavaRDD<String> binders = 
                sqlc.range(pssmhcCfg.start, pssmhcCfg.end, 1, appCfg.partitions)
                .map(new Impl.PeptideGenSparkFunc(), Encoders.STRING())
                .toJavaRDD()
                .map(new Impl.PSSMHCpanSparkFunc(Xml.Utils.firstOrDef(args)))
                .filter(new Impl.ScoreFilterSparkFunc(pssmhcCfg.ic50Threshold))
                .map(scp -> {return scp.peptide;});
            
            binders.persist(StorageLevel.MEMORY_AND_DISK());
            System.out.format("Binders qnty %d\n", binders.count());
            JavaRDD<String> clusterCenters = jsc.textFile(appCfg.peptidesFilename, appCfg.partitions);

            JavaPairRDD<String, String> pairs = clusterCenters.cartesian(binders);

            JavaRDD<Tuple3<String, String, Double>> simTriples = 
                pairs.map(new PepSimSparkFunc())
                     .filter(triple -> { return triple._3() >= 0.8; });
            simTriples.persist(StorageLevel.MEMORY_AND_DISK());
            System.out.format("Filtered triples qnty %d \n", simTriples.count());
            
            // {Bn -> {(Ck, Snk)}}
            JavaPairRDD<String, Impl.ScoredPeptide> simPairs = simTriples.mapToPair(triple -> {
                Impl.ScoredPeptide centerWithSim = new Impl.ScoredPeptide(triple._1(), triple._3());
                return new Tuple2(triple._2(), centerWithSim); 
            });
            
            JavaPairRDD<String, Impl.ScoredPeptide> simMaxPairs = simPairs.reduceByKey(
                (scp1, scp2) -> { return (scp1.score > scp2.score) ? scp1 : scp2; } 
            );

            JavaPairRDD<String, Impl.ScoredPeptide> simMaxPairsInv = simMaxPairs.mapToPair(tuple -> {
                Impl.ScoredPeptide binderWithSim = new Impl.ScoredPeptide(tuple._1, tuple._2.score);
                return new Tuple2<>(tuple._2.peptide, binderWithSim);
            });
            simMaxPairsInv = simMaxPairsInv.coalesce(appCfg.partitions);
            simMaxPairsInv.persist(StorageLevel.MEMORY_AND_DISK());
                        
            Map<String, Long> clustersCount = simMaxPairsInv.countByKey();
            System.out.format("Clusters qnty %d\n", clustersCount.size());
            FormatMap(clustersCount);

            JavaPairRDD<String, Iterable<Impl.ScoredPeptide>> simMaxPairsGrp = simMaxPairsInv.groupByKey();
            simMaxPairsGrp.saveAsTextFile("output-clusters");
            
            jsc.close();
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
}

