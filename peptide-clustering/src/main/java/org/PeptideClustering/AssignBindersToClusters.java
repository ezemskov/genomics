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
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.sql.Encoders;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.storage.StorageLevel;
import scala.Tuple2;
import scala.Tuple3;

public class AssignBindersToClusters
{
    static class PepSimSparkFunc
                    extends PeptideSimilarity
                    implements Serializable,
                    PairFunction<Tuple2<String, String>, String, Impl.ScoredPeptide>
    {
        private final double dissimMax;
        
        PepSimSparkFunc(double dissimMax_)
        {
            dissimMax = dissimMax_;
        }
            
        @Override
        public Tuple2<String, Impl.ScoredPeptide> call(Tuple2<String, String> pepPair)
        {
            double score = dissimilarity(pepPair._1, pepPair._2, dissimMax);
            return new Tuple2<>(pepPair._1, new Impl.ScoredPeptide(pepPair._2, score));
        }
    }

    static class TupleScoreFilterSparkFunc 
                            implements Serializable,
                            Function<Tuple2<String, Impl.ScoredPeptide>, Boolean>
    {
        public TupleScoreFilterSparkFunc(double scoreThreshold_)
        {
            scoreThreshold = scoreThreshold_;
        }

        @Override
        public Boolean call(Tuple2<String, Impl.ScoredPeptide> tuple)
        {
            return !Double.isNaN(tuple._2.score) && 
                   (tuple._2.score < scoreThreshold);
        }

        double scoreThreshold;
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
            Impl.PSSMHCpanSparkFunc   pssmhcSparkFunc = new Impl.PSSMHCpanSparkFunc(Xml.Utils.firstOrDef(args));
            Impl.ScoreFilterSparkFunc ic50FilterSparkFunc = new Impl.ScoreFilterSparkFunc(pssmhcCfg.ic50Threshold);
            
            JavaRDD<String> binders = 
                sqlc.range(pssmhcCfg.start, pssmhcCfg.end, 1, appCfg.partitions)
                .map(new Impl.PeptideGenSparkFunc(), Encoders.STRING())
                .toJavaRDD()
                .map(pssmhcSparkFunc)
                .filter(ic50FilterSparkFunc)
                .map(scp -> {return scp.peptide;});
            
            binders.persist(StorageLevel.MEMORY_AND_DISK());
            System.out.format("Binders with IC50<%d qnty %d\n", pssmhcCfg.ic50Threshold, binders.count());
            
            JavaRDD<String> clusterCenters = 
                jsc.textFile(appCfg.peptidesFilename, appCfg.partitions)
                   .map(pssmhcSparkFunc)
                   .filter(ic50FilterSparkFunc)
                   .map(scp -> { return scp.peptide; });

            System.out.format("Cluster centers with IC50<%d qnty %d\n", pssmhcCfg.ic50Threshold, clusterCenters.count());

            JavaPairRDD<String, String> pairs = binders.cartesian(clusterCenters);

            PepSimSparkFunc simFunc = new PepSimSparkFunc(1/appCfg.minSimilarity);
            simFunc.SetMatrix(SubstMatrices.get(appCfg.matrix));
                
            // {Bn -> {(Ck, Snk)}}
            JavaPairRDD<String, Impl.ScoredPeptide> simPairs = 
                pairs.mapToPair(simFunc)
                     .filter(new TupleScoreFilterSparkFunc(1/appCfg.minSimilarity));
            simPairs.persist(StorageLevel.MEMORY_AND_DISK());
            System.out.format("Pairs with similarity>%.2f qnty %d\n", appCfg.minSimilarity, simPairs.count());
            
            // {Bn -> (Ck_max, Snk_max)}
            JavaPairRDD<String, Impl.ScoredPeptide> simMaxPairs = simPairs.reduceByKey(
                (scp1, scp2) -> { return (scp1.score < scp2.score) ? scp1 : scp2; } 
            );

            // {Ck -> (Bn, Snk)}
            JavaPairRDD<String, Impl.ScoredPeptide> simMaxPairsInv = simMaxPairs.mapToPair(tuple -> {
                Impl.ScoredPeptide binderWithSim = new Impl.ScoredPeptide(tuple._1, 1.0/tuple._2.score);
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

