package org.PeptideClustering;

import java.io.Serializable;
import java.util.Map;
import org.PSSMHC.Impl;
import org.PSSMHC.ScoredPeptide;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.sql.Encoders;
import org.apache.spark.sql.SQLContext;
import scala.Tuple2;
import scala.Tuple3;

public class AssignBindersToClusters
{
    static class SimFunc
                    implements Serializable,
                    Function<Tuple2<String, String>, Tuple3<String, String, Double>>
    {
        PeptideSimilarity simCalc = new PeptideSimilarity();

        @Override
        public Tuple3<String, String,Double> call(Tuple2<String, String> pepPair)
        {
            return new Tuple3(pepPair._1, pepPair._2,
                simCalc.similarity(pepPair._1, pepPair._2)
            );
        }
    }

    static class CmdlineCfg
    {
        static String CmdlineHelpStr = "Usage : java org.PeptideClustering.AssignBindersToClusters cluster_centers_list.txt partitions\n";

        int partitions;
        String peptidesFilename = "";

        public CmdlineCfg(String[] args, int firstArgIdx)
        {
            if (args.length < firstArgIdx+2)
            {
                throw new RuntimeException(CmdlineHelpStr);
            }

            peptidesFilename = args[firstArgIdx];
            partitions = Integer.parseInt(args[firstArgIdx + 1]);
        }
    }
    
    static <T> void FormatMap(Map<String, T> map)
    {
        for (Map.Entry<String, T> elem : map.entrySet())
        {
            System.out.format("%s %s\n", elem.getKey(), elem.getValue().toString());
        }
    }
    
    public static void main(String[] args) throws Exception 
    {
        try
        {

            SparkConf conf = new SparkConf()
                    .setAppName("Spark peptide clusterization : cluster generated binders around real centers");
            JavaSparkContext jsc = new JavaSparkContext(conf);
            SQLContext sqlc = new SQLContext(jsc);

            org.PSSMHC.Impl.ScoreFunc pssmhc = new org.PSSMHC.Impl.ScoreFunc();
            int nextArgIdx = pssmhc.InitFromCmdline(args);
            org.PSSMHC.Impl.CmdlineCfg pssmhcCfg = new org.PSSMHC.Impl.CmdlineCfg(args, nextArgIdx);
            CmdlineCfg appCfg = new CmdlineCfg(args, pssmhcCfg.NextArgIdx());

            PeptideGenFunc gen = new PeptideGenFunc();            
            org.PSSMHC.Impl.Ic50FilterFunc filterFunc = new org.PSSMHC.Impl.Ic50FilterFunc(pssmhcCfg.ic50Threshold);

            //todo : avoid storing scores !
            JavaRDD<String> binders = 
                sqlc.range(pssmhcCfg.start, pssmhcCfg.end, 1, appCfg.partitions)
                .map(gen, Encoders.STRING())
                .toJavaRDD()
                .map(pssmhc)
                .filter(filterFunc)
                .map(scp -> {return scp.peptide;});
            
            System.out.format("Binders qnty %d first %s\n", binders.count(), binders.first());
            JavaRDD<String> clusterCenters = jsc.textFile(appCfg.peptidesFilename, appCfg.partitions);

            System.out.format("Centers qnty %d distinct %d first %s\n", 
                clusterCenters.count(), clusterCenters.distinct().count(), clusterCenters.first());

            JavaPairRDD<String, String> pairs = clusterCenters.cartesian(binders);
            System.out.format("Pairs qnty %d first %s\n", pairs.count(), pairs.first().toString());

            SimFunc sim = new SimFunc();
            JavaRDD<Tuple3<String, String, Double>> simTriples = 
                pairs.map(sim)
                     .filter(triple -> { return triple._3() >= 0.2; });
            System.out.format("Filtered triples qnty %d first %s\n", simTriples.count(), simTriples.first().toString());
            
            // {Bn -> {(Ck, Snk)}}
            JavaPairRDD<String, ScoredPeptide> simPairs = simTriples.mapToPair(triple -> {
                ScoredPeptide centerWithSim = new ScoredPeptide(triple._1(), triple._3());
                return new Tuple2(triple._2(), centerWithSim); 
            });
            System.out.format("{Bn -> {(Ck, Snk)}} : qnty %d first %s\n", simPairs.count(), simPairs.first().toString());
            
            // {Bn -> (Ci_max, Sni_max)}
            JavaPairRDD<String, ScoredPeptide> simMaxPairs = simPairs.reduceByKey(
                (scp1, scp2) -> { return (scp1.ic50 > scp2.ic50) ? scp1 : scp2; } 
            );
            System.out.format("Binter to its cluster : qnty %d first %s\n", simMaxPairs.count(), simMaxPairs.take(2).toString());

            // {Ci_max -> (Bn, Sni_max)}
            JavaPairRDD<String, ScoredPeptide> simMaxPairsInv = simMaxPairs.mapToPair(tuple -> {
                ScoredPeptide binderWithSim = new ScoredPeptide(tuple._1, tuple._2.ic50);
                return new Tuple2<>(tuple._2.peptide, binderWithSim);
            });
            System.out.format("Clusters to binders : qnty %d first %s\n", simMaxPairsInv.count(), simMaxPairsInv.take(2).toString());
            
            // {Ck -> {Bi, Ski)}
            JavaPairRDD<String, Iterable<ScoredPeptide>> simMaxPairsGrp = simMaxPairsInv.groupByKey();
            simMaxPairsGrp.cache();

            System.out.format("Clusters to set of binders : qnty %d first %s\n", simMaxPairsGrp.count(), simMaxPairsGrp.take(2).toString());
            
            Map<String, Long> clustersCount = simMaxPairsGrp.countByKey();

            System.out.println("Cluster sizes");
            FormatMap(clustersCount);
            
            Map<String, Double> clusterSimCount = simMaxPairsGrp.mapValues(scpIt -> {
                double sumSim = 0.0;
                for (ScoredPeptide scp : scpIt) { sumSim += scp.ic50; }
                return sumSim;
            }).collectAsMap();

            System.out.println("Cluster similarities : ");
            FormatMap(clusterSimCount);
            
            jsc.close();
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
}

