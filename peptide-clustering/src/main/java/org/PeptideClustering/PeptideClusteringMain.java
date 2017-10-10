package org.PeptideClustering;

import info.debatty.spark.kmedoids.Clusterer;
import info.debatty.spark.kmedoids.Solution;
import info.debatty.spark.kmedoids.SolutionClusters;
import info.debatty.spark.kmedoids.Cluster;
import info.debatty.spark.kmedoids.budget.TrialsBudget;
import info.debatty.spark.kmedoids.neighborgenerator.ClaransNeighborGenerator;
import java.io.Serializable;
import java.util.Iterator;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SQLContext;

import org.PSSMHC.PeptideGen;
import org.apache.spark.sql.Encoders;

class PeptideGenFunc extends PeptideGen
                     implements Serializable, 
                     MapFunction<Row, String>
{
    public String call(Row idx)
    {
        return Generate(idx.getLong(0));
    }
}

public class PeptideClusteringMain
{
    public static void main(String[] args) throws Exception 
    {
        SparkConf conf = new SparkConf()
                .setAppName("Spark k-medoids clusterer");
        JavaSparkContext jsc = new JavaSparkContext(conf);
        SQLContext sqlc = new SQLContext(jsc);

        PeptideGenFunc gen = new PeptideGenFunc();
        JavaRDD<String> pepts = sqlc.range(0, 100, 1, 2)
                .map(gen, Encoders.STRING())
                .toJavaRDD();

        Clusterer<String> clusterer = new Clusterer<>();
        clusterer.setK(5);
        clusterer.setSimilarity(new PeptideSimilarity());
        clusterer.setNeighborGenerator(new ClaransNeighborGenerator<>());
        clusterer.setBudget(new TrialsBudget(30));
        Solution<String> res = clusterer.cluster(pepts);

        SolutionClusters<String> res2 = new SolutionClusters<>();
        res2.setSimilarity(new PeptideSimilarity());
        res2.AssignElemsToClusters(pepts.collect(), res.getMedoids());
            
        System.out.println(res.toString());
        for (Cluster<String> cluster : res2.getClusters())
        {
            System.out.format("\nMedoid : %s\n", cluster.medoid);
            for (Cluster.ElemSim<String> elem : cluster.elems)
            {
                System.out.format("\tElement : %s %f.0\n", elem.elem, elem.sim);
            }
        }
        jsc.close();
    }
}
