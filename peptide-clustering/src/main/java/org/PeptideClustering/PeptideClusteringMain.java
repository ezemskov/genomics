package org.PeptideClustering;

import info.debatty.spark.kmedoids.Clusterer;
import info.debatty.spark.kmedoids.Solution;
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
    static String CmdlineHelpStr = "Usage : java org.PeptideClustering.PeptideClusteringMain peptides_list.txt <partitions>\n";

    int partitions = 2;
    String peptidesFilename = "";
    
    public int InitFromCmdline(String[] args)
    {
        if (args.length < 2)
        {
            throw new RuntimeException(CmdlineHelpStr);
        }
        
        peptidesFilename = args[0];
        partitions = Integer.parseInt(args[1]);
        return 2;
    }
    
    public static void main(String[] args) throws Exception 
    {
        try
        {
            PeptideClusteringMain appCfg = new PeptideClusteringMain();
            appCfg.InitFromCmdline(args);

            SparkConf conf = new SparkConf()
                    .setAppName("Spark peptide clusterization");
            JavaSparkContext jsc = new JavaSparkContext(conf);
            JavaRDD<String> pepts = jsc.textFile(appCfg.peptidesFilename, appCfg.partitions);

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
            for (SolutionClusters.Cluster<String> cluster : res2.getClusters())
            {
                System.out.format("\nMedoid : %s\n", cluster.medoid);
                for (SolutionClusters.ElemSim<String> elem : cluster.elems)
                {
                    System.out.format("\tElement : %s %f.0\n", elem.elem, elem.sim);
                }
            }
            jsc.close();
        }
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
    }
}
