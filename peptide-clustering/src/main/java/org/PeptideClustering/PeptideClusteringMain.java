package org.PeptideClustering;

import info.debatty.spark.kmedoids.Clusterer;
import info.debatty.spark.kmedoids.Solution;
import info.debatty.spark.kmedoids.budget.TrialsBudget;
import info.debatty.spark.kmedoids.neighborgenerator.ClaransNeighborGenerator;
import java.io.Serializable;
//import org.apache.log4j.Level;
//#import org.apache.log4j.Logger;
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
        //try
        {
            SparkConf conf = new SparkConf()
                    .setAppName("Spark k-medoids clusterer");
            JavaSparkContext jsc = new JavaSparkContext(conf);
            SQLContext sqlc = new SQLContext(jsc);

            PeptideGenFunc gen = new PeptideGenFunc();
            JavaRDD<String> pepts = sqlc.range(0, 1000, 1, 2)
                    .map(gen, Encoders.STRING())
                    .toJavaRDD();

            Clusterer<String> clusterer = new Clusterer<>();
            clusterer.setK(10);
            clusterer.setSimilarity(new PeptideSimilarity());
            clusterer.setNeighborGenerator(new ClaransNeighborGenerator<>());
            clusterer.setBudget(new TrialsBudget(30));
            Solution<String> res = clusterer.cluster(pepts);
            System.out.println(res.toString());
            for (String medoid : res.getMedoids())
            {
                System.out.format("Medoid : %s\n", medoid);
            }
            jsc.close();
        }
    }
}
