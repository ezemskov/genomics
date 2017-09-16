package org.PSSMHC;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SparkSession;

import java.util.ArrayList;
import java.util.List;

public final class PSSMHCSparkTest
{
    
    public static void main(String[] args) throws Exception 
    {
        SparkSession spark = SparkSession
          .builder()
          .appName("PSSMHCSparkTest")
          .getOrCreate();

        JavaSparkContext jsc = new JavaSparkContext(spark.sparkContext());
        
        PSSMHCPan app = new PSSMHCpan(args);        
        int partitions = 2;
        
        app.ScoreAllPeptides();
        /*
        JavaRDD<Integer> dataSet = jsc.parallelize(l, slices);

        int count = dataSet.map(integer -> {
          double x = Math.random() * 2 - 1;
          double y = Math.random() * 2 - 1;
          return (x * x + y * y <= 1) ? 1 : 0;
        }).reduce((integer, integer2) -> integer + integer2);

        System.out.println("Pi is roughly " + 4.0 * count / n);
    */
        spark.stop();
    }
}
