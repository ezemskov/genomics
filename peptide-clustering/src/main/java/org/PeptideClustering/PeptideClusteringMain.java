package org.PeptideClustering;

import info.debatty.spark.kmedoids.Clusterer;
import info.debatty.spark.kmedoids.Solution;
import info.debatty.spark.kmedoids.budget.TrialsBudget;
import info.debatty.spark.kmedoids.neighborgenerator.ClaransNeighborGenerator;
import org.PSSMHC.Impl;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.w3c.dom.Element;

public class PeptideClusteringMain
{
    public static class XmlCfg
    {
        int partitions, clustersQnty, maxTrials;
        String peptidesFilename = "";
        
        public XmlCfg(String xmlFilename) throws Exception
        {
            Element root = Impl.XmlUtils.parseXml(xmlFilename);
            Element elem = Impl.XmlUtils.getChildElem(root, "clustering");
            if (elem == null) { return; }
            
            peptidesFilename = elem.getAttribute("peptidesFilePath");
            clustersQnty     = Integer.parseInt(elem.getAttribute("clustersQnty"));
            maxTrials        = Integer.parseInt(elem.getAttribute("maxTrials"));
            partitions       = Integer.parseInt(Impl.XmlUtils.getChildAttr(root, "spark", "partitions"));
        }
    }
    
    public static void main(String[] args) throws Exception 
    {
        try
        {
            XmlCfg appCfg = new XmlCfg(Impl.XmlUtils.firstOrDef(args));

            SparkConf conf = new SparkConf()
                    .setAppName("SparkPeptideClustering");
            JavaSparkContext jsc = new JavaSparkContext(conf);
            JavaRDD<String> pepts = jsc.textFile(appCfg.peptidesFilename, appCfg.partitions);

            Clusterer<String> clusterer = new Clusterer<>();
            clusterer.setK(appCfg.clustersQnty);
            clusterer.setSimilarity(new PeptideSimilarity());
            clusterer.setNeighborGenerator(new ClaransNeighborGenerator<>());
            clusterer.setBudget(new TrialsBudget(appCfg.maxTrials));
            Solution<String> res = clusterer.cluster(pepts);

            SolutionClusters<String> res2 = new SolutionClusters<>();
            res2.setSimilarity(new PeptideSimilarity());
            res2.AssignElemsToClusters(pepts.collect(), res.getMedoids());

            System.out.println(res.toString());
            for (SolutionClusters.Cluster<String> cluster : res2.getClusters())
            {
                System.out.format("\nMedoid : %s totalSim %.2f avgSim %.2f\n", 
                    cluster.medoid, cluster.totalSim, cluster.totalSim/cluster.elems.size());
                for (SolutionClusters.ElemSim<String> elem : cluster.elems)
                {
                    System.out.format("\tElement : %s %.2f\n", elem.elem, elem.sim);
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
