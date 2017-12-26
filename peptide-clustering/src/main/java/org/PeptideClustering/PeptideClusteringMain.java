package org.PeptideClustering;

import info.debatty.spark.kmedoids.Clusterer;
import info.debatty.spark.kmedoids.Solution;
import info.debatty.spark.kmedoids.budget.TrialsBudget;
import info.debatty.spark.kmedoids.neighborgenerator.ClaransNeighborGenerator;
import org.PSSMHC.Impl;
import org.PSSMHC.PeptideGen;
import org.PSSMHC.Xml;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;

public class PeptideClusteringMain
{
    public static void main(String[] args) throws Exception 
    {
        //try
        {
            XmlCfg appCfg = new XmlCfg(Xml.Utils.firstOrDef(args));

            SparkConf conf = new SparkConf()
                    .setAppName("SparkPeptideClustering");
            JavaSparkContext jsc = new JavaSparkContext(conf);

            Xml.Cfg pssmhcCfg = new Xml.Cfg(Xml.Utils.firstOrDef(args));
            Xml.PSSMCfg pssmConfig = pssmhcCfg.getSinglePSSMCfg();
            pssmConfig.peptideLength = PeptideGen.pepLen; //NB! should be some extra parameter in xml

            Impl.PSSMHCpanSparkFunc   pssmhcSparkFunc = new Impl.PSSMHCpanSparkFunc(pssmConfig);
            Impl.ScoreFilterSparkFunc ic50FilterSparkFunc = new Impl.ScoreFilterSparkFunc(pssmhcCfg.ic50Threshold);
            
            JavaRDD<String> binders = 
                jsc.textFile(appCfg.peptidesFilename, appCfg.partitions)
                   .map(pssmhcSparkFunc)
                   .filter(ic50FilterSparkFunc)
                   .map(scp -> { return scp.peptide; });
            System.out.format("Binders with IC50<%d qnty %d\n", pssmhcCfg.ic50Threshold, binders.count());

            PeptideSimilarity simCalc = new PeptideSimilarity().SetMatrix(SubstMatrices.get(appCfg.matrix));
            Clusterer<String> clusterer = new Clusterer<>();
            clusterer.setK(appCfg.clustersQnty);
            clusterer.setSimilarity(simCalc);
            clusterer.setNeighborGenerator(new ClaransNeighborGenerator<>());
            clusterer.setBudget(new TrialsBudget(appCfg.maxTrials));
            Solution<String> res = clusterer.cluster(binders);

            SolutionClusters<String> res2 = new SolutionClusters<>();
            res2.setSimilarity(simCalc);
            res2.AssignElemsToClusters(binders.collect(), res.getMedoids());

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
        /*
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
        */
    }
}
