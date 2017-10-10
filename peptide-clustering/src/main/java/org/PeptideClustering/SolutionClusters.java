package org.PeptideClustering;

import info.debatty.spark.kmedoids.Similarity;
import java.util.ArrayList;
import java.util.List;

//Assigns elements to pre-calculated medoids
public class SolutionClusters<T>
{
    public static class ElemSim<T>
    {
        public ElemSim(T elem_, double sim_)
        {
            elem = elem_;
            sim = sim_;
        }

        public T elem;
        public double sim = Double.NaN;
    }
    
    public static class Cluster<T>
    {
        public T medoid;
        public List<ElemSim<T>> elems = new ArrayList<>();
        public double totalSim = 0.0;
    }
    
    public final void setSimilarity(final Similarity<T> simCalc_) 
    {
        this.simCalc = simCalc_;
    }
    
    public List<Cluster<T>> getClusters()
    {
        return clusters;
    }
    
    public void AssignElemsToClusters(List<T> elements, List<T> medoids)
    {
        if (simCalc == null) {
            throw new IllegalStateException("Similarity is not defined!");
        }

        clusters.clear();
        for (T med : medoids)
        {
            Cluster c = new Cluster();
            c.medoid = med;
            clusters.add(c);
        }
        
        for (T elem : elements)
        {
            double currSims[] = new double[medoids.size()];
            for (int i = 0; i < medoids.size(); ++i)
            {
                currSims[i] = simCalc.similarity(elem, medoids.get(i));
            }
            int cluster_index = argmax(currSims);
            if (elem.equals(medoids.get(cluster_index))) //don't add medoid itself to the cluster
            {
                continue;
            }
            
            Cluster<T> cluster = clusters.get(cluster_index);
            cluster.elems.add(new ElemSim<>(elem, currSims[cluster_index]));
            cluster.totalSim += currSims[cluster_index];
            clusters.set(cluster_index, cluster);
        }
    }        

    static int argmax(final double[] values) {
        double max = -Double.MAX_VALUE;
        int max_index = -1;
        for (int i = 0; i < values.length; i++) {
            if (values[i] > max) {
                max = values[i];
                max_index = i;
            }
        }

        return max_index;
    }
    
    private ArrayList<Cluster<T>> clusters = new ArrayList<>();
    private Similarity<T> simCalc = null;
}
