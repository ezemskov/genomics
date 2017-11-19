package org.PeptideClustering;

import info.debatty.spark.kmedoids.Similarity;
import org.PSSMHC.Impl;

import java.io.Serializable;
import java.util.HashMap;
import scala.Tuple2;

class PeptideSimilarity implements Similarity<String> 
{
    private SubstMatrix SM = null;
    private AminoPairSet pairs = null;
    
    public <T extends PeptideSimilarity> T SetMatrix(SubstMatrix SM_)
    {
        this.SM = SM_;
        this.pairs = new AminoPairSet(this.SM);
        //System.out.format("pairs %s\n", this.pairs.toString());
        return (T)this;
    }

    private void CheckParams(String p1, String p2)
    {
        assert(p1.matches(Impl.Consts.alphabetRegex));
        assert(p2.matches(Impl.Consts.alphabetRegex));
        assert(p1.length() == p2.length());
        assert(SM != null);
    }
    
    @Override
    public double similarity(String p1, String p2)
    {        
        CheckParams(p1, p2);
        double sc12 = 0.0, sc11 = 0.0, sc22 = 0.0;
        for (int i=0; i<p1.length(); ++i)
        {
            char ch1 = p1.charAt(i);
            char ch2 = p2.charAt(i);
            sc12 += SM.get(ch1).get(ch2);
            sc11 += SM.get(ch1).get(ch1);
            sc22 += SM.get(ch2).get(ch2);
        }
            
        return 2*sc12/(sc11 + sc22);
    }
    
    public double dissimilarity(String p1, String p2, double resMax)
    {
        CheckParams(p1, p2);
        double sc12 = 0.0, sc11_22 = 0.0;

        for (int i=0; i<p1.length(); ++i)
        {
            sc12 += SM.get(p1.charAt(i)).get(p2.charAt(i));
        }
        sc12 *= 2;
        
        double res = 0.0;
        double sumMax = resMax * sc12;
        for (int i=0; i<p1.length(); ++i)
        {
            char ch1 = p1.charAt(i);
            char ch2 = p2.charAt(i);
            sc11_22 += (SM.get(ch1).get(ch1) + SM.get(ch2).get(ch2));            
            if (sc11_22 >= sumMax) 
            { 
                return Double.NaN; 
            }
        }

        return sc11_22 / sc12;
    }

    public int posDiff(String p1, String p2)
    {
        assert(p1.length() == p2.length());
        int res = 0;
        for (int i=0; i<p1.length(); ++i)
        {
            Character c1 = p1.charAt(i);
            Character c2 = p2.charAt(i);
            if (!c1.equals(c2) &&
                !pairs.contains(new Tuple2<>(c1, c2)))
            {
                res += 1;
            }
        }
        return res;
    }
    
    //Maximum (empirical) amount of different amino acids present in 9-meer peptide pair
    //so it has similarity over similarityBound
    //Todo : generalize for different matrices and peptide lengths ?
    public static int maxPosDiff(double similarityBound)
    {
        if (similarityBound >= 0.9) { return 3; }
        if (similarityBound >= 0.8) { return 4; }
        if (similarityBound >= 0.7) { return 5; }
        if (similarityBound >= 0.5) { return 6; }
        return 9;
    }
}
