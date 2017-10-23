package org.PeptideClustering;

import info.debatty.spark.kmedoids.Similarity;

import java.io.Serializable;
import java.util.HashMap;

class Consts
{
    public static final String alphabet = "ACDEFGHIKLMNPQRSTVWY";
    public static final int    aLen = alphabet.length();
    public static final String alphabetRegex = "[" + alphabet + "]+";
    
    public static final double[] PosWeights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
}
class SubstMatrixRow extends HashMap<Character, Double>         implements Serializable
{
    SubstMatrixRow() { super(Consts.aLen, 1.0f); }
}

class SubstMatrix    extends HashMap<Character, SubstMatrixRow> implements Serializable 
{
    public SubstMatrix(double[][] vals)
    {
        super(Consts.aLen, 1.0f);
        
        final String ExcMsg =  "Wrong substitution matrix size";
        if (vals.length != Consts.aLen) 
        {
            throw new RuntimeException(ExcMsg);
        }

        for (int i=0; i<vals.length; ++i) 
        {
            if (vals[i].length != Consts.aLen) 
            {
                throw new RuntimeException(ExcMsg);
            }
            
            SubstMatrixRow row = new SubstMatrixRow();
            for (int j=0; j<vals.length; ++j) 
            {
                row.put(Consts.alphabet.charAt(j), vals[i][j]);
            }
            
            put(Consts.alphabet.charAt(i), row);
        }
    }
}

class PeptideSimilarity implements Similarity<String> 
{
    private SubstMatrix SM = null;
    
    public <T extends PeptideSimilarity> T SetMatrix(SubstMatrix SM_)
    {
        this.SM = SM_;
        return (T)this;
    }

    private void CheckParams(String p1, String p2)
    {
        assert(p1.matches(Consts.alphabetRegex));
        assert(p2.matches(Consts.alphabetRegex));
        assert(p1.length() == p2.length());
        assert(Consts.PosWeights.length >= p1.length());
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
            sc12 += Consts.PosWeights[i] * SM.get(ch1).get(ch2);
            sc11 += Consts.PosWeights[i] * SM.get(ch1).get(ch1);
            sc22 += Consts.PosWeights[i] * SM.get(ch2).get(ch2);
        }
            
        return 2*sc12/(sc11 + sc22);
    }
    
    public double dissimilarity(String p1, String p2, double resMax)
    {
        CheckParams(p1, p2);
        double sc12 = 0.0, sc11_22 = 0.0;

        for (int i=0; i<p1.length(); ++i)
        {
            sc12 += Consts.PosWeights[i] * SM.get(p1.charAt(i)).get(p2.charAt(i));
        }
        sc12 *= 2;
        
        double res = 0.0;
        double sumMax = resMax * sc12;
        for (int i=0; i<p1.length(); ++i)
        {
            char ch1 = p1.charAt(i);
            char ch2 = p2.charAt(i);
            sc11_22 += Consts.PosWeights[i] * (SM.get(ch1).get(ch1) + SM.get(ch2).get(ch2));            
            if (sc11_22 >= sumMax) 
            { 
                return Double.NaN; 
            }
        }

        return sc11_22 / sc12;
    }

    public static int posDiff(String p1, String p2)
    {
        assert(p1.length() == p2.length());
        int res = 0;
        for (int i=0; i<p1.length(); ++i)
        {
            if (p1.charAt(i) != p2.charAt(i))
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
        if (similarityBound >= 0.9) { return 2; }
        if (similarityBound >= 0.8) { return 4; }
        if (similarityBound >= 0.5) { return 6; }
        return 9;
    }
}
