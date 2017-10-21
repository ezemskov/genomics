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

    @Override
    public double similarity(String p1, String p2)
    {
        assert(p1.matches(Consts.alphabetRegex) && 
               p2.matches(Consts.alphabetRegex) &&
               (p1.length() == p2.length()) &&
               (Consts.PosWeights.length >= p1.length()) && 
               (SM != null));
        
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
}
