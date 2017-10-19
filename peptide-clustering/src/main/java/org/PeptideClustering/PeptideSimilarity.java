package org.PeptideClustering;

import info.debatty.spark.kmedoids.Similarity;

import java.io.Serializable;
import java.util.HashMap;

class Consts
{
    public static String alphabet = "ACDEFGHIKLMNPQRSTVWY";
}

class SubstMatrixRow extends HashMap<Character, Double>         implements Serializable {}
class SubstMatrix    extends HashMap<Character, SubstMatrixRow> implements Serializable
{
    public SubstMatrix(double[][] vals)
    {
        String ExcMsg =  "Wrong substitution matrix size";
        if (vals.length != Consts.alphabet.length()) 
        {
            throw new RuntimeException(ExcMsg);
        }

        for (int i=0; i<vals.length; ++i) 
        {
            if (vals[i].length != Consts.alphabet.length()) 
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

    public double get(Character aa1, Character aa2)
    {
        return get(aa1).get(aa2);
    }
}

public class PeptideSimilarity implements Similarity<String> 
{
    private static String alphabetRegex = "[" + Consts.alphabet + "]+";
    private static double[] PosWeights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    private static SubstMatrix SM = new SubstMatrices.Blosum62();
        
    public double similarity(String p1, String p2)
    {
        assert(p1.matches(alphabetRegex) && 
               p2.matches(alphabetRegex) &&
               (p1.length() == p2.length()) &&
               (PosWeights.length >= p1.length()));
        
        double sc12 = 0.0, sc11 = 0.0, sc22 = 0.0;
        for (int i=0; i<p1.length(); ++i)
        {
            char ch1 = p1.charAt(i);
            char ch2 = p2.charAt(i);
            sc12 += PosWeights[i] * SM.get(ch1, ch2);
            sc11 += PosWeights[i] * SM.get(ch1, ch1);
            sc22 += PosWeights[i] * SM.get(ch2, ch2);
        }
            
        return 2*sc12/(sc11 + sc22);
    }
}
