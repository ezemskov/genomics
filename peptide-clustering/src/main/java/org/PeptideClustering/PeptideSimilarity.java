package org.PeptideClustering;

import info.debatty.spark.kmedoids.Similarity;

class Consts
{
    public static String alphabet = "ACDEFGHIKLMNPQRSTVWY";
}

public class PeptideSimilarity implements Similarity<String> 
{
    private static String alphabetRegex = "[" + Consts.alphabet + "]+";

    public double similarity(String p1, String p2)
    {
        if (!p1.matches(alphabetRegex) || 
            !p2.matches(alphabetRegex) || 
             p1.length() != p2.length())
        {
            return 0.0;
        }
        
        double res = 0.0;

        for (int i=0; i<p1.length(); ++i)
        {
            if (p1.charAt(i) == p2.charAt(i))
            {
                res += 1.0;
            }
        }
            
        return res;
    }
}
