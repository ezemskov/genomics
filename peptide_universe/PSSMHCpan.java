package org.PSSMHC;

import java.lang.Math;
import java.util.*;
import java.lang.System;
import java.lang.String;
import java.io.*;

class PSSMHCParser
{
}

class PSSMHCpan
{
    static String alphabet = "ACDEFGHIKLMNPQRSTVWY";
    static String peptide = "CLMPCGRRQ";

    static double score_max = 0.8;
    static double score_min = 0.8 * (1 - Math.log(50000) / Math.log(500));
    static double score_range = score_max - score_min;

    public static class WeightMatrixRow extends HashMap<Character, Double> {}
    public static class WeightMatrix extends ArrayList<WeightMatrixRow> {}

    public static class AllelePair { public String name; public int length; }
    public static class WeightMatrices extends HashMap<AllelePair, WeightMatrix> {}


    WeightMatrix weight_matrix = new WeightMatrix();

//BufferedReader in = new BufferedReader(new FileReader("foo.in"));

    private static WeightMatrixRow RandomRow()
    {
        WeightMatrixRow res = new WeightMatrixRow();
        for (int i=0; i<alphabet.length(); i++)
        {
            res.put(new Character(alphabet.charAt(i)), 2 * Math.random());
        }
        return res;
    }

    PSSMHCpan()
    {
        for (int i=0; i<peptide.length(); i++)
        {
            this.weight_matrix.add(RandomRow());
        }
    }

    double score_one_peptide()
    {
        double score = 0;
        for (int ch_pos=0; ch_pos<peptide.length(); ch_pos++)
        {
            Map<Character, Double> weightRow = weight_matrix.get(ch_pos);
            Double weight = weightRow.get(peptide.charAt(ch_pos));
            score += weight;
        }

        score = score / peptide.length();
        score = Math.max(Math.min(score, score_max), score_min);

        return Math.pow(50000, (score_max - score)/score_range);
    }

    public static void main(String[] args) 
    {
        PSSMHCpan b = new PSSMHCpan();

        long startTime = System.currentTimeMillis();
        try
        {
            PrintWriter writer = new PrintWriter("res.txt", "UTF-8");

            for (int i=0; i<500000; i++)
            {
                double ic50 = b.score_one_peptide();
            }

            for (int i=0; i<50000; i++)
            {
                writer.format("%s %f\n", peptide, 5.5555);
            }
            writer.close();
        } 
        catch (IOException e) {
           // do something
        }        
        long duration = System.currentTimeMillis() - startTime;
        System.out.format("Duration %d msec\n", duration);

    }
}