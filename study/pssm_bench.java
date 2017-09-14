import java.lang.Math;
import java.util.*;
import java.lang.System;
import java.lang.String;
import java.io.*;

class PSSMbench 
{
    static char alphabet[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    static String peptide = "CLMPCGRRQ";

    static double score_max = 0.8;
    static double score_min = 0.8 * (1 - Math.log(50000) / Math.log(500));
    static double score_range = score_max - score_min;

    List<Map<Character, Double>> weight_matrix = new ArrayList<Map<Character, Double>>();

    private static Map<Character, Double> RandomRow()
    {
        //Map<Character, Double> res = new HashMap<Character, Double>(alphabet.length);
        Map<Character, Double> res = new HashMap<Character, Double>();
        for (int i=0; i<alphabet.length; i++)
        {
            res.put(new Character(alphabet[i]), 2 * Math.random());
        }
        return res;
    }

    PSSMbench()
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
        PSSMbench b = new PSSMbench();

        long startTime = System.currentTimeMillis();
        try
        {
            PrintWriter writer = new PrintWriter("res.txt", "UTF-8");

            for (int i=0; i<1000000; i++)
            {
                double ic50 = b.score_one_peptide();
                writer.format("%s %f\n", peptide, ic50);
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