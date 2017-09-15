//package PSSMHC;

import java.lang.Math;
import java.util.*;
import java.lang.System;
import java.lang.String;
import java.io.*;

class WeightMatrixColumn extends HashMap<Character, Double> {}
class WeightMatrix extends ArrayList<WeightMatrixColumn> {}

class AllelePair
{
    public int hashCode()
    {
        return String.format("%s_%d", alName, pepLength).hashCode();
    }

    public boolean equals(Object obj)
    {
        AllelePair ap = (AllelePair)obj;
        return (ap != null) && (alName.equals(ap.alName)) && (pepLength == ap.pepLength);
    }

    public String alName; 
    public int pepLength; 
}
class WeightMatrices extends HashMap<AllelePair, WeightMatrix> {}

class PSSMParser
{
    public static String alphabet = "ACDEFGHIKLMNPQRSTVWY";

    //filePath is pssm_file.list
    public static WeightMatrices ParsePSSMfileList(String filePath)
    {
        WeightMatrices res = new WeightMatrices();
        try
        {
            BufferedReader reader = new BufferedReader(new FileReader(filePath));
            while (reader.ready())
            {
                String rowStr = reader.readLine();
                String[] row = rowStr.split("\t");
                if (row.length != 3)
                {
                    System.err.format("Skip invalid line PSSM list line '%s'\n", rowStr);
                    continue;
                }

                AllelePair al = new AllelePair();
                al.alName = row[1];
                al.pepLength = Integer.parseInt(row[2]);   

                String pssmFilePath = row[0];
                WeightMatrix pssm  = ParseOnePSSM(pssmFilePath);
                res.put(al, pssm);
            }
        }
        catch (Exception e) 
        {
            System.err.format("Error reading %s : %s\n", filePath, e.getMessage());
        }        
        return res;
    }

    //filePath is database/PSSM/HLA-xxxxx_N.pssm
    private static WeightMatrix ParseOnePSSM(String filePath)
    {
        WeightMatrix res = new WeightMatrix();
        try
        {
            BufferedReader reader = new BufferedReader(new FileReader(filePath));
            
            String headerStr = reader.readLine();
            String[] header = headerStr.split(" ");
            if (header.length != 2)
            {
                System.err.format("Invalid header '%s' of %s\n", headerStr, filePath);
            }
            
            int colsQnty = Integer.parseInt(header[1]);    //matrix columns count (peptide length)
            for (int iCol=0; iCol<colsQnty; iCol++)
            {
                res.add(new WeightMatrixColumn());
            }

            int iRow = 0;
            while (reader.ready() && (iRow < alphabet.length()))
            {
                String[] row = reader.readLine().split("\t");
                if (row.length > colsQnty)
                {
                    System.err.format("Invalid PSSM size : [%d by %d]\n", iRow, row.length);
                    return res;
                }

                for (int iCol=0; iCol<colsQnty; iCol++)
                {
                    res.get(iCol).put(new Character(alphabet.charAt(iRow)), Double.parseDouble(row[iCol]));
                }

                iRow +=1;
            }

            if (iRow != alphabet.length())
            {
                System.err.format("Invalid PSSM of %d rows\n", iRow);
            }
        } 
        catch (Exception e) 
        {
            System.err.format("Error reading %s : %s\n", filePath, e.getMessage());
        }        
        return res;
    }
}

class PSSMHCpan
{
    static String peptide = "CLMPCGRRQ";

    static double score_max = 0.8;
    static double score_min = 0.8 * (1 - Math.log(50000) / Math.log(500));
    static double score_range = score_max - score_min;

    
    private String peptidesFilename, alleleName, PSSMlistFilename;
    private int peptideLength;
    private AllelePair al = new AllelePair();


    WeightMatrices pssms = null;

    PSSMHCpan(String[] args)
    {
        if (args.length != 4)
        {
            System.out.format("Usage : java PSSMHCpan peptides_list.fa <peptide_length> <allele name> database/PSSM/pssm_file.list\n");
            System.exit(1);
        }
        
        peptidesFilename = args[0];
        al.pepLength = Integer.parseInt(args[1]);
        al.alName = args[2];
        PSSMlistFilename = args[3];

        pssms = PSSMParser.ParsePSSMfileList(PSSMlistFilename);
    }

    double score_one_peptide()
    {
        WeightMatrix pssm = pssms.get(al);

        double score = 0;
        for (int ch_pos=0; ch_pos<peptide.length(); ch_pos++)
        {
            Map<Character, Double> weightRow = pssm.get(ch_pos);
            Double weight = weightRow.get(peptide.charAt(ch_pos));
            score += weight;
        }

        score = score / peptide.length();
        score = Math.max(Math.min(score, score_max), score_min);

        return Math.pow(50000, (score_max - score)/score_range);
    }

    public static void main(String[] args) 
    {
        PSSMHCpan app = new PSSMHCpan(args);

        try
        {
            PrintWriter writer = new PrintWriter("res.txt", "UTF-8");
            double ic50 = app.score_one_peptide();
            writer.format("%s %f\n", peptide, ic50);
            writer.close();
        } 
        catch (IOException e) {
           // do something
        }
    }
}