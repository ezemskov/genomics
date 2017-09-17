package org.PSSMHC;

import java.lang.Math;
import java.util.*;
import java.lang.System;
import java.lang.String;
import java.lang.RuntimeException;
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

class WeightMatrixPathMap extends HashMap<AllelePair, String> {}
class WeightMatrices extends HashMap<AllelePair, WeightMatrix> {}

class ScoredPeptide
{
    ScoredPeptide(String peptide, double ic50) 
    { 
        this.peptide = peptide; 
        this.ic50 = ic50;
    }
    
    public String peptide; 
    public double ic50; 
}

class PSSMParser
{
    public static String alphabet = "ACDEFGHIKLMNPQRSTVWY";
    private static String alphabetRegex = "[" + alphabet + "]+";

    public static ArrayList<ScoredPeptide> ParseFasta(String filePath)
    {
        ArrayList<ScoredPeptide> res = new ArrayList<>();
        try
        {
            BufferedReader reader = new BufferedReader(new FileReader(filePath));
            while (reader.ready())
            {
                String recHeader = reader.readLine();
                if (recHeader.charAt(0) != '>')
                {
                    System.err.format("Wrong Fasta record header '%s' in %s, ignore rest of the file\n", recHeader, filePath);
                    return res;
                }
                
                String peptide = reader.readLine();
                if (peptide.matches(alphabetRegex))
                {
                    res.add(new ScoredPeptide(peptide, -1.0));
                }
                else
                {
                    System.err.format("Skip wrong Fasta record '%s' in %s\n", peptide, filePath);
                }
            }
        }
        catch (Exception e) 
        {
            System.err.format("Error reading %s : %s\n", filePath, e.getMessage());
        }        
        return res;
    }
        
    //filePath is pssm_file.list
    //ap is a 'allele name/length' pair
    public static WeightMatrix FindAndParsePSSM(String pssmListPath, AllelePair ap)
    {
        WeightMatrix res = null;
        try
        {
            BufferedReader reader = new BufferedReader(new FileReader(pssmListPath));
            while (reader.ready())
            {
                String rowStr = reader.readLine();
                String[] row = rowStr.split("\t");
                if (row.length != 3)
                {
                    System.err.format("Skip invalid line PSSM list line '%s'\n", rowStr);
                    continue;
                }
                
                if (ap.alName.equals(row[1]) && 
                    (ap.pepLength == Integer.parseInt(row[2])))
                {
                    String pssmFilePath = row[0];
                    res = ParsePSSM(pssmFilePath);
                }
            }
        }
        catch (Exception e) 
        {
            System.err.format("Error reading %s : %s\n", filePath, e.getMessage());
        }        
        return res;
    }

    //filePath is database/PSSM/HLA-xxxxx_N.pssm
    private static WeightMatrix ParsePSSM(String filePath)
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

public class PSSMHCpan
{
    private static double score_max = 0.8;
    private static double score_min = 0.8 * (1 - Math.log(50000) / Math.log(500));
    private static double score_range = score_max - score_min;
    
    private WeightMatrix pssm = null;
    public ArrayList<ScoredPeptide> peptides = new ArrayList<ScoredPeptide>();
    public ArrayList<String> peptides2 = new ArrayList<String>();
    
    public void InitFromCmdline(String[] args)
    {
        if (args.length != 4)
        {
            throw new RuntimeException("Usage : java PSSMHCpan peptides_list.fa <peptide_length> <allele name> database/PSSM/pssm_file.list\n");
        }
        
        String peptidesFilename = args[0];
        AllelePair ap = new AllelePair();
        ap.pepLength = Integer.parseInt(args[1]);
        ap.alName = args[2];
        String PSSMlistFilename = args[3];

        pssm = PSSMParser.FindAndParsePSSM(PSSMlistFilename, ap);

        peptides = PSSMParser.ParseFasta(peptidesFilename);
        for (ScoredPeptide scp : peptides)
        {
            peptides2.add(scp.peptide);
        }
    }
    
    public double ScoreOnePeptide(String peptide)
    {
        if (pssm == null)
        {
            throw new RuntimeException("PSSM is not initialized");
        }
        double score = 0;
        for (int ch_pos=0; ch_pos<peptide.length(); ch_pos++)
        {
            Map<Character, Double> weightRow = pssm.get(ch_pos);
            if(weightRow == null) {
                throw new RuntimeException("Invalid PSSM width");
            }
            Double weight = weightRow.get(peptide.charAt(ch_pos));
            assert(weight != null); //assumes PSSM is parsed and peptide is validated
            score += weight;
        }

        score = score / peptide.length();
        score = Math.max(Math.min(score, score_max), score_min);

        return Math.pow(50000, (score_max - score)/score_range);
    }

    public void ScoreAllPeptides()
    {
        for (ScoredPeptide scPep : peptides)
        {
            scPep.ic50 = ScoreOnePeptide(scPep.peptide);
            System.out.format("%s %f\n", scPep.peptide, scPep.ic50);
        }
    }
}

class MainWrapper
{    
    public static void main(String[] args) 
    {
        try
        {
            PSSMHCpan app = new PSSMHCpan();
            app.InitFromCmdline(args);
            app.ScoreAllPeptides();
        }
        catch (Exception ex)
        {
            if (!ex.getMessage().isEmpty())
            {
                System.err.print(ex.getMessage());
            }
        }
    }
}