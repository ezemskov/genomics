package org.PSSMHC;

import java.lang.Math;
import java.util.*;
import java.lang.System;
import java.lang.String;
import java.lang.RuntimeException;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.Serializable;

class WeightMatrixColumn extends HashMap<Character, Double> implements Serializable {}
class WeightMatrix extends ArrayList<WeightMatrixColumn> implements Serializable {}

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

class ScoredPeptide implements Serializable
{
    ScoredPeptide(String peptide, double ic50) 
    { 
        this.peptide = peptide; 
        this.ic50 = ic50;
    }
    
    public String toString()
    {
        return new String().format("%s,%.0f", peptide, ic50);
    }
        
    public String peptide; 
    public double ic50; 
}

class Consts
{
    public static String alphabet = "ACDEFGHIKLMNPQRSTVWY";
    public static int aLen = alphabet.length();
}

class PSSMParser
{
    private static String alphabetRegex = "[" + Consts.alphabet + "]+";

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
            System.err.format("Error reading %s : %s\n", pssmListPath, e.getMessage());
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
                while (reader.ready() && (iRow < Consts.aLen))
            {
                String[] row = reader.readLine().split("\t");
                if (row.length > colsQnty)
                {
                    System.err.format("Invalid PSSM size : [%d by %d]\n", iRow, row.length);
                    return res;
                }

                for (int iCol=0; iCol<colsQnty; iCol++)
                {
                    res.get(iCol).put(new Character(Consts.alphabet.charAt(iRow)), Double.parseDouble(row[iCol]));
                }

                iRow +=1;
            }

            if (iRow != Consts.aLen)
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

class PSSMHCpan implements Serializable
{
    public static String CmdlineHelpStr = "Usage : "
        + "spark-submit --class org.PSSMHC.PSSMHCSpark --master spark://<IP>:7077 /path/to/PSSMHCpan-1.0.jar <peptide_length> <allele_name> path/to/pssm_file.list <peptide start idx> <peptide qnty> <partitions> doSrcPersistCount doBinderPersist doBinderStore doBinderCount\n"
        + "e.g. spark-submit --class org.PSSMHC.PSSMHCSpark --master spark://192.168.56.1:7077 ./PSSMHCpan-1.0.jar 9 HLA-A0201 database/PSSM/pssm_file.list 1 1000 4 0 1 1 1";

    
    private static double score_max = 0.8;
    private static double score_min = 0.8 * (1 - Math.log(50000) / Math.log(500));
    private static double score_range = score_max - score_min;
    
    protected WeightMatrix pssm = null;
    
    public int InitFromCmdline(String[] args)
    {
        if (args.length < 3)
        {
            throw new RuntimeException(CmdlineHelpStr);
        }
        
        AllelePair ap = new AllelePair();
        ap.pepLength = Integer.parseInt(args[0]);
        ap.alName = args[1];
        String PSSMlistFilename = args[2];

        pssm = PSSMParser.FindAndParsePSSM(PSSMlistFilename, ap);
        
        return 3;
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
}

class PSSMHCpanFasta extends PSSMHCpan
{
    public static String CmdlineHelpStr = "Usage : java org.PSSMHC.PSSMHCpanJava peptides_list.fa <peptide_length> <allele name> database/PSSM/pssm_file.list\n";

    public transient ArrayList<String> peptides = new ArrayList<String>();

    public int InitFromCmdline(String[] args)
    {
        if (args.length < 4)
        {
            throw new RuntimeException(CmdlineHelpStr);
        }
        
        String peptidesFilename = args[0];
        AllelePair ap = new AllelePair();
        ap.pepLength = Integer.parseInt(args[1]);
        ap.alName = args[2];
        String PSSMlistFilename = args[3];

        pssm = PSSMParser.FindAndParsePSSM(PSSMlistFilename, ap);

        ArrayList<ScoredPeptide> scPeptides = PSSMParser.ParseFasta(peptidesFilename);
        for (ScoredPeptide scp : scPeptides)
        {
            peptides.add(scp.peptide);
        }
        
        return 4;
    }

    public void ScoreAllPeptides()
    {
        for (String pep : peptides)
        {
            double ic50 = ScoreOnePeptide(pep);
            System.out.println(new ScoredPeptide(pep, ic50));
        }
    }
}

final class PSSMHCpanJava
{    
    public static void main(String[] args) 
    {
        try
        {
            PSSMHCpanFasta app = new PSSMHCpanFasta();
            app.InitFromCmdline(args);
            app.ScoreAllPeptides();
        }
        catch (Exception ex)
        {
            System.err.print(ex.getMessage() + "\n");
        }
    }
}