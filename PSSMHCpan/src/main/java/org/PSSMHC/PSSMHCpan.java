package org.PSSMHC;

import java.util.*;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.Serializable;
import org.w3c.dom.Element;

class WeightMatrixColumn extends HashMap<Character, Double> implements Serializable {}
class WeightMatrix extends ArrayList<WeightMatrixColumn> implements Serializable {}

class AllelePair
{
    @Override
    public int hashCode()
    {
        return String.format("%s_%d", alName, pepLength).hashCode();
    }

    @Override
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

class Consts
{
    public static final String alphabet = "ACDEFGHIKLMNPQRSTVWY";
    public static int aLen = alphabet.length();
}

class PSSMParser
{
    private static final String alphabetRegex = "[" + Consts.alphabet + "]+";
        
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
    private static final double ScoreMax = 0.8;
    private static final double ScoreMin = 0.8 * (1 - Math.log(50000) / Math.log(500));
    private static final double ScoreRange = ScoreMax - ScoreMin;
    
    protected WeightMatrix pssm = null;
    
    PSSMHCpan(String xmlFilename) throws Exception
    {
        Element root = Xml.Utils.parseXml(xmlFilename);
        Element elem = Xml.Utils.getChildElem(root, "PSSMHC");
        if (elem == null) { return; }
        
        AllelePair ap = new AllelePair();
        ap.pepLength = Integer.parseInt(elem.getAttribute("peptideLength"));
        ap.alName = elem.getAttribute("alleleName");

        String PSSMlistFilename = elem.getAttribute("fileList");
        pssm = PSSMParser.FindAndParsePSSM(PSSMlistFilename, ap);
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
        score = Math.max(Math.min(score, ScoreMax), ScoreMin);

        return Math.pow(50000, (ScoreMax - score)/ScoreRange);
    }
}
