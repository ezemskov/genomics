package org.PSSMHC;

import java.util.*;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.Serializable;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import org.w3c.dom.Element;

class WeightMatrixCol extends HashMap<Character, Double> implements Serializable
{
    WeightMatrixCol() { super(Impl.Consts.aLen, 1.0f); }
}
class WeightMatrix extends ArrayList<WeightMatrixCol> implements Serializable {}

//{allele name, peptide size} -> path to weight matrix file
class WeightMatrixPathMap extends HashMap<Xml.StringIntPair, String> {}
//{allele name, peptide size} -> weight matrix
class WeightMatrices extends HashMap<Xml.StringIntPair, WeightMatrix> {}

class PSSMParser
{
    public static WeightMatrix FindAndParsePSSM(Xml.PSSMCfg cfg) throws Exception
    {
        final Path pssmListPath = FileSystems.getDefault().getPath(cfg.pathPrefix, cfg.pssmListRelPath);
        BufferedReader reader = new BufferedReader(new FileReader(pssmListPath.toFile()));
        while (reader.ready())
        {
            String rowStr = reader.readLine();
            String[] row = rowStr.split("\t");
            if (row.length != 3)
            {
                System.err.format("Skip invalid line PSSM list line '%s'\n", rowStr);
                continue;
            }

            if (cfg.allele.equals(row[1]) && 
               (cfg.peptideLength == Integer.parseInt(row[2])))
            {
                final String pssmFilePathStr = row[0];
                final Path pssmPath = FileSystems.getDefault().getPath(cfg.pathPrefix, pssmFilePathStr);
                return ParsePSSM(pssmPath);
            }
        }
        
        throw new Exception(String.format("PSSM for allele %s length %d was not found in %s", 
            cfg.allele, cfg.peptideLength, pssmListPath.toString()));
    }

    //filePath is database/PSSM/HLA-xxxxx_N.pssm
    private static WeightMatrix ParsePSSM(Path filePath) throws Exception
    {
        WeightMatrix res = new WeightMatrix();
        BufferedReader reader = new BufferedReader(new FileReader(filePath.toFile()));

        String headerStr = reader.readLine();
        String[] header = headerStr.split(" ");
        if (header.length != 2)
        {
            throw new Exception(String.format(
                "Invalid header '%s' of %s", headerStr, filePath.toString()));
        }

        int colsQnty = Integer.parseInt(header[1]);    //matrix columns count (peptide length)
        for (int iCol=0; iCol<colsQnty; iCol++)
        {
            res.add(new WeightMatrixCol());
        }

        int iRow = 0;
        while (reader.ready() && (iRow < Impl.Consts.aLen))
        {
            String[] row = reader.readLine().split("\t");
            if (row.length != colsQnty)
            {
                throw new Exception(String.format(
                    "Invalid PSSM of %d columns in %s", row.length, filePath.toString()));
            }

            for (int iCol=0; iCol<colsQnty; iCol++)
            {
                res.get(iCol).put(Impl.Consts.alphabet.charAt(iRow), Double.parseDouble(row[iCol]));
            }

            iRow +=1;
        }

        if (iRow != Impl.Consts.aLen)
        {
            throw new Exception(String.format(
                "Invalid PSSM of %d rows in %s", iRow, filePath.toString()));
        }

        return res;
    }
}

class PSSMHCpan implements Serializable
{
    
    private static final double ScoreMax = 0.8;
    private static final double ScoreMin = 0.8 * (1 - Math.log(50000) / Math.log(500));
    private static final double ScoreRange = ScoreMax - ScoreMin;
    
    protected WeightMatrix _pssm = null;
    
    PSSMHCpan(Xml.PSSMCfg cfg) throws Exception
    {
        _pssm = PSSMParser.FindAndParsePSSM(cfg);
        assert(_pssm != null);  //assured in PSSMParser
    }
        
    public double ScoreOnePeptide(String peptide)
    {
        if (_pssm.size() != peptide.length())
        {
            return Double.NaN;
        }
        
        double score = 0;
        for (int ch_pos=0; ch_pos<peptide.length(); ch_pos++)
        {
            Map<Character, Double> weightCol = _pssm.get(ch_pos);
            Double weight = weightCol.get(peptide.charAt(ch_pos));
            if (weight == null)
            {
                return Double.NaN;
            }
            score += weight;
        }

        score = score / peptide.length();
        score = Math.max(Math.min(score, ScoreMax), ScoreMin);

        return Math.pow(50000, (ScoreMax - score)/ScoreRange);
    }
}
