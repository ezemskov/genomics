package org.PSSMHC;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.sql.Row;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;
import org.w3c.dom.Document; 
import org.w3c.dom.Element; 
import org.w3c.dom.NodeList; 
import org.w3c.dom.Node; 

public class Impl
{
    public static class PSSMHCpanSparkFunc 
                            extends PSSMHCpan
                            implements Serializable,
                            Function<String,ScoredPeptide>
    {
        public PSSMHCpanSparkFunc(String xmlFilename) throws Exception
        {
            super(xmlFilename);
        }
            
        @Override
        public ScoredPeptide call(String peptide)
        {
            return new ScoredPeptide(peptide, ScoreOnePeptide(peptide));
        }
    }

    public static class PeptideBloomSparkFunc 
                            extends PeptideBloomFilter
                            implements Serializable,
                            VoidFunction<ScoredPeptide>
    {
        @Override
        public void call(ScoredPeptide scPep)
        {
            Add(scPep.peptide);
        }
    }

    public static class PeptideGenSparkFunc 
                            extends PeptideGen
                            implements Serializable, 
                            MapFunction<Row, String>
    {
        @Override
        public String call(Row idx)
        {
            return Generate(idx.getLong(0));
        };
    }

    public static class ScoreFilterSparkFunc 
                            implements Serializable,
                            Function<ScoredPeptide,Boolean>
    {
        public ScoreFilterSparkFunc(double scoreThreshold_)
        {
            scoreThreshold = scoreThreshold_;
        }

        @Override
        public Boolean call(ScoredPeptide scPep)
        {
            return (scPep.score < scoreThreshold);
        }

        double scoreThreshold;
    }

    public static class ScoredPeptide implements Serializable
    {
        public ScoredPeptide(String peptide, double score) 
        { 
            this.peptide = peptide; 
            this.score = score;
        }

        @Override
        public String toString()
        {
            return String.format("%s,%.2f", peptide, score);
        }

        public String peptide; 
        public double score; 
    }
}
