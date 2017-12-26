package org.PSSMHC;

import java.io.Serializable;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.sql.Row;

public class Impl
{
    public static class Consts
    {
        public static final String alphabet = "ACDEFGHIKLMNPQRSTVWY";
        public static final int aLen = alphabet.length();
        public static final String alphabetRegex = "[" + Consts.alphabet + "]+";
    }

    public static class PSSMHCpanSparkFunc 
                            extends PSSMHCpan
                            implements Serializable,
                            Function<String,ScoredPeptide>
    {
        public PSSMHCpanSparkFunc(Xml.PSSMCfg cfg) throws Exception
        {
            super(cfg);
        }
            
        @Override
        public ScoredPeptide call(String peptide)
        {
            return new ScoredPeptide(peptide, ScoreOnePeptide(peptide));
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
            return (scoreThreshold > 0)  ? (scPep.score < scoreThreshold) : true;
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
