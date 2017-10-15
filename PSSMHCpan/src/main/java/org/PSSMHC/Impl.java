package org.PSSMHC;

import java.io.Serializable;
import java.util.HashMap;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.sql.Row;

public class Impl
{
    public static class ScoreFunc extends PSSMHCpan
                    implements Serializable,
                    Function<String,ScoredPeptide>
    {
        @Override
        public ScoredPeptide call(String peptide)
        {
            return new ScoredPeptide(peptide, ScoreOnePeptide(peptide));
        }
    }

    public static class BloomFunc extends PeptideBloomFilter
                    implements Serializable,
                    VoidFunction<ScoredPeptide>
    {
        @Override
        public void call(ScoredPeptide scPep)
        {
            Add(scPep.peptide);
        }
    }

    public static class GenFunc extends PeptideGen
                    implements Serializable, 
                    MapFunction<Row, String>
    {
        public String call(Row idx)
        {
            return Generate(idx.getLong(0));
        }
    }

    public static class Ic50FilterFunc implements Serializable,
                         Function<ScoredPeptide,Boolean>
    {
        public Ic50FilterFunc(double ic50Threshold_)
        {
            ic50Threshold = ic50Threshold_;
        }

        public Boolean call(ScoredPeptide scPep)
        {
            return (scPep.ic50 < ic50Threshold);
        }

        double ic50Threshold;
    }

    public static class CmdlineCfg
    {
        public long start;
        public long end;
        public int partitions;
        public boolean doScore, doBinderPersist, doBinderStore, doBinderCount;
        public int ic50Threshold;
        private int firstArgIdx;
        
        public CmdlineCfg(String[] args, int firstArgIdx_)
        {
            firstArgIdx = firstArgIdx_;
            if (args.length < firstArgIdx+8)
            {
                throw new RuntimeException(PSSMHCpan.CmdlineHelpStr);
            }

            start =       ParseLongWithSuffix(args[firstArgIdx]);
            end = start + ParseLongWithSuffix(args[firstArgIdx+1]);    
            partitions =  Integer.parseUnsignedInt(args[firstArgIdx+2]);

            doScore     = (Integer.parseInt(args[firstArgIdx+3]) == 1);
            doBinderPersist   = (Integer.parseInt(args[firstArgIdx+4]) == 1);
            doBinderStore     = (Integer.parseInt(args[firstArgIdx+5]) == 1);
            doBinderCount     = (Integer.parseInt(args[firstArgIdx+6]) == 1);
            ic50Threshold = Integer.parseInt(args[firstArgIdx+7]);
        }
        
        public int NextArgIdx() { return firstArgIdx + 8; }
            
        private static long ParseLongWithSuffix(String val)
        {
            int suffixPos = val.length()-1;
            Character suffix = val.toUpperCase().charAt(suffixPos);
            Long multiplier = suffixes.get(suffix);
            if ((val.length() < 2) || (multiplier == null))
            {
                return Long.parseUnsignedLong(val);
            }

            return multiplier * Long.parseUnsignedLong(val.substring(0, suffixPos));
        }

        private static final HashMap<Character, Long> suffixes;
        static
        {
            suffixes = new HashMap<>();
            suffixes.put('K', (long)1E3);
            suffixes.put('M', (long)1E6);
            suffixes.put('G', (long)1E9);
            suffixes.put('T', (long)1E12);
        }
    }    
}
