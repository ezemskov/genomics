package org.PeptideClustering;

import org.junit.Before;
import org.junit.Test;
import java.util.Random;
import org.PSSMHC.PeptideGen;
import java.lang.System;

public class PeptideSimilarityBench
{
    PeptideSimilarity simCalc = null;
    
    @Before
    public void setUp() throws Exception
    {
        simCalc = new PeptideSimilarity().SetMatrix(new SubstMatrices.Blosum62());
    }

    @Test
    public void benchmarkEmptyCall()
    {
        System.out.println("empty");
        benchmark((p1, p2) -> { return 0.0; }, false);
    }

    @Test
    public void benchmarkGen()
    {
        System.out.println("gen only");
        benchmark((p1, p2) -> { return 0.0; }, true);
    }

    @Test
    public void benchmarkPosDiff()
    {
        System.out.println("gen and posDiff");
        benchmark((p1, p2) -> { return simCalc.posDiff(p1, p2); }, true);
    }
    
    @Test
    public void benchmarkPosDiffOnly()
    {
        System.out.println("posDiff only");
        benchmark((p1, p2) -> { return simCalc.posDiff(p1, p2); }, false);
    }

    @Test
    public void benchmarkSimilarity()
    {
        System.out.println("similarity only");
        benchmark((p1, p2) -> { return simCalc.similarity(p1, p2); }, false);
    }

    @Test
    public void benchmarkNW()
    {
        System.out.println("NeedlemanWunsch");
        benchmark((p1, p2) -> { return NeedlemanWunsch.calcNWscore(p1, p2).score; }, false);
    }

    
    interface IPeptideDistance {
       public double call(String p1, String p2);
    }
    
    private void benchmark(IPeptideDistance func, boolean doRegenerate)
    {
        Random randGen = new Random();
        PeptideGen pepGen = PeptideGen.CreateDefLen();
        
        final int IterationsCount = 1000000;
        final long startTime = System.currentTimeMillis();
        int interruptedCalcCount = 0;
        
        String p1 = pepGen.Generate(Math.abs(randGen.nextInt()));
        String p2 = pepGen.Generate(Math.abs(randGen.nextInt()));
        for (int i=0; i<IterationsCount; ++i)
        {
            if (doRegenerate)
            {
                p1 = pepGen.Generate(Math.abs(randGen.nextInt()));
                p2 = pepGen.Generate(Math.abs(randGen.nextInt()));
            }
            double res = func.call(p1, p2);
            if (doRegenerate && Double.isNaN(res)) {
                ++interruptedCalcCount;
            }
        }

        final long durationMsec = System.currentTimeMillis() - startTime;
        final double durationAvgUsec = (double)(1E3 * durationMsec/IterationsCount);
        final double interruptedCalcFrac = interruptedCalcCount/(double)IterationsCount;
        System.out.format("Duration total %d msec avg %.3f usec interrupted %.6f\n", 
            durationMsec, durationAvgUsec, interruptedCalcFrac);
    }
}
