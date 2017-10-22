package org.PeptideClustering;

import org.junit.Before;
import org.junit.Test;
import java.util.Random;
import java.lang.AssertionError;
import org.PSSMHC.PeptideGen;
import static org.junit.Assert.*;
import java.util.stream.LongStream;
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
    public void benchmarkGenOnly()
    {
        System.out.println("generate only");
        benchmark((p1, p2) -> { return Double.NaN; } );
    }

    @Test
    public void benchmarkSimilarity()
    {
        System.out.println("similarity");
        benchmark((p1, p2) -> { return simCalc.similarity(p1, p2); } );
    }

    @Test
    public void benchmarkDissimilarity1()
    {
        System.out.println("dissimilarity 0.1");
        benchmark((p1, p2) -> { return simCalc.dissimilarity(p1, p2, 1/0.1); } );
    }

    @Test
    public void benchmarkDissimilarity2()
    {
        System.out.println("dissimilarity 0.2");
        benchmark((p1, p2) -> { return simCalc.dissimilarity(p1, p2, 1/0.2); } );
    }

    @Test
    public void benchmarkDissimilarity5()
    {
        System.out.println("dissimilarity 0.5");
        benchmark((p1, p2) -> { return simCalc.dissimilarity(p1, p2, 1.0/0.5); } );
    }

    @Test
    public void benchmarkDissimilarity7()
    {
        System.out.println("dissimilarity 0.7");
        benchmark((p1, p2) -> { return simCalc.dissimilarity(p1, p2, 1.0/0.7); } );
    }
    
    @Test
    public void benchmarkDissimilarity8()
    {
        System.out.println("dissimilarity 0.8");
        benchmark((p1, p2) -> { return simCalc.dissimilarity(p1, p2, 1.0/0.8); } );
    }

    @Test
    public void benchmarkDissimilarity9()
    {
        System.out.println("dissimilarity 0.9");
        benchmark((p1, p2) -> { return simCalc.dissimilarity(p1, p2, 1.0/0.9); } );
    }
    
    interface IPeptideDistance {
       public double call(String p1, String p2);
    }
    
    private void benchmark(IPeptideDistance func)
    {
        Random randGen = new Random();
        PeptideGen pepGen = new PeptideGen();
        
        final int IterationsCount = 1000000;
        final long startTime = System.currentTimeMillis();
        int interruptedCalcCount = 0;
        
        for (int i=0; i<IterationsCount; ++i)
        {
            String p1 = pepGen.Generate(Math.abs(randGen.nextInt()));
            String p2 = pepGen.Generate(Math.abs(randGen.nextInt()));
            double res = func.call(p1, p2);
            if (Double.isNaN(res)) {
                ++interruptedCalcCount;
            }
        }

        final long durationMsec = System.currentTimeMillis() - startTime;
        final double durationAvgUsec = (double)(1E3 * durationMsec/IterationsCount);
        final double interruptedCalcFrac = interruptedCalcCount/(double)IterationsCount;
        System.out.format("Duration total %d msec avg %.2f usec interrupted %.6f\n", 
            durationMsec, durationAvgUsec, interruptedCalcFrac);
    }
}
