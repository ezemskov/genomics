package org.PeptideClustering;

import org.junit.Before;
import org.junit.Test;
import java.util.Random;
import java.lang.AssertionError;
import org.PSSMHC.PeptideGen;
import static org.junit.Assert.*;

public class PeptideSimilarityTest
{
    PeptideSimilarity simCalc = null;
    //wide threshold accounts for difference between integer and real subst matrices
    static double Threshold = 0.02;
    static String oxytocin = "CYIQNCPLG";
    
    @Before
    public void setUp() throws Exception
    {
        //alternatively : 
        //simCalc = new PeptideSimilarity().SetMatrix(SubstMatrices.get("blosum62"));
        simCalc = new PeptideSimilarity().SetMatrix(new SubstMatrices.Blosum62());
    }

    @Test(expected = AssertionError.class)
    public void testDefaultInit()
    {
        assertEquals(1.0, new PeptideSimilarity().similarity(oxytocin, oxytocin), Threshold);
    }

    @Test
    public void testSimilarityToSelf()
    {
        assertEquals(1.0, simCalc.similarity(oxytocin, oxytocin), Threshold);
    }

    @Test
    public void testSimilarityToOther()
    {
        String other = "SQPLSALWG";
        assertEquals(-0.114, simCalc.similarity(oxytocin, other), 0.02);
        assertEquals(0.504, simCalc.similarity("RPVDQNTQS", "KAFDRNTES"), Threshold);
    }

    @Test
    public void testAvgRandSimilarity()
    {
        for (int i=0; i<100; ++i)
        {
            testAvgRandSimilarityOnce();
        }
    }

    @Test
    public void testSubstMatrix() throws Exception
    {
        SubstMatrix matrix = SubstMatrices.get("blosum62");
        
        assertEquals(matrix.get('C').get('M'), -1, 0.1);
        assertEquals(matrix.get('K').get('R'),  2, 0.1);
    }
    
    private void testAvgRandSimilarityOnce()
    {
        Random randGen = new Random();
        PeptideGen pepGen = PeptideGen.CreateDefLen();
        int qnty = 1000;
        double totalSim = 0.0;
        double randomThreshold = 0.05;
        
        for (int i=0; i<qnty; ++i)
        {
            String p1 = pepGen.Generate(Math.abs(randGen.nextInt()));
            String p2 = pepGen.Generate(Math.abs(randGen.nextInt()));
            totalSim += simCalc.similarity(p1, p2);
        }
        
        double avgSim = totalSim/qnty;
        assertEquals(0.1, avgSim, randomThreshold);
    }

    @Test
    public void testMapPosDiffBlosum62()
    {
        MaxPosDiff m = new MaxPosDiff(new SubstMatrices.Blosum62());

        assertEquals(PeptideGen.PepLenDefault, m.lowerBound(-1.0));
        assertEquals(PeptideGen.PepLenDefault, m.lowerBound(-0.15));
        assertEquals(PeptideGen.PepLenDefault, m.lowerBound(-0.1));
        assertEquals(PeptideGen.PepLenDefault, m.lowerBound(0.0));

        assertEquals(8, m.lowerBound(0.1));
        assertEquals(8, m.lowerBound(0.16));
        assertEquals(8, m.lowerBound(0.2));
        assertEquals(8, m.lowerBound(0.3));
        assertEquals(8, m.lowerBound(0.35));
        assertEquals(7, m.lowerBound(0.4));
        assertEquals(7, m.lowerBound(0.5));
        assertEquals(7, m.lowerBound(0.6));
        assertEquals(7, m.lowerBound(0.65));
        assertEquals(6, m.lowerBound(0.7));
        assertEquals(6, m.lowerBound(0.8));
        assertEquals(6, m.lowerBound(0.9));
        assertEquals(5, m.lowerBound(1.0));
        assertEquals(5, m.lowerBound(1.05));
        assertEquals(5, m.lowerBound(2.0));
    }
}
