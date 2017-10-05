package org.PSSMHC;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class PeptideGenTest
{
    PeptideGen _gen;
    
    @Before
    public void setUp()
    {
        _gen = new PeptideGen();
   }

    public void testGenerateRange(int start, int qnty, String expected)
    {
        String[] expectedArray = expected.split(",");
        for (int i=start; i<start+qnty; ++i)
        {
            assertEquals(expectedArray[i-start], _gen.Generate(i));
        }
    }

    @Test
    public void testRanges1() 
    {
        testGenerateRange(18, 5, 
            "AAAAAAAAW," + 
            "AAAAAAAAY," + 
            "AAAAAAACA," + 
            "AAAAAAACC," + 
            "AAAAAAACD,");
    }

    @Test
    public void testRanges2() 
    {
        testGenerateRange(7998, 4, 
            "AAAAAAYYW," + 
            "AAAAAAYYY," + 
            "AAAAACAAA," + 
            "AAAAACAAC,");
    }

    @Test
    public void testRanges3() 
    {
        testGenerateRange(399, 4, 
            "AAAAAAAYY," + 
            "AAAAAACAA," + 
            "AAAAAACAC," + 
            "AAAAAACAD,");
    }

    @Test
    public void testInterruptedRanges() 
    {
        testGenerateRange(1199, 2, 
            "AAAAAADYY," + 
            "AAAAAAEAA,");

        testGenerateRange(2000, 2, 
            "AAAAAAGAA," + 
            "AAAAAAGAC,");

        testGenerateRange(19*(160000 + 20 + 1), 3, 
            "AAAAYAAYY," + 
            "AAAAYACAA," + 
            "AAAAYACAC,");

    }
    
    @Test(expected = IndexOutOfBoundsException.class)
    public void testIdxOutOfRange1() 
    {
        _gen.Generate(-10);
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testIdxOutOfRange2() 
    {
        _gen.Generate((long)512E9 + 1);
    }
    
}
