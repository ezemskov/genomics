package org.PSSMHC;

import org.junit.Test;
import static org.junit.Assert.*;

public class PeptideGenTest
{
    PeptideGen _gen = null;

    void setUp(int pepLen)    //NB called manually
    {
        _gen = new PeptideGen(pepLen);
    }
        
    void testGenerateRange(int start, int qnty, String expected)
    {
        String[] expectedArray = expected.split(",");
        for (int i=start; i<start+qnty; ++i)
        {
            assertEquals(expectedArray[i-start], _gen.Generate(i));
        }
    }
    
    @Test
    public void test8meers() 
    {
        setUp(8);
        testGenerateRange(19, 3, 
            "AAAAAAAY," + 
            "AAAAAACA," + 
            "AAAAAACC,");
        
        testGenerateRange(7999, 2, 
            "AAAAAYYY," + 
            "AAAACAAA,");
    }
    
    @Test
    public void testRanges1() 
    {
        setUp(9);
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
        setUp(9);
        testGenerateRange(7998, 4, 
            "AAAAAAYYW," + 
            "AAAAAAYYY," + 
            "AAAAACAAA," + 
            "AAAAACAAC,");
    }

    @Test
    public void testRanges3() 
    {
        setUp(9);
        testGenerateRange(399, 4, 
            "AAAAAAAYY," + 
            "AAAAAACAA," + 
            "AAAAAACAC," + 
            "AAAAAACAD,");
    }

    @Test
    public void testInterruptedRanges() 
    {
        setUp(9);
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
        new PeptideGen().Generate(-10);
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testIdxOutOfRange2() 
    {
        new PeptideGen(9).Generate((long)512E9 + 1);
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testIdxOutOfRange3() 
    {
        new PeptideGen(8).Generate((long)256E8 + 1);
    }
    
    @Test(expected = IndexOutOfBoundsException.class)
    public void testIdxOutOfRange4() 
    {
        new PeptideGen(10).Generate((long)2048E11 + 1);
    }
}
