package org.PSSMHC;

import java.lang.StringBuilder;

public class PeptideGen
{
    public static int pepLen = 9;
    private static long maxIndex = (long)Math.pow(Impl.Consts.aLen, pepLen) - 1; //20^9 - 1
    
    private StringBuilder lastPepStr = null;
    private int[]  lastPepCharPos = new int[pepLen];
    private long   lastIdxPep = -2;
    
    /**
     * @param index number of peptide in lexicographical order, 
     * i.e. AAAAAAAAA is #0, AAAAAAAAC is #1, and YYYYYYYYY is #20^9-1
     * @return peptide 
     */
    public String Generate(long idxPep)
    {
        if ((idxPep < 0) || (idxPep > maxIndex))
        {
            throw new IndexOutOfBoundsException(String.format("Peptide index %d is out of bounds", idxPep));
        }

        return (idxPep == lastIdxPep+1) ? GenerateNext() : GenerateFirst(idxPep);
    }
    
    private String GenerateNext()
    {
        lastIdxPep += 1;
        
        int idxChar = pepLen-1;
        //find left bound of '....YYY' chars sequence, starting from right, 
        //and set it to '....AAAA' (i.e. perform rollover)
        for (; (idxChar >= 0) && (lastPepCharPos[idxChar] == Impl.Consts.aLen-1); --idxChar)
        {
            lastPepCharPos[idxChar] = 0;
            lastPepStr.setCharAt(idxChar, Impl.Consts.alphabet.charAt(0));
        }
        
        lastPepCharPos[idxChar] += 1;
        lastPepStr.setCharAt(idxChar, Impl.Consts.alphabet.charAt(lastPepCharPos[idxChar]));
        
        return lastPepStr.toString();
    }
    
    private String GenerateFirst(long idxPep)
    {
        lastIdxPep = idxPep;
        
        //index of leftmost char to be incremented, starting from right ('least significant char')
        int idxAm = (int)Math.floor(Math.log(idxPep)/Math.log(Impl.Consts.aLen));  
        long idxPepRemainder = idxPep;
        int charPosInAlphabet = 0;
        for (; idxAm>=0; idxAm -= 1)
        {
            long charPosPow = (long)Math.pow(Impl.Consts.aLen, idxAm);
            charPosInAlphabet = (int)(idxPepRemainder / charPosPow);
            lastPepCharPos[pepLen-idxAm-1] = charPosInAlphabet;
            idxPepRemainder -= (charPosInAlphabet * charPosPow);
        }
        
        lastPepStr = ConvertCharcodesToString(lastPepCharPos);
        return lastPepStr.toString(); 
    }

    private static StringBuilder ConvertCharcodesToString(int[] charcodes)
    {
        StringBuilder pepConverter = new StringBuilder("AAAAAAAAA");
        for (int i=0; i<charcodes.length; ++i)
        {
            pepConverter.setCharAt(i, Impl.Consts.alphabet.charAt(charcodes[i]));
        }
        return pepConverter;
    }
    
    public static void main(String[] args) throws Exception 
    {
        try
        {
            if (args.length != 2)
            {
                throw new RuntimeException("Usage : java org.PSSMHC.PeptideGen index count");
            }

            long start = Long.parseLong(args[0]);
            long count = Long.parseLong(args[1]);
            PeptideGen pepGen = new PeptideGen();
            for (long i = start; i<start+count; ++i)
            {
                System.out.format(">\n%s\n", pepGen.Generate(i));
            }
        }
        catch (Exception ex)
        {
            System.err.print(ex + "\n");
        }
    }         
}
