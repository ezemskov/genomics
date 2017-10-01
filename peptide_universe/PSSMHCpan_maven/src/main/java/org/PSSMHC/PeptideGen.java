package org.PSSMHC;

import java.lang.StringBuilder;

public class PeptideGen
{
    private static int pepLen = 9;
    private static double lnALen = Math.log(Consts.aLen);         //ln(20)
    private static long maxIndex = (long)Math.pow(Consts.aLen, pepLen);  //20^9      
    
    private int[]  lastPepCharPos = new int[pepLen];
    private long   lastIdxPep = -2;
    
    private StringBuilder pepConverter = new StringBuilder("AAAAAAAAA");

    /**
     * @param index number of peptide in lexicographical order, 
     * i.e. AAAAAAAAA is #0, AAAAAAAAC is #1, and YYYYYYYYY is #20^9-1
     * @return peptide 
     */
    public String Generate(long idxPep)
    {
        //return (idxPep == lastIdxPep+1) ? GenerateNext() : GenerateFirst(idxPep);
        return GenerateFirst(idxPep);
    }
    
    private int FindLastIdxOfLastAmino()
    {
        int i = lastPepCharPos.length-1;
        for (; (i >= 0) && (lastPepCharPos[i] == Consts.aLen-1); --i) {}
        return i;
    }
    
    private String ConvertLastCharPosToString()
    {
        for (int i=0; i<lastPepCharPos.length; ++i)
        {
            pepConverter.setCharAt(i, Consts.alphabet.charAt(lastPepCharPos[i]));
        }
        return pepConverter.toString();
    }
    
    private String GenerateNext()
    {
        lastIdxPep += 1;
        
        if (lastPepCharPos[pepLen-1] == Consts.aLen-1)
        {
            int idx=FindLastIdxOfLastAmino();
            for (int i=idx+1; i<lastPepCharPos.length; ++i)
            {
                lastPepCharPos[i] = 0;
            }
            lastPepCharPos[idx] += 1; 
        }
        else
        {
            lastPepCharPos[pepLen-1] += 1;
        }
        
        return ConvertLastCharPosToString();
    }
    
    private String GenerateFirst(long idxPep)
    {
        lastIdxPep = idxPep;
        if ((idxPep < 0) || (idxPep > maxIndex))
        {
            throw new RuntimeException(String.format("Peptide index %d is out of bounds", idxPep));
        }
        
        //index of leftmost amino to be incremented, starting from right ('least significant amino')
        int idxAm = (int)Math.floor(Math.log(idxPep)/lnALen);  
        long idxPepRemainder = idxPep;
        int charPosInAlphabet = 0;
        for (; idxAm>=0; idxAm -= 1)
        {
            long charPosPow = (long)Math.pow(Consts.aLen, idxAm); //todo : integer powers here !
            charPosInAlphabet = (int)(idxPepRemainder / charPosPow);
            lastPepCharPos[pepLen-idxAm-1] = charPosInAlphabet;
            idxPepRemainder -= (charPosInAlphabet * charPosPow);
        }
        
        return ConvertLastCharPosToString();
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
