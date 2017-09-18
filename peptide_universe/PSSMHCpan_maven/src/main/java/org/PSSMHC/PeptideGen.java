package org.PSSMHC;

public class PeptideGen
{
    private static int pepLen = 9;
    private static double lnALen = Math.log(Consts.aLen);         //ln(20)
    private static long maxIndex = (long)Math.pow(Consts.aLen, pepLen);  //20^9      
    
    /**
     * @param index number of peptide in lexicographical order, 
     * i.e. AAAAAAAAA is #0, AAAAAAAAC is #1, and YYYYYYYYY is #20^9-1
     * @return peptide 
     */
    public static String Generate(long idxPep)
    {
        if ((idxPep < 0) || (idxPep > maxIndex))
        {
            throw new RuntimeException(String.format("Peptide index %d is out of bounds", idxPep));
        }
        
        StringBuffer res = new StringBuffer("AAAAAAAAA");
        
        //index of leftmost amino to be incremented, starting from right ('least significant amino')
        int idxAmMax = (int)Math.floor(Math.log(idxPep)/lnALen);  
        //for (int idxAm=idxAmMax; idxAm < pepLen; idxAm += 1)
        long idxPepRemainder = idxPep;
        for (int idxAm = idxAmMax; idxAm>=0; idxAm -= 1)
        {
            long charPosPow = (long)Math.pow(Consts.aLen, idxAm); //todo : integer powers here !
            int charPosInAlphabet = (int)(idxPepRemainder / charPosPow);
            res.setCharAt(pepLen-idxAm-1, Consts.alphabet.charAt(charPosInAlphabet));
            idxPepRemainder -= (charPosInAlphabet * charPosPow);
        }
        
        return res.toString();
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
            for (long i = start; i<start+count; ++i)
            {
                System.out.format(">\n%s\n", PeptideGen.Generate(i));
            }
        }
        catch (Exception ex)
        {
            System.err.print(ex + "\n");
        }
    }         
}
