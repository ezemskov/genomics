package org.PSSMHC;

import java.lang.StringBuilder;

public class PeptideGen
{
    public static int PepLenDefault = 9;
    private int _pepLen = -1;
    private long _maxIndex = -1; 
    
    private StringBuilder _lastPepStr = null;
    private int[]  _lastPepCharPos = null;
    private long   _lastIdxPep = -2;
    
    public PeptideGen(int pepLen)
    {
        _pepLen = pepLen;
        _maxIndex = (long)Math.pow(Impl.Consts.aLen, _pepLen) - 1; //20^N - 1
        _lastPepCharPos = new int[_pepLen];
    }

    public PeptideGen() 
    { 
        this(PepLenDefault); 
    }

    /**
     * @param index number of peptide in lexicographical order, 
     * i.e. AAAAAAAAA is #0, AAAAAAAAC is #1, and YYYYYYYYY is #20^9-1
     * @return peptide 
     */
    public String Generate(long idxPep)
    {
        if ((idxPep < 0) || (idxPep > _maxIndex))
        {
            throw new IndexOutOfBoundsException(String.format("Peptide index %d is out of bounds", idxPep));
        }

        return (idxPep == _lastIdxPep+1) ? GenerateNext() : GenerateFirst(idxPep);
    }
    
    private String GenerateNext()
    {
        _lastIdxPep += 1;
        
        int idxChar = _pepLen-1;
        //find left bound of '....YYY' chars sequence, starting from right, 
        //and set it to '....AAAA' (i.e. perform rollover)
        for (; (idxChar >= 0) && (_lastPepCharPos[idxChar] == Impl.Consts.aLen-1); --idxChar)
        {
            _lastPepCharPos[idxChar] = 0;
            _lastPepStr.setCharAt(idxChar, Impl.Consts.alphabet.charAt(0));
        }
        
        _lastPepCharPos[idxChar] += 1;
        _lastPepStr.setCharAt(idxChar, Impl.Consts.alphabet.charAt(_lastPepCharPos[idxChar]));
        
        return _lastPepStr.toString();
    }
    
    private String GenerateFirst(long idxPep)
    {
        _lastIdxPep = idxPep;
        
        //index of leftmost char to be incremented, starting from right ('least significant char')
        int idxAm = (int)Math.floor(Math.log(idxPep)/Math.log(Impl.Consts.aLen));  
        long idxPepRemainder = idxPep;
        int charPosInAlphabet = 0;
        for (; idxAm>=0; idxAm -= 1)
        {
            long charPosPow = (long)Math.pow(Impl.Consts.aLen, idxAm);
            charPosInAlphabet = (int)(idxPepRemainder / charPosPow);
            _lastPepCharPos[_pepLen-idxAm-1] = charPosInAlphabet;
            idxPepRemainder -= (charPosInAlphabet * charPosPow);
        }
        
        _lastPepStr = ConvertCharcodesToString(_lastPepCharPos);
        return _lastPepStr.toString(); 
    }

    private static StringBuilder ConvertCharcodesToString(int[] charcodes)
    {
        StringBuilder pepConverter = new StringBuilder(charcodes.length);
        for (int i=0; i<charcodes.length; ++i)
        {
            pepConverter.append(Impl.Consts.alphabet.charAt(charcodes[i]));
        }
        return pepConverter;
    }
}
