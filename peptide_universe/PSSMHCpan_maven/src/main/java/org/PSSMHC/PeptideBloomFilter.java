package org.PSSMHC;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.Exception;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnels;
import com.google.common.hash.Funnel;
import java.nio.charset.Charset;

class PeptideBloomFilter
{
    private BloomFilter<CharSequence> filter = null;
    private static Funnel<CharSequence> funnel = Funnels.stringFunnel(Charset.defaultCharset()); //todo : 20-char charset ?
    private static int size = Integer.MAX_VALUE;
    private static double fpp = 0.01;
    
    PeptideBloomFilter()
    {
        filter = BloomFilter.create(funnel, size, fpp);
    }

    PeptideBloomFilter(String filename)
    {
        try
        {
            FileInputStream fs = new FileInputStream(filename);
            ObjectInputStream objs = new ObjectInputStream(fs);
            filter = (BloomFilter<CharSequence>)objs.readObject();
            fs.close();
        }
        catch (Exception e) 
        {
            System.err.format("Error reading Bloom filter from %s : %s", filename, e);
        }
    }
    
    public boolean Save(String filename)
    {
        try
        {
            FileOutputStream fs = new FileOutputStream(filename);
            ObjectOutputStream objs = new ObjectOutputStream(fs);
            objs.writeObject(filter);
            fs.close();
            return true;
        }
        catch (Exception e) 
        {
            System.err.format("Error writing to %s : %s", filename, e);
            return false;
        }        
    }

    public void Add(String peptide)
    {
        filter.put(peptide);
    }
    
    public boolean MightContain(String peptide)
    {
        return filter.mightContain(peptide);
    }
}
