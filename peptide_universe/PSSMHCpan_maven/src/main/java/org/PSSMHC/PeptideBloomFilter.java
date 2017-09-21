package org.PSSMHC;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.Exception;

import org.apache.spark.util.sketch.BloomFilter;

import java.nio.charset.Charset;

class PeptideBloomFilter
{
    private BloomFilter filter = null;
    private static long size = (long)1E8;
    private static double fpp = 0.01;
    
    PeptideBloomFilter()
    {
        filter = BloomFilter.create(size, fpp);
    }

    PeptideBloomFilter(String filename)
    {
        try
        {
            FileInputStream fs = new FileInputStream(filename);
            ObjectInputStream objs = new ObjectInputStream(fs);
            filter = (BloomFilter)objs.readObject();
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
        filter.putString(peptide);
    }
    
    public boolean MightContain(String peptide)
    {
        return filter.mightContainString(peptide);
    }
}
