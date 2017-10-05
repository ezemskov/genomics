package org.PSSMHC;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;

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

    final class PeptideBloomTest
{
    private static ArrayList<String> ParsePeptideResFile(String filename) throws FileNotFoundException, IOException
    {
        ArrayList<String> res = new ArrayList<String>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        while (reader.ready())
        {
            String rowStr = reader.readLine();
            res.add(rowStr.split(",")[0]);
        }
        
        return res;
    }
        
    public static void main(String[] args) 
    {
        try
        {
            if (args.length < 1)
            {
                throw new RuntimeException("Usage : java org.PSSMHC.PeptideBloomTest bloom_filter peptides_list");
            }

            PeptideBloomFilter bloom = new PeptideBloomFilter(args[0]);
            
            int count = 0;
            ArrayList<String> peptides = ParsePeptideResFile(args[1]);
            for (String p : peptides)
            {
                if (bloom.MightContain(p))
                {
                    count += 1;
                }
            }
            System.out.format("Bloom filter returns true for %d of %d peptides\n", count, peptides.size());
        }
        catch (Exception e) 
        {
            System.err.format("Error reading %s\n", e.getMessage());
        }        
    }
}
