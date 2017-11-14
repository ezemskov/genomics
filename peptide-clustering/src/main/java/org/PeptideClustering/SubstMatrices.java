package org.PeptideClustering;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import org.PSSMHC.Impl;

class SubstMatrixRow extends HashMap<Character, Double>         implements Serializable
{
    SubstMatrixRow() { super(SubstMatrix.matrixAlphabet.length(), 1.0f); }
}

class SubstMatrix    extends HashMap<Character, SubstMatrixRow> implements Serializable 
{
    //order of amino acids in matrices published by NCBI
    public static final String matrixAlphabet = "ARNDCQEGHILKMFPSTWYVBZX";
    
    public SubstMatrix(double[][] vals)
    {
        super(Impl.Consts.aLen, 1.0f);
        
        final String ExcMsg =  "Wrong substitution matrix size";
        if (vals.length != matrixAlphabet.length()) 
        {
            throw new RuntimeException(ExcMsg);
        }

        for (int i=0; i<vals.length; ++i) 
        {
            if (vals[i].length != matrixAlphabet.length()) 
            {
                throw new RuntimeException(ExcMsg);
            }
            
            SubstMatrixRow row = new SubstMatrixRow();
            for (int j=0; j<vals.length; ++j) 
            {
                row.put(matrixAlphabet.charAt(j), vals[i][j]);
            }
            
            put(matrixAlphabet.charAt(i), row);
        }
    }
}

public class SubstMatrices
{
    static final protected HashMap<String, SubstMatrix> matrices = new HashMap<String, SubstMatrix>();

    public static SubstMatrix get(String name) throws Exception
    {
        SubstMatrix res = matrices.get(name.toLowerCase());
        if (res == null)
        {
            throw new Exception("Valid matrix names are " + matrices.keySet().toString());
        }
        return res;
    }

    static 
    { 
        matrices.put("blosum50", new Blosum50()); 
        matrices.put("blosum62", new Blosum62()); 
        matrices.put("blosum80", new Blosum80()); 
        matrices.put("pam120", new Pam120()); 
        //matrices.put("pam150", new Pam150()); 
        //matrices.put("pam200", new Pam200()); 
    }    

    static class Blosum50 extends SubstMatrix implements Serializable
    {
        public Blosum50() { super(vals); }
        
        static double[][] vals = new double[][] {
            { 5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1},
            {-2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1},
            {-1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  4,  0, -1},
            {-2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  5,  1, -1},
            {-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -2},
            {-1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1},
            {-1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1},
            { 0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -2},
            {-2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1},
            {-1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1},
            {-2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1},
            {-1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1},
            {-1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1},
            {-3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -2},
            {-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -2},
            { 1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1},
            { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1,  0},
            {-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -3},
            {-2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1},
            { 0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -4, -3, -1},
            {-2, -1,  4,  5, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -4,  5,  2, -1},
            {-1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  2,  5, -1},
            {-1, -1, -1, -1, -2, -1, -1, -2, -1, -1, -1, -1, -1, -2, -2, -1,  0, -3, -1, -1, -1, -1, -1}
        };
    }
            
    static class Blosum62 extends SubstMatrix implements Serializable
    {
        public Blosum62() { super(vals); }

        static double[][] vals = new double[][] {
            { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0},
            {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1},
            {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1},
            {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1},
            { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2},
            {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1},
            {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1},
            { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1},
            {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1},
            {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1},
            {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1},
            {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1},
            {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1},
            {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1},
            {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2},
            { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0},
            { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0},
            {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2},
            {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1},
            { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1},
            {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1},
            {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1},
            { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1}        
        };
    }

    static class Blosum80 extends SubstMatrix implements Serializable
    {
        public Blosum80() { super(vals); }

        static double[][] vals = new double[][] {
            { 7, -3, -3, -3, -1, -2, -2,  0, -3, -3, -3, -1, -2, -4, -1,  2,  0, -5, -4, -1, -3, -2, -1},
            {-3,  9, -1, -3, -6,  1, -1, -4,  0, -5, -4,  3, -3, -5, -3, -2, -2, -5, -4, -4, -2,  0, -2},
            {-3, -1,  9,  2, -5,  0, -1, -1,  1, -6, -6,  0, -4, -6, -4,  1,  0, -7, -4, -5,  5, -1, -2},
            {-3, -3,  2, 10, -7, -1,  2, -3, -2, -7, -7, -2, -6, -6, -3, -1, -2, -8, -6, -6,  6,  1, -3},
            {-1, -6, -5, -7, 13, -5, -7, -6, -7, -2, -3, -6, -3, -4, -6, -2, -2, -5, -5, -2, -6, -7, -4},
            {-2,  1,  0, -1, -5,  9,  3, -4,  1, -5, -4,  2, -1, -5, -3, -1, -1, -4, -3, -4, -1,  5, -2},
            {-2, -1, -1,  2, -7,  3,  8, -4,  0, -6, -6,  1, -4, -6, -2, -1, -2, -6, -5, -4,  1,  6, -2},
            { 0, -4, -1, -3, -6, -4, -4,  9, -4, -7, -7, -3, -5, -6, -5, -1, -3, -6, -6, -6, -2, -4, -3},
            {-3,  0,  1, -2, -7,  1,  0, -4, 12, -6, -5, -1, -4, -2, -4, -2, -3, -4,  3, -5, -1,  0, -2},
            {-3, -5, -6, -7, -2, -5, -6, -7, -6,  7,  2, -5,  2, -1, -5, -4, -2, -5, -3,  4, -6, -6, -2},
            {-3, -4, -6, -7, -3, -4, -6, -7, -5,  2,  6, -4,  3,  0, -5, -4, -3, -4, -2,  1, -7, -5, -2},
            {-1,  3,  0, -2, -6,  2,  1, -3, -1, -5, -4,  8, -3, -5, -2, -1, -1, -6, -4, -4, -1,  1, -2},
            {-2, -3, -4, -6, -3, -1, -4, -5, -4,  2,  3, -3,  9,  0, -4, -3, -1, -3, -3,  1, -5, -3, -2},
            {-4, -5, -6, -6, -4, -5, -6, -6, -2, -1,  0, -5,  0, 10, -6, -4, -4,  0,  4, -2, -6, -6, -3},
            {-1, -3, -4, -3, -6, -3, -2, -5, -4, -5, -5, -2, -4, -6, 12, -2, -3, -7, -6, -4, -4, -2, -3},
            { 2, -2,  1, -1, -2, -1, -1, -1, -2, -4, -4, -1, -3, -4, -2,  7,  2, -6, -3, -3,  0, -1, -1},
            { 0, -2,  0, -2, -2, -1, -2, -3, -3, -2, -3, -1, -1, -4, -3,  2,  8, -5, -3,  0, -1, -2, -1},
            {-5, -5, -7, -8, -5, -4, -6, -6, -4, -5, -4, -6, -3,  0, -7, -6, -5, 16,  3, -5, -8, -5, -5},
            {-4, -4, -4, -6, -5, -3, -5, -6,  3, -3, -2, -4, -3,  4, -6, -3, -3,  3, 11, -3, -5, -4, -3},
            {-1, -4, -5, -6, -2, -4, -4, -6, -5,  4,  1, -4,  1, -2, -4, -3,  0, -5, -3,  7, -6, -4, -2},
            {-3, -2,  5,  6, -6, -1,  1, -2, -1, -6, -7, -1, -5, -6, -4,  0, -1, -8, -5, -6,  6,  0, -3},
            {-2,  0, -1,  1, -7,  5,  6, -4,  0, -6, -5,  1, -3, -6, -2, -1, -2, -5, -4, -4,  0,  6, -1},
            {-1, -2, -2, -3, -4, -2, -2, -3, -2, -2, -2, -2, -2, -3, -3, -1, -1, -5, -3, -2, -3, -1, -2}
        };
    }

    static class Pam120 extends SubstMatrix implements Serializable
    {
        public Pam120() { super(vals); }

        static double[][] vals = new double[][] {
            { 2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0},
            {-2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1},
            { 0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0},
            { 0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1},
            {-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3},
            { 0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1},
            { 0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1},
            { 1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1},
            {-1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1},
            {-1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1},
            {-2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1},
            {-1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1},
            {-1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1},
            {-3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2},
            { 1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1},
            { 1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0},
            { 1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0},
            {-6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4},
            {-3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2},
            { 0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1},
            { 0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1},
            { 0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1},
            { 0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1}
        };
    }

    static class Pam150 extends SubstMatrix implements Serializable
    {
        public Pam150() { super(vals); }

        static double[][] vals = new double[][] {
            { 3, -2,  0,  0, -2, -1,  0,  1, -2, -1, -2, -2, -1, -4,  1,  1,  1, -6, -3,  0,  0,  0, -1},
            {-2,  6, -1, -2, -4,  1, -2, -3,  1, -2, -3,  3, -1, -4, -1, -1, -2,  1, -4, -3, -2,  0, -1},
            { 0, -1,  3,  2, -4,  0,  1,  0,  2, -2, -3,  1, -2, -4, -1,  1,  0, -4, -2, -2,  3,  1, -1},
            { 0, -2,  2,  4, -6,  1,  3,  0,  0, -3, -5, -1, -3, -6, -2,  0, -1, -7, -4, -3,  3,  2, -1},
            {-2, -4, -4, -6,  9, -6, -6, -4, -3, -2, -6, -6, -5, -5, -3,  0, -3, -7,  0, -2, -5, -6, -3},
            {-1,  1,  0,  1, -6,  5,  2, -2,  3, -3, -2,  0, -1, -5,  0, -1, -1, -5, -4, -2,  1,  4, -1},
            { 0, -2,  1,  3, -6,  2,  4, -1,  0, -2, -4, -1, -2, -6, -1, -1, -1, -7, -4, -2,  2,  4, -1},
            { 1, -3,  0,  0, -4, -2, -1,  4, -3, -3, -4, -2, -3, -5, -1,  1, -1, -7, -5, -2,  0, -1, -1},
            {-2,  1,  2,  0, -3,  3,  0, -3,  6, -3, -2, -1, -3, -2, -1, -1, -2, -3,  0, -3,  1,  1, -1},
            {-1, -2, -2, -3, -2, -3, -2, -3, -3,  5,  1, -2,  2,  0, -3, -2,  0, -5, -2,  3, -2, -2, -1},
            {-2, -3, -3, -5, -6, -2, -4, -4, -2,  1,  5, -3,  3,  1, -3, -3, -2, -2, -2,  1, -4, -3, -2},
            {-2,  3,  1, -1, -6,  0, -1, -2, -1, -2, -3,  4,  0, -6, -2, -1,  0, -4, -4, -3,  0,  0, -1},
            {-1, -1, -2, -3, -5, -1, -2, -3, -3,  2,  3,  0,  7, -1, -3, -2, -1, -5, -3,  1, -3, -2, -1},
            {-4, -4, -4, -6, -5, -5, -6, -5, -2,  0,  1, -6, -1,  7, -5, -3, -3, -1,  5, -2, -5, -5, -3},
            { 1, -1, -1, -2, -3,  0, -1, -1, -1, -3, -3, -2, -3, -5,  6,  1,  0, -6, -5, -2, -2, -1, -1},
            { 1, -1,  1,  0,  0, -1, -1,  1, -1, -2, -3, -1, -2, -3,  1,  2,  1, -2, -3, -1,  0, -1,  0},
            { 1, -2,  0, -1, -3, -1, -1, -1, -2,  0, -2,  0, -1, -3,  0,  1,  4, -5, -3,  0,  0, -1, -1},
            {-6,  1, -4, -7, -7, -5, -7, -7, -3, -5, -2, -4, -5, -1, -6, -2, -5, 12, -1, -6, -5, -6, -4},
            {-3, -4, -2, -4,  0, -4, -4, -5,  0, -2, -2, -4, -3,  5, -5, -3, -3, -1,  8, -3, -3, -4, -3},
            { 0, -3, -2, -3, -2, -2, -2, -2, -3,  3,  1, -3,  1, -2, -2, -1,  0, -6, -3,  4, -2, -2, -1},
            { 0, -2,  3,  3, -5,  1,  2,  0,  1, -2, -4,  0, -3, -5, -2,  0,  0, -5, -3, -2,  3,  2, -1},
            { 0,  0,  1,  2, -6,  4,  4, -1,  1, -2, -3,  0, -2, -5, -1, -1, -1, -6, -4, -2,  2,  4, -1},
            {-1, -1, -1, -1, -3, -1, -1, -1, -1, -1, -2, -1, -1, -3, -1,  0, -1, -4, -3, -1, -1, -1, -1}
        };
    }

    static class Pam250 extends SubstMatrix implements Serializable
    {
        public Pam250() { super(vals); }

        static double[][] vals = new double[][] {
            { 2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0},
            {-2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1},
            { 0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0},
            { 0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1},
            {-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3},
            { 0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1},
            { 0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1},
            { 1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1},
            {-1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1},
            {-1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1},
            {-2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1},
            {-1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1},
            {-1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1},
            {-3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2},
            { 1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1},
            { 1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0},
            { 1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0},
            {-6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4},
            {-3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2},
            { 0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1},
            { 0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1},
            { 0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1},
            { 0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1}
        };
    }
}
,
,
,
,
,
,
,
,
,
,
,
,
,
,
,
,
,
,
,
,