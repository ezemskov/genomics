package org.PSSMHC;

import java.io.Serializable;

public class ScoredPeptide implements Serializable
{
    ScoredPeptide(String peptide, double ic50) 
    { 
        this.peptide = peptide; 
        this.ic50 = ic50;
    }
    
    public String toString()
    {
        return String.format("%s,%.0f", peptide, ic50);
    }
        
    public String peptide; 
    public double ic50; 
}


