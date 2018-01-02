//--------------------------------------------------------------------------------------------------
// 
// Description: Implementation of Needleman-Wunsch global alignment.
// This code is written by Ren-Xiang Yan in China Agricultural University and is originally based on 
// the fortran implementation from Dr. Yang Zhang (http://zhanglab.ccmb.med.umich.edu/NW-align/).
// Last update is in 2010/08/14. 
//-----------------x-------------------x-------------------x--------------------------------------------

package org.PeptideClustering;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.PSSMHC.Impl;
import org.PSSMHC.Xml;
import static org.PeptideClustering.AssignBindersToClusters.FormatMap;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.Encoders;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.storage.StorageLevel;
import scala.Tuple2;

public class NeedlemanWunsch 
{
    static final int gap_open = -11, gap_extn = -1; 
    
    static int[][] imut = new int[24][24];
    
    public static class AlignResult implements Serializable
    {
        public double identity = Double.NaN;
        public int score = 0;
        public int alignedLength = 0;
        public int identicalLength = 0;
        
        @Override
        public String toString()
        {
            return String.format("score %d alignLen %d identLen %d ident %.2f", 
                score, alignedLength, identicalLength, identity);
        }
    }
    
    public static AlignResult calcNWscore(String f1,String f2)  
    {
        String seqW = "*ARNDCQEGHILKMFPSTWYVBZX";	       // Amino acide order in the BLAST's scoring matrix (e.g.,Blosum62). 
        f1 = "*" + f1;                                     // Add a '*' character in the head of a sequence and this can make java code much more consistent with orginal fortran code.   
        f2 = "*" + f2;                                     // Use 1 to represent the first position of the sequence in the original fortran code,and 1 stand for the second position in java code. Here, add a '*' character in the head of a sequence could make 1 standard for the first postion of thse sequence in java code.     
        int[] seq1 = new int[f1.length()];                 
        int[] seq2 = new int[f2.length()];		           // seq1 and seq2 are arrays that store the amino acid order numbers of sequence1 and sequence2.
        int i,j;		                                   // For example, 1 stand for A, 2 represent R and etc.
        for(i=1;i<f1.length();i++)
        {
                for(j=1;j<seqW.length();j++)
                {
                        if(f1.charAt(i)==seqW.charAt(j))
                        {
                                seq1[i]=j;
                        }
                }
        }		

        for(i=1;i<f2.length();i++)
        {
                for(j=1;j<seqW.length();j++)
                {
                        if(f2.charAt(i)==seqW.charAt(j))
                        {
                                seq2[i]=j;
                        }
                }
        }

         int[][] score = new int[f1.length()][f2.length()];		// score[i][j] stard for the alignment score that align ith position of the first sequence to the jth position of the second sequence.
         for(i=1;i<f1.length();i++)
         {
                  for(j=1;j<f2.length();j++)
                  {
                          score[i][j] = imut[seq1[i]][seq2[j]];
                  }
         }	

        int[] j2i = new int[f2.length()+1];
        for(j=1;j<f2.length();j++)
        {
                j2i[j] = -1; // !all are not aligned
        }		



        int[][] val = new int[f1.length()+1][f2.length()+1];  // val[][] was assigned as a global variable, and the value could be printed in the final.
        int[][] idir = new int[f1.length()+1][f2.length()+1];
        int[][] preV = new int[f1.length()+1][f2.length()+1];
        int[][] preH = new int[f1.length()+1][f2.length()+1];
        int D,V,H;
        boolean standard = true;
        if(standard)    // If you want to use alternative implementation of Needleman-Wunsch dynamic program , you can assign "false" value to the "standard" variable.  
        {			
                ////////////////////////////////////////////////////////////////////////////////
                //	     This is a standard Needleman-Wunsch dynamic program (by Y. Zhang 2005).
                //	     1. Count multiple-gap.
                //	     2. The gap penality W(k)=Go+Ge*k1+Go+Ge*k2 if gap open on both sequences			
                //	     idir[i][j]=1,2,3, from diagonal, horizontal, vertical
                //	     val[i][j] is the cumulative score of (i,j)
                ////////////////////////////////////////////////////////////////////////////////			

                int[][] jpV = new int[f1.length()+1][f2.length()+1];
                int[][] jpH = new int[f1.length()+1][f2.length()+1];								
                val[0][0]=0;
                val[1][0] =gap_open;	
                for(i=2;i<f1.length();i++){
                        val[i][0] = val[i-1][0]+gap_extn;
                }
                for(i=1;i<f1.length();i++)
                {

                        preV[i][0] = val[i][0]; // not use preV at the beginning
                idir[i][0] = 0;         // useless
                jpV[i][0] = 1;          // useless
                jpH[i][0] = i;          // useless
                }
                val[0][1]=gap_open;
                for(j=2;j<f2.length();j++){
                        val[0][j]=val[0][j-1]+gap_extn;
                }
                for(j=1;j<f2.length();j++)
                {
                 preH[0][j] = val[0][j];
                 idir[0][j] = 0;
                 jpV[0][j] = j;
                 jpH[0][j] = 1;
                }

                // DP ------------------------------>
                for(j=1;j<f2.length();j++)
                {			
                        for(i=1;i<f1.length();i++)
                        {			
                                // D=VAL(i-1,j-1)+SCORE(i,j)--------------->
                                D = val[i-1][j-1] + score[i][j];	// from diagonal, val(i,j) is val(i-1,j-1)			

                                //	H=H+gap_open ------->				
                                jpH[i][j] = 1;				
                                int val1 = val[i-1][j] + gap_open;  // gap_open from both D and V
                                int val2 = preH[i-1][j] + gap_extn; // gap_extn from horizontal
                                if(val1>val2)   // last step from D or V
                                {
                                        H = val1;
                                }
                                else            // last step from H
                                {
                                        H = val2;
                                        if(i > 1)
                                        {
                                                jpH[i][j] = jpH[i-1][j] + 1;  // record long-gap
                                        }
                                }

                                //	V=V+gap_open --------->					
                                jpV[i][j] = 1;
                                val1 = val[i][j-1] + gap_open;
                                val2 = preV[i][j-1] + gap_extn;
                                if(val1>val2)
                                {
                                        V = val1;
                                }
                                else
                                {
                                        V = val2;
                                        if(j > 1)
                                        {
                                                jpV[i][j] = jpV[i][j-1] + 1;   // record long-gap   
                                        }
                                }

                                preH[i][j] = H; // unaccepted H
                                preV[i][j] = V;	// unaccepted V				
                                if((D>H)&&(D>V))
                                {
                                        idir[i][j]=1;
                                        val[i][j]=D;
                                }
                                else if(H > V)
                                {   
                                        idir[i][j] = 2;
                        val[i][j] = H;
                                }
                    else
                    {
                         idir[i][j] = 3;
                              val[i][j] = V;
                    }
                        }			
                }			

                //  tracing back the pathway			
                i = f1.length()-1;
                j = f2.length()-1;		
                while((i>0)&&(j>0))  
                {			 
                         if(idir[i][j]==1)       // from diagonal
                         {
                                j2i[j] = i;
                    i=i-1;
                    j=j-1;
                         }
                         else if(idir[i][j]==2)  // from horizonal
                 { 	        	 
                                 int temp1 = jpH[i][j];                  //  
                         for(int me=1;me<=temp1;me++)            //  In the point view of a programer, 
                         {                                       //  you should not use the  "for(int me=1;me<=jpH[i][j];me++)".
                                if(i>0)	        	                 //  If you use up sentence,the value of jpH[i][j] is changed when variable i changes. 
                        {                                    //  So the value of jpH[i][j] was assigned to the value temp1 and use the setence "for(int me=1;me<=temp1;me++)" here. 
                           i=i-1;                            // 
                        }	                                 //
                          }	        	                                
                 }
                         else
                 {	 
                                int temp2 = jpV[i][j]; 
                    for(int me=1;me<=temp2;me++)             //  In the point view of a programer,
                    {                                        //  you should not use the  "for(int me=1;me<=jpV[i][j];me++)".
                       if(j>0)                               //  Because when variable i change, the jpV[i][j] employed here is also change. 
                       {                                     //  So the value of jpV[i][j] was assigned to the value temp2 and use the setence "for(int me=1;me<=temp2;me++)" here.
                          j=j-1;                             //
                       }
                    }	           
                 }			 
                }	
        }
        else   
        {			
                /////////////////////////////////////////////////////////////////////////////////
                //		     This is an alternative implementation of Needleman-Wunsch dynamic program 
                //		     (by Y. Zhang 2005)
                //		     1. Count two-layer iteration and multiple-gaps
                //		     2. The gap penality W(k)=Go+Ge*k1+Ge*k2 if gap open on both sequences
                //		
                //		     idir[i][j]=1,2,3, from diagonal, horizontal, vertical
                //		     val[i][j] is the cumulative score of (i,j)
                ////////////////////////////////////////////////////////////////////////////////			

                int[][] preD = new int[f1.length()+1][f2.length()+1];
                int[][] idirH = new int[f1.length()+1][f2.length()+1];
                int[][] idirV = new int[f1.length()+1][f2.length()+1];						
                val[0][0] = 0;		
                for(i=1;i<f1.length();i++)
                {
                        val[i][0] = 0;
                        idir[i][0] = 0;
                        preD[i][0] = 0;
                        preH[i][0] = -1000;
                        preV[i][0] = -1000;
                }

                for(j=1;j<f2.length();j++)
                {
                        val[0][j] = 0;
                        idir[0][j] = 0;
                        preD[0][j] = 0;
                        preH[0][j] = -1000;
                        preV[0][j] = -1000;
                }

                // DP ------------------------------>
                for(j=1;j<f2.length();j++)
                {			
                        for(i=1;i<f1.length();i++)
                        {			
                                // preD=VAL(i-1,j-1)+SCORE(i,j)--------------->
                                preD[i][j] = val[i-1][j-1] + score[i][j];
                                // preH: pre-accepted H----------------------->
                                D = preD[i-1][j] + gap_open;
                                H = preH[i-1][j] + gap_extn;
                                V = preV[i-1][j] + gap_extn;
                                if((D>H)&&(D>V))
                                {
                                        preH[i][j] = D;
                                        idirH[i-1][j] = 1;
                                }
                                else if(H>V)
                                {
                                        preH[i][j] = H;
                                        idirH[i-1][j] = 2;
                                }
                                else
                                {
                                        preH[i][j] = V;
                                        idirH[i-1][j] = 3;
                                }

                                // preV: pre-accepted V----------------------->
                                D = preD[i][j-1] + gap_open;
                                H = preH[i][j-1] + gap_extn;
                                V = preV[i][j-1] + gap_extn;				
                                if((D>H)&&(D>V))
                                {
                                        preV[i][j] = D;
                                        idirV[i][j-1] = 1;
                                }
                                else if(H>V)
                                {
                                        preV[i][j] = H;
                                        idirV[i][j-1] = 2;
                                }
                                else
                                {
                                        preV[i][j] = V;
                                        idirV[i][j-1] = 3;
                                }

                                // decide idir(i,j)----------->
                                if((preD[i][j]>preH[i][j])&&(preD[i][j]>preV[i][j]))
                                {
                                        idir[i][j] = 1;
                                        val[i][j] = preD[i][j];
                                }			
                                else if(preH[i][j]>preV[i][j])
                                {
                                        idir[i][j] = 2;
                                        val[i][j] = preH[i][j];
                                }
                                else
                                {
                                        idir[i][j] = 3;
                                        val[i][j] = preV[i][j];
                                }				
                        }			
                }		

                //  tracing back the pathway			
                i = f1.length()-1;
                j = f2.length()-1;		
                while((i>0)&&(j>0))  
                {			 
                         if(idir[i][j]==1)
                         {
                                 j2i[j] = i;
                                 i = i - 1;
                                 j = j - 1;
                         }
                         else if(idir[i][j]==2)
                         {
                                 i = i - 1;
                                 idir[i][j] = idirH[i][j];
                         }
                         else
                         {
                                 j = j - 1;
                                 idir[i][j] = idirV[i][j];
                         }
                }			
        }		

        // calculate sequence identity			
        int L_id=0;
        int L_ali=0;
        for(j=1;j<f2.length();j++)
        {	    		
                    if(j2i[j]>0)
                {
                            i=j2i[j];
                            L_ali=L_ali+1;
                        if(seq1[i]==seq2[j])
                        {	            	
                            L_id=L_id+1;
                        }
                } 	        
        }	   
        
        AlignResult res = new AlignResult();
        res.alignedLength = L_ali;
        res.identicalLength = L_id;
        res.identity = L_id*1.0/(f2.length()-1);	    
        res.score = val[f1.length()-1][f2.length()-1];
        return res;
    }
    
    static  //BLOSUM62 used in BLAST.
    {
        imut[1][1]=4;
        imut[1][2]=-1;
        imut[1][3]=-2;
        imut[1][4]=-2;
        imut[1][5]=0;
        imut[1][6]=-1;
        imut[1][7]=-1;
        imut[1][8]=0;
        imut[1][9]=-2;
        imut[1][10]=-1;
        imut[1][11]=-1;
        imut[1][12]=-1;
        imut[1][13]=-1;
        imut[1][14]=-2;
        imut[1][15]=-1;
        imut[1][16]=1;
        imut[1][17]=0;
        imut[1][18]=-3;
        imut[1][19]=-2;
        imut[1][20]=0;
        imut[1][21]=-2;
        imut[1][22]=-1;
        imut[1][23]=0;
        imut[2][1]=-1;
        imut[2][2]=5;
        imut[2][3]=0;
        imut[2][4]=-2;
        imut[2][5]=-3;
        imut[2][6]=1;
        imut[2][7]=0;
        imut[2][8]=-2;
        imut[2][9]=0;
        imut[2][10]=-3;
        imut[2][11]=-2;
        imut[2][12]=2;
        imut[2][13]=-1;
        imut[2][14]=-3;
        imut[2][15]=-2;
        imut[2][16]=-1;
        imut[2][17]=-1;
        imut[2][18]=-3;
        imut[2][19]=-2;
        imut[2][20]=-3;
        imut[2][21]=-1;
        imut[2][22]=0;
        imut[2][23]=-1;
        imut[3][1]=-2;
        imut[3][2]=0;
        imut[3][3]=6;
        imut[3][4]=1;
        imut[3][5]=-3;
        imut[3][6]=0;
        imut[3][7]=0;
        imut[3][8]=0;
        imut[3][9]=1;
        imut[3][10]=-3;
        imut[3][11]=-3;
        imut[3][12]=0;
        imut[3][13]=-2;
        imut[3][14]=-3;
        imut[3][15]=-2;
        imut[3][16]=1;
        imut[3][17]=0;
        imut[3][18]=-4;
        imut[3][19]=-2;
        imut[3][20]=-3;
        imut[3][21]=3;
        imut[3][22]=0;
        imut[3][23]=-1;
        imut[4][1]=-2;
        imut[4][2]=-2;
        imut[4][3]=1;
        imut[4][4]=6;
        imut[4][5]=-3;
        imut[4][6]=0;
        imut[4][7]=2;
        imut[4][8]=-1;
        imut[4][9]=-1;
        imut[4][10]=-3;
        imut[4][11]=-4;
        imut[4][12]=-1;
        imut[4][13]=-3;
        imut[4][14]=-3;
        imut[4][15]=-1;
        imut[4][16]=0;
        imut[4][17]=-1;
        imut[4][18]=-4;
        imut[4][19]=-3;
        imut[4][20]=-3;
        imut[4][21]=4;
        imut[4][22]=1;
        imut[4][23]=-1;
        imut[5][1]=0;
        imut[5][2]=-3;
        imut[5][3]=-3;
        imut[5][4]=-3;
        imut[5][5]=9;
        imut[5][6]=-3;
        imut[5][7]=-4;
        imut[5][8]=-3;
        imut[5][9]=-3;
        imut[5][10]=-1;
        imut[5][11]=-1;
        imut[5][12]=-3;
        imut[5][13]=-1;
        imut[5][14]=-2;
        imut[5][15]=-3;
        imut[5][16]=-1;
        imut[5][17]=-1;
        imut[5][18]=-2;
        imut[5][19]=-2;
        imut[5][20]=-1;
        imut[5][21]=-3;
        imut[5][22]=-3;
        imut[5][23]=-2;
        imut[6][1]=-1;
        imut[6][2]=1;
        imut[6][3]=0;
        imut[6][4]=0;
        imut[6][5]=-3;
        imut[6][6]=5;
        imut[6][7]=2;
        imut[6][8]=-2;
        imut[6][9]=0;
        imut[6][10]=-3;
        imut[6][11]=-2;
        imut[6][12]=1;
        imut[6][13]=0;
        imut[6][14]=-3;
        imut[6][15]=-1;
        imut[6][16]=0;
        imut[6][17]=-1;
        imut[6][18]=-2;
        imut[6][19]=-1;
        imut[6][20]=-2;
        imut[6][21]=0;
        imut[6][22]=3;
        imut[6][23]=-1;
        imut[7][1]=-1;
        imut[7][2]=0;
        imut[7][3]=0;
        imut[7][4]=2;
        imut[7][5]=-4;
        imut[7][6]=2;
        imut[7][7]=5;
        imut[7][8]=-2;
        imut[7][9]=0;
        imut[7][10]=-3;
        imut[7][11]=-3;
        imut[7][12]=1;
        imut[7][13]=-2;
        imut[7][14]=-3;
        imut[7][15]=-1;
        imut[7][16]=0;
        imut[7][17]=-1;
        imut[7][18]=-3;
        imut[7][19]=-2;
        imut[7][20]=-2;
        imut[7][21]=1;
        imut[7][22]=4;
        imut[7][23]=-1;
        imut[8][1]=0;
        imut[8][2]=-2;
        imut[8][3]=0;
        imut[8][4]=-1;
        imut[8][5]=-3;
        imut[8][6]=-2;
        imut[8][7]=-2;
        imut[8][8]=6;
        imut[8][9]=-2;
        imut[8][10]=-4;
        imut[8][11]=-4;
        imut[8][12]=-2;
        imut[8][13]=-3;
        imut[8][14]=-3;
        imut[8][15]=-2;
        imut[8][16]=0;
        imut[8][17]=-2;
        imut[8][18]=-2;
        imut[8][19]=-3;
        imut[8][20]=-3;
        imut[8][21]=-1;
        imut[8][22]=-2;
        imut[8][23]=-1;
        imut[9][1]=-2;
        imut[9][2]=0;
        imut[9][3]=1;
        imut[9][4]=-1;
        imut[9][5]=-3;
        imut[9][6]=0;
        imut[9][7]=0;
        imut[9][8]=-2;
        imut[9][9]=8;
        imut[9][10]=-3;
        imut[9][11]=-3;
        imut[9][12]=-1;
        imut[9][13]=-2;
        imut[9][14]=-1;
        imut[9][15]=-2;
        imut[9][16]=-1;
        imut[9][17]=-2;
        imut[9][18]=-2;
        imut[9][19]=2;
        imut[9][20]=-3;
        imut[9][21]=0;
        imut[9][22]=0;
        imut[9][23]=-1;
        imut[10][1]=-1;
        imut[10][2]=-3;
        imut[10][3]=-3;
        imut[10][4]=-3;
        imut[10][5]=-1;
        imut[10][6]=-3;
        imut[10][7]=-3;
        imut[10][8]=-4;
        imut[10][9]=-3;
        imut[10][10]=4;
        imut[10][11]=2;
        imut[10][12]=-3;
        imut[10][13]=1;
        imut[10][14]=0;
        imut[10][15]=-3;
        imut[10][16]=-2;
        imut[10][17]=-1;
        imut[10][18]=-3;
        imut[10][19]=-1;
        imut[10][20]=3;
        imut[10][21]=-3;
        imut[10][22]=-3;
        imut[10][23]=-1;
        imut[11][1]=-1;
        imut[11][2]=-2;
        imut[11][3]=-3;
        imut[11][4]=-4;
        imut[11][5]=-1;
        imut[11][6]=-2;
        imut[11][7]=-3;
        imut[11][8]=-4;
        imut[11][9]=-3;
        imut[11][10]=2;
        imut[11][11]=4;
        imut[11][12]=-2;
        imut[11][13]=2;
        imut[11][14]=0;
        imut[11][15]=-3;
        imut[11][16]=-2;
        imut[11][17]=-1;
        imut[11][18]=-2;
        imut[11][19]=-1;
        imut[11][20]=1;
        imut[11][21]=-4;
        imut[11][22]=-3;
        imut[11][23]=-1;
        imut[12][1]=-1;
        imut[12][2]=2;
        imut[12][3]=0;
        imut[12][4]=-1;
        imut[12][5]=-3;
        imut[12][6]=1;
        imut[12][7]=1;
        imut[12][8]=-2;
        imut[12][9]=-1;
        imut[12][10]=-3;
        imut[12][11]=-2;
        imut[12][12]=5;
        imut[12][13]=-1;
        imut[12][14]=-3;
        imut[12][15]=-1;
        imut[12][16]=0;
        imut[12][17]=-1;
        imut[12][18]=-3;
        imut[12][19]=-2;
        imut[12][20]=-2;
        imut[12][21]=0;
        imut[12][22]=1;
        imut[12][23]=-1;
        imut[13][1]=-1;
        imut[13][2]=-1;
        imut[13][3]=-2;
        imut[13][4]=-3;
        imut[13][5]=-1;
        imut[13][6]=0;
        imut[13][7]=-2;
        imut[13][8]=-3;
        imut[13][9]=-2;
        imut[13][10]=1;
        imut[13][11]=2;
        imut[13][12]=-1;
        imut[13][13]=5;
        imut[13][14]=0;
        imut[13][15]=-2;
        imut[13][16]=-1;
        imut[13][17]=-1;
        imut[13][18]=-1;
        imut[13][19]=-1;
        imut[13][20]=1;
        imut[13][21]=-3;
        imut[13][22]=-1;
        imut[13][23]=-1;
        imut[14][1]=-2;
        imut[14][2]=-3;
        imut[14][3]=-3;
        imut[14][4]=-3;
        imut[14][5]=-2;
        imut[14][6]=-3;
        imut[14][7]=-3;
        imut[14][8]=-3;
        imut[14][9]=-1;
        imut[14][10]=0;
        imut[14][11]=0;
        imut[14][12]=-3;
        imut[14][13]=0;
        imut[14][14]=6;
        imut[14][15]=-4;
        imut[14][16]=-2;
        imut[14][17]=-2;
        imut[14][18]=1;
        imut[14][19]=3;
        imut[14][20]=-1;
        imut[14][21]=-3;
        imut[14][22]=-3;
        imut[14][23]=-1;
        imut[15][1]=-1;
        imut[15][2]=-2;
        imut[15][3]=-2;
        imut[15][4]=-1;
        imut[15][5]=-3;
        imut[15][6]=-1;
        imut[15][7]=-1;
        imut[15][8]=-2;
        imut[15][9]=-2;
        imut[15][10]=-3;
        imut[15][11]=-3;
        imut[15][12]=-1;
        imut[15][13]=-2;
        imut[15][14]=-4;
        imut[15][15]=7;
        imut[15][16]=-1;
        imut[15][17]=-1;
        imut[15][18]=-4;
        imut[15][19]=-3;
        imut[15][20]=-2;
        imut[15][21]=-2;
        imut[15][22]=-1;
        imut[15][23]=-2;
        imut[16][1]=1;
        imut[16][2]=-1;
        imut[16][3]=1;
        imut[16][4]=0;
        imut[16][5]=-1;
        imut[16][6]=0;
        imut[16][7]=0;
        imut[16][8]=0;
        imut[16][9]=-1;
        imut[16][10]=-2;
        imut[16][11]=-2;
        imut[16][12]=0;
        imut[16][13]=-1;
        imut[16][14]=-2;
        imut[16][15]=-1;
        imut[16][16]=4;
        imut[16][17]=1;
        imut[16][18]=-3;
        imut[16][19]=-2;
        imut[16][20]=-2;
        imut[16][21]=0;
        imut[16][22]=0;
        imut[16][23]=0;
        imut[17][1]=0;
        imut[17][2]=-1;
        imut[17][3]=0;
        imut[17][4]=-1;
        imut[17][5]=-1;
        imut[17][6]=-1;
        imut[17][7]=-1;
        imut[17][8]=-2;
        imut[17][9]=-2;
        imut[17][10]=-1;
        imut[17][11]=-1;
        imut[17][12]=-1;
        imut[17][13]=-1;
        imut[17][14]=-2;
        imut[17][15]=-1;
        imut[17][16]=1;
        imut[17][17]=5;
        imut[17][18]=-2;
        imut[17][19]=-2;
        imut[17][20]=0;
        imut[17][21]=-1;
        imut[17][22]=-1;
        imut[17][23]=0;
        imut[18][1]=-3;
        imut[18][2]=-3;
        imut[18][3]=-4;
        imut[18][4]=-4;
        imut[18][5]=-2;
        imut[18][6]=-2;
        imut[18][7]=-3;
        imut[18][8]=-2;
        imut[18][9]=-2;
        imut[18][10]=-3;
        imut[18][11]=-2;
        imut[18][12]=-3;
        imut[18][13]=-1;
        imut[18][14]=1;
        imut[18][15]=-4;
        imut[18][16]=-3;
        imut[18][17]=-2;
        imut[18][18]=11;
        imut[18][19]=2;
        imut[18][20]=-3;
        imut[18][21]=-4;
        imut[18][22]=-3;
        imut[18][23]=-2;
        imut[19][1]=-2;
        imut[19][2]=-2;
        imut[19][3]=-2;
        imut[19][4]=-3;
        imut[19][5]=-2;
        imut[19][6]=-1;
        imut[19][7]=-2;
        imut[19][8]=-3;
        imut[19][9]=2;
        imut[19][10]=-1;
        imut[19][11]=-1;
        imut[19][12]=-2;
        imut[19][13]=-1;
        imut[19][14]=3;
        imut[19][15]=-3;
        imut[19][16]=-2;
        imut[19][17]=-2;
        imut[19][18]=2;
        imut[19][19]=7;
        imut[19][20]=-1;
        imut[19][21]=-3;
        imut[19][22]=-2;
        imut[19][23]=-1;
        imut[20][1]=0;
        imut[20][2]=-3;
        imut[20][3]=-3;
        imut[20][4]=-3;
        imut[20][5]=-1;
        imut[20][6]=-2;
        imut[20][7]=-2;
        imut[20][8]=-3;
        imut[20][9]=-3;
        imut[20][10]=3;
        imut[20][11]=1;
        imut[20][12]=-2;
        imut[20][13]=1;
        imut[20][14]=-1;
        imut[20][15]=-2;
        imut[20][16]=-2;
        imut[20][17]=0;
        imut[20][18]=-3;
        imut[20][19]=-1;
        imut[20][20]=4;
        imut[20][21]=-3;
        imut[20][22]=-2;
        imut[20][23]=-1;
        imut[21][1]=-2;
        imut[21][2]=-1;
        imut[21][3]=3;
        imut[21][4]=4;
        imut[21][5]=-3;
        imut[21][6]=0;
        imut[21][7]=1;
        imut[21][8]=-1;
        imut[21][9]=0;
        imut[21][10]=-3;
        imut[21][11]=-4;
        imut[21][12]=0;
        imut[21][13]=-3;
        imut[21][14]=-3;
        imut[21][15]=-2;
        imut[21][16]=0;
        imut[21][17]=-1;
        imut[21][18]=-4;
        imut[21][19]=-3;
        imut[21][20]=-3;
        imut[21][21]=4;
        imut[21][22]=1;
        imut[21][23]=-1;
        imut[22][1]=-1;
        imut[22][2]=0;
        imut[22][3]=0;
        imut[22][4]=1;
        imut[22][5]=-3;
        imut[22][6]=3;
        imut[22][7]=4;
        imut[22][8]=-2;
        imut[22][9]=0;
        imut[22][10]=-3;
        imut[22][11]=-3;
        imut[22][12]=1;
        imut[22][13]=-1;
        imut[22][14]=-3;
        imut[22][15]=-1;
        imut[22][16]=0;
        imut[22][17]=-1;
        imut[22][18]=-3;
        imut[22][19]=-2;
        imut[22][20]=-2;
        imut[22][21]=1;
        imut[22][22]=4;
        imut[22][23]=-1;
        imut[23][1]=0;
        imut[23][2]=-1;
        imut[23][3]=-1;
        imut[23][4]=-1;
        imut[23][5]=-2;
        imut[23][6]=-1;
        imut[23][7]=-1;
        imut[23][8]=-1;
        imut[23][9]=-1;
        imut[23][10]=-1;
        imut[23][11]=-1;
        imut[23][12]=-1;
        imut[23][13]=-1;
        imut[23][14]=-1;
        imut[23][15]=-2;
        imut[23][16]=0;
        imut[23][17]=0;
        imut[23][18]=-2;
        imut[23][19]=-1;
        imut[23][20]=-1;
        imut[23][21]=-1;
        imut[23][22]=-1;
        imut[23][23]=-1;
    }

    public static void main(String[] args) throws Exception 
    {
        //try
        {
            SparkConf conf = new SparkConf()
                    .setAppName("PeptideNeedlemanWunschAlignment");
            JavaSparkContext jsc = new JavaSparkContext(conf);
            SQLContext sqlc = new SQLContext(jsc);

            Xml.Cfg pssmhcCfg = new Xml.Cfg(Xml.Utils.firstOrDef(args));
            XmlCfg appCfg = new XmlCfg(Xml.Utils.firstOrDef(args));
            Impl.PSSMHCpanSparkFunc   pssmhcSparkFunc = new Impl.PSSMHCpanSparkFunc(pssmhcCfg.pssmConfigs.get(0));
            Impl.ScoreFilterSparkFunc ic50FilterSparkFunc = new Impl.ScoreFilterSparkFunc(pssmhcCfg.ic50Threshold);
            
            Xml.PeptideGenCfg genCfg = pssmhcCfg.genCfg;
            JavaRDD<String> binders = 
                sqlc.range(genCfg.start, genCfg.end, 1, appCfg.partitions)
                .map(new Impl.PeptideGenSparkFunc(genCfg.peptideLength), Encoders.STRING())
                .toJavaRDD()
                .map(pssmhcSparkFunc)
                .filter(ic50FilterSparkFunc)
                .map(scp -> {return scp.peptide;});
            
            //binders.persist(StorageLevel.MEMORY_AND_DISK());
            //System.out.format("Binders with IC50<%d qnty %d of %d\n",  pssmhcCfg.ic50Threshold, binders.count(), (pssmhcCfg.end - pssmhcCfg.start));
            
            List<String> clusterCenters = 
                jsc.textFile(appCfg.peptidesFilename, appCfg.partitions)
                   .map(pssmhcSparkFunc)
                   .filter(ic50FilterSparkFunc)
                   .map(scp -> { return scp.peptide; })
                   .collect();
            jsc.broadcast(clusterCenters);
            
            System.out.format("Cluster centers with IC50<%d qnty %d\n", 
                pssmhcCfg.ic50Threshold, clusterCenters.size());

            AssignBindersToClusters.PepSimSparkFunc simFunc = new AssignBindersToClusters.PepSimSparkFunc(1/appCfg.minSimilarity);
            simFunc.SetMatrix(SubstMatrices.get(appCfg.matrix));
            
            // {Bn -> Ck}
            final double minSimilarity = appCfg.minSimilarity;
            JavaPairRDD<String, Impl.ScoredPeptide> simPairs = binders.flatMapToPair(binder -> {
                        ArrayList<Tuple2<String, String>> res = new ArrayList<>();
                        for (String center : clusterCenters) {
                            res.add(new Tuple2<>(binder, center));
                        }
                        return res.iterator();
                    }).mapToPair(tuple -> {
                        final double score = calcNWscore(tuple._1, tuple._2).score;
                        return new Tuple2<>(tuple._1, new Impl.ScoredPeptide(tuple._2, score));
                    }).filter(tuple -> {
                        return (tuple._2.score > minSimilarity); 
                    }).persist(StorageLevel.MEMORY_AND_DISK());
            
            System.out.format("Pairs with similarity>%.2f qnty %d\n", appCfg.minSimilarity, simPairs.count());
           
            // {Bn -> (Ck_max, Snk_max)}
            JavaPairRDD<String, Impl.ScoredPeptide> simMaxPairs = simPairs.reduceByKey(
                (scp1, scp2) -> { return (scp1.score > scp2.score) ? scp1 : scp2; } 
            );

            // {Ck -> (Bn, Snk)}
            JavaPairRDD<String, Impl.ScoredPeptide> simMaxPairsInv = simMaxPairs.mapToPair(tuple -> {
                Impl.ScoredPeptide binderWithSim = new Impl.ScoredPeptide(tuple._1, tuple._2.score);
                return new Tuple2<>(tuple._2.peptide, binderWithSim);
            });
            simMaxPairsInv = simMaxPairsInv.coalesce(appCfg.partitions);
            simMaxPairsInv.persist(StorageLevel.MEMORY_AND_DISK());
                        
            Map<String, Long> clustersCount = simMaxPairsInv.countByKey();
            System.out.format("Clusters qnty %d\n", clustersCount.size());
            FormatMap(clustersCount);

            JavaPairRDD<String, Iterable<Impl.ScoredPeptide>> simMaxPairsGrp = 
                simMaxPairsInv.groupByKey()
                              .coalesce(1);
            simMaxPairsGrp.saveAsTextFile("output-clusters");
            
            jsc.close();
        }
        /*
        catch (Exception ex)
        {
            System.err.print(ex.toString() + "\n");
        }
        */
    }

}
