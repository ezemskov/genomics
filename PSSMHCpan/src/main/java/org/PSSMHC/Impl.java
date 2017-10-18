package org.PSSMHC;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.sql.Row;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;
import org.w3c.dom.Document; 
import org.w3c.dom.Element; 
import org.w3c.dom.NodeList; 
import org.w3c.dom.Node; 

public class Impl
{
    public static class PSSMHCpanSparkFunc 
                            extends PSSMHCpan
                            implements Serializable,
                            Function<String,ScoredPeptide>
    {
        public PSSMHCpanSparkFunc(String xmlFilename) throws Exception
        {
            super(xmlFilename);
        }
            
        @Override
        public ScoredPeptide call(String peptide)
        {
            return new ScoredPeptide(peptide, ScoreOnePeptide(peptide));
        }
    }

    public static class PeptideBloomSparkFunc 
                            extends PeptideBloomFilter
                            implements Serializable,
                            VoidFunction<ScoredPeptide>
    {
        @Override
        public void call(ScoredPeptide scPep)
        {
            Add(scPep.peptide);
        }
    }

    public static class PeptideGenSparkFunc 
                            extends PeptideGen
                            implements Serializable, 
                            MapFunction<Row, String>
    {
        @Override
        public String call(Row idx)
        {
            return Generate(idx.getLong(0));
        };
    }

    public static class ScoreFilterSparkFunc 
                            implements Serializable,
                            Function<ScoredPeptide,Boolean>
    {
        public ScoreFilterSparkFunc(double scoreThreshold_)
        {
            scoreThreshold = scoreThreshold_;
        }

        @Override
        public Boolean call(ScoredPeptide scPep)
        {
            return (scPep.score < scoreThreshold);
        }

        double scoreThreshold;
    }

    public static class ScoredPeptide implements Serializable
    {
        public ScoredPeptide(String peptide, double score) 
        { 
            this.peptide = peptide; 
            this.score = score;
        }

        @Override
        public String toString()
        {
            return String.format("%s,%.2f", peptide, score);
        }

        public String peptide; 
        public double score; 
    }

    public static class XmlUtils
    {
        public static String DefaultXmlFilename = "PeptideCfg.xml";
        
        public static String firstOrDef(String args[])
        {
            return (args.length > 0) ? args[0] : DefaultXmlFilename;
        }
        
        public static Element parseXml(String filename) 
            throws ParserConfigurationException, SAXException, IOException
        {
            String filenameOrDefault = filename.isEmpty() ? DefaultXmlFilename : filename;
            DocumentBuilder builder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            Document domDoc = builder.parse(new File(filenameOrDefault));
            assert(domDoc != null);    //should throw instead
            return domDoc.getDocumentElement();
        }

        public static Element getChildElem(Node parent, String childName)
        {
            NodeList children = parent.getChildNodes();
            for (int i=0; i<children.getLength(); ++i)
            {
                Node child = children.item(i);
                if (!child.getNodeName().equals(childName))
                {
                    continue;
                }
                
                if (child.getNodeType() != Node.ELEMENT_NODE)
                {
                    return null;
                }
                    
                return (Element)child;
            }
            
            return null;
        }

        public static String getChildAttr(Node parent, String childName, String childAttr)
        {
            Element child = getChildElem(parent, childName);
            return (child != null) ? child.getAttribute(childAttr) : "";
        }
    }

    public static class XmlCfg
    {
        public long start;
        public long end;
        public int partitions;
        public boolean doScore, doBinderPersist, doBinderStore, doBinderCount;
        public int ic50Threshold;
        
        public XmlCfg(String xmlFilename) throws Exception
        {
            Element root = XmlUtils.parseXml(xmlFilename);
            Element elem = XmlUtils.getChildElem(root, "generator");
            if (elem == null) { return; }
            
            start            = ParseLongWithSuffix(elem.getAttribute("start"));
            end      = start + ParseLongWithSuffix(elem.getAttribute("qnty"));
            doScore          = elem.getAttribute("doScore").equals("1");
            doBinderPersist  = elem.getAttribute("doBinderPersist").equals("1");
            doBinderStore    = elem.getAttribute("doBinderStore").equals("1");
            doBinderCount    = elem.getAttribute("doBinderCount").equals("1");
            ic50Threshold    = Integer.parseInt(elem.getAttribute("ic50Threshold"));            
            partitions       = Integer.parseInt(XmlUtils.getChildAttr(root, "spark", "partitions"));
        }

        private static long ParseLongWithSuffix(String val)
        {
            int suffixPos = val.length()-1;
            Character suffix = val.toUpperCase().charAt(suffixPos);
            Long multiplier = suffixes.get(suffix);
            if ((val.length() < 2) || (multiplier == null))
            {
                return Long.parseUnsignedLong(val);
            }

            return multiplier * Long.parseUnsignedLong(val.substring(0, suffixPos));
        }

        private static final HashMap<Character, Long> suffixes;
        static
        {
            suffixes = new HashMap<>();
            suffixes.put('K', (long)1E3);
            suffixes.put('M', (long)1E6);
            suffixes.put('G', (long)1E9);
            suffixes.put('T', (long)1E12);
        }
    }
}
