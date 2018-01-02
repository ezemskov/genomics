package org.PSSMHC;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.nio.file.FileSystems;
import java.nio.file.FileSystem;
import java.util.HashMap;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class Xml
{
    public static class StringIntPair
    {
        StringIntPair(String first_, int second_)
        {
            first = first_;
            second = second_;
        }
            
        @Override
        public int hashCode()
        {
            return String.format("%s_%d", first, second).hashCode();
        }

        @Override
        public boolean equals(Object obj)
        {
            StringIntPair ap = (StringIntPair)obj;
            return (ap != null) && (first.equals(ap.first)) && (second == ap.second);
        }

        public String first; 
        public int second; 
    }
    
    public static class Utils
    {
        public static String DefaultXmlFilename = "PeptideCfg.xml";

        public static String firstOrDef(String[] args)
        {
            return (args.length > 0) ? args[0] : DefaultXmlFilename;
        }

        public static Element parseXml(String filename) throws ParserConfigurationException, SAXException, IOException
        {
            String filenameOrDefault = filename.isEmpty() ? DefaultXmlFilename : filename;
            DocumentBuilder builder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            Document domDoc = builder.parse(new File(filenameOrDefault));
            assert (domDoc != null); //should throw instead
            return domDoc.getDocumentElement();
        }

        public static List<Element> getAllChildNodes(Node parent, String childName)
        {
            ArrayList<Element> res = new ArrayList<>();
            if (parent == null)
            {
                return res;
            }
            
            NodeList children = parent.getChildNodes();
            for (int i = 0; i < children.getLength(); ++i)
            {
                Node child = children.item(i);
                if (child.getNodeName().equals(childName) &&
                    (child.getNodeType() == Node.ELEMENT_NODE))
                {
                    res.add((Element)child);
                }
            }
            
            return res;
        }
        
        public static Element getChildNode(Node parent, String childName) throws Exception
        {
            List<Element> res = getAllChildNodes(parent, childName);
            if (res.isEmpty())
            {
                throw new Exception("<" + parent.getNodeName() + "><" + childName + "> is not found");
            }
            
            return res.get(0);
        }

        public static String getChildAttr(Node parent, String childName, String childAttr) throws Exception
        {
            final String res = getChildNode(parent, childName).getAttribute(childAttr);
            if (res.isEmpty())
            {
                throw new Exception("<" + parent.getNodeName() + "><" + 
                    childName + " " + childAttr + "=\"...\"/> is not found or empty");
            }
            return res;
        }
    }

    public static class PSSMCfg
    {
        public String pathPrefix = "";
        public String pssmListRelPath = "";
        public String allele = "";
        public int    peptideLength = -1;
    }

    public static class PeptideGenCfg
    {
        public long start = -1;
        public long end = -1;
        public int  peptideLength = -1;
    }
    
    public static class Cfg
    {
        public PeptideGenCfg genCfg = new PeptideGenCfg();
        public int partitions;
        public boolean doBinderPersist, doBinderStore, doBinderCount;
        public int ic50Threshold;
        
        public ArrayList<StringIntPair> peptideFiles = new ArrayList<>();
        public ArrayList<PSSMCfg> pssmConfigs = new ArrayList<>();
        public String pssmListPath = "";
        
        public Cfg(String xmlFilename) throws Exception
        {
            final Element root = Utils.parseXml(xmlFilename);
            
            genCfg.start     = ParseLongWithSuffix(Utils.getChildAttr(root, "generator", "start"));
            genCfg.end       = genCfg.start + ParseLongWithSuffix(Utils.getChildAttr(root, "generator", "qnty"));
            genCfg.peptideLength = Integer.parseInt(Utils.getChildAttr(root, "generator", "peptideLength"));
            doBinderPersist  = Utils.getChildAttr(root, "binders", "doBinderPersist").equals("1");
            doBinderStore    = Utils.getChildAttr(root, "binders", "doBinderStore").equals("1");
            doBinderCount    = Utils.getChildAttr(root, "binders", "doBinderCount").equals("1");
            ic50Threshold    = Integer.parseInt(Utils.getChildAttr(root, "binders", "ic50Threshold"));            
            partitions       = Integer.parseInt(Utils.getChildAttr(root, "spark", "partitions"));
                        
            final FileSystem fs     = FileSystems.getDefault();
            final String pathPrefix = Utils.getChildAttr(root, "system", "pathPrefix");
            
            {
                final Element elem = Utils.getChildNode(root, "PSSMHC");
                List<Element> elemChildren = Utils.getAllChildNodes(elem, "allele");
                for(Element child : elemChildren)
                {
                    PSSMCfg pssmCfg = new PSSMCfg();
                    pssmCfg.pathPrefix = pathPrefix;
                    pssmCfg.pssmListRelPath = Utils.getChildAttr(root, "PSSMHC", "fileList");
                    pssmCfg.allele = child.getTextContent();
                    pssmConfigs.add(pssmCfg);
                }
            }{
                final Element elem = Utils.getChildNode(root, "binders");
                List<Element> elemChildren = Utils.getAllChildNodes(elem, "peptidesFile");
                for(Element child : elemChildren)
                {
                    final String pepRelPath = child.getAttribute("path");
                    peptideFiles.add(new StringIntPair(
                        fs.getPath(pathPrefix, pepRelPath).toString(), 
                        Integer.parseInt(child.getAttribute("peptideLength"))   //todo : proper error message is attribute not found
                    ));
                }
            }
        }

        public PSSMCfg getSinglePSSMCfg() throws Exception
        {
            if (pssmConfigs.size() != 1)
            {
                throw new Exception("<PSSMHC><allele> is not found or not unique");
            }
            
            return pssmConfigs.get(0);
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