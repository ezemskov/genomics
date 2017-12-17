package org.PSSMHC;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystems;
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

        public static Element getChildElem(Node parent, String childName)
        {
            NodeList children = parent.getChildNodes();
            for (int i = 0; i < children.getLength(); ++i)
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
                return (Element) child;
            }
            return null;
        }

        public static String getChildAttr(Node parent, String childName, String childAttr)
        {
            Element child = getChildElem(parent, childName);
            return (child != null) ? child.getAttribute(childAttr) : "";
        }
    }

    public static class Cfg
    {
        public long start;
        public long end;
        public int partitions;
        public boolean doScore, doBinderPersist, doBinderStore, doBinderCount;
        public int ic50Threshold;
        String peptidesFilename = "";
        int peptideLength = -1;
        String alleleName = "";
        
        public Cfg(String xmlFilename) throws Exception
        {
            Element root = Utils.parseXml(xmlFilename);
            Element elem = Utils.getChildElem(root, "generator");
            if (elem == null) { return; }
            
            start            = ParseLongWithSuffix(elem.getAttribute("start"));
            end      = start + ParseLongWithSuffix(elem.getAttribute("qnty"));
            doScore          = elem.getAttribute("doScore").equals("1");
            doBinderPersist  = elem.getAttribute("doBinderPersist").equals("1");
            doBinderStore    = elem.getAttribute("doBinderStore").equals("1");
            doBinderCount    = elem.getAttribute("doBinderCount").equals("1");
            ic50Threshold    = Integer.parseInt(elem.getAttribute("ic50Threshold"));            
            partitions       = Integer.parseInt(Utils.getChildAttr(root, "spark", "partitions"));
            peptideLength    = Integer.parseInt(Utils.getChildAttr(root, "PSSMHC", "peptideLength"));
            alleleName       = Utils.getChildAttr(root, "PSSMHC", "alleleName");

            String pathPrefix     = Xml.Utils.getChildAttr(root, "system", "pathPrefix");
            String pepRelFilename = elem.getAttribute("peptidesFilePath");
            
            if (!pepRelFilename.isEmpty())
            {
                peptidesFilename = FileSystems.getDefault().getPath(pathPrefix, pepRelFilename).toString();
            }
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