package org.PeptideClustering;

import java.nio.file.FileSystems;
import org.PSSMHC.Xml;
import org.w3c.dom.Element;

public class XmlCfg
{
    int partitions = 0, clustersQnty = 0, maxTrials = 0;
    String peptidesFilename = "";
    String matrix = "";
    double minSimilarity = Double.NaN;
    int centersIc50Threshold = -1;

    public XmlCfg(String xmlFilename) throws Exception
    {
        Element root = Xml.Utils.parseXml(xmlFilename);
        Element elem = Xml.Utils.getChildNode(root, "clustering");
        if (elem == null) { return; }

        String pathPrefix       = Xml.Utils.getChildAttr(root, "system", "pathPrefix");
        String pepRelFilename = elem.getAttribute("peptidesFilePath");
        
        peptidesFilename =      FileSystems.getDefault().getPath(pathPrefix, pepRelFilename).toString();
        centersIc50Threshold = Integer.parseInt(elem.getAttribute("centersIc50Threshold"));
        matrix          =                       elem.getAttribute("matrix");
        clustersQnty    = Integer.parseInt(     elem.getAttribute("clustersQnty"));
        maxTrials       = Integer.parseInt(     elem.getAttribute("maxTrials"));
        minSimilarity   = Double.parseDouble(   elem.getAttribute("minSimilarity"));
        partitions      = Integer.parseInt(Xml.Utils.getChildAttr(root, "spark", "partitions"));
    }
}
