package org.PeptideClustering;

import org.PSSMHC.Xml;
import org.w3c.dom.Element;

public class XmlCfg
{
    int partitions = 0, clustersQnty = 0, maxTrials = 0;
    String peptidesFilename = "";
    String matrix = "";
    double minSimilarity = Double.NaN;

    public XmlCfg(String xmlFilename) throws Exception
    {
        Element root = Xml.Utils.parseXml(xmlFilename);
        Element elem = Xml.Utils.getChildElem(root, "clustering");
        if (elem == null) { return; }

        peptidesFilename =                    elem.getAttribute("peptidesFilePath");
        matrix           =                    elem.getAttribute("matrix");
        clustersQnty     = Integer.parseInt(  elem.getAttribute("clustersQnty"));
        maxTrials        = Integer.parseInt(  elem.getAttribute("maxTrials"));
        minSimilarity    = Double.parseDouble(elem.getAttribute("minSimilarity"));
        partitions       = Integer.parseInt(Xml.Utils.getChildAttr(root, "spark", "partitions"));
    }
}
