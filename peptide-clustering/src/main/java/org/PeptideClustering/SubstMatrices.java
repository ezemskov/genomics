package org.PeptideClustering;

import java.io.Serializable;
import java.util.HashMap;
import org.PSSMHC.Impl;

class SubstMatrixRow extends HashMap<Character, Double>         implements Serializable
{
    SubstMatrixRow() { super(Impl.Consts.aLen, 1.0f); }
}

class SubstMatrix    extends HashMap<Character, SubstMatrixRow> implements Serializable 
{
    public SubstMatrix(double[][] vals)
    {
        super(Impl.Consts.aLen, 1.0f);
        
        final String ExcMsg =  "Wrong substitution matrix size";
        if (vals.length != Impl.Consts.aLen) 
        {
            throw new RuntimeException(ExcMsg);
        }

        for (int i=0; i<vals.length; ++i) 
        {
            if (vals[i].length != Impl.Consts.aLen) 
            {
                throw new RuntimeException(ExcMsg);
            }
            
            SubstMatrixRow row = new SubstMatrixRow();
            for (int j=0; j<vals.length; ++j) 
            {
                row.put(Impl.Consts.alphabet.charAt(j), vals[i][j]);
            }
            
            put(Impl.Consts.alphabet.charAt(i), row);
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
    }    

    static class Blosum50 extends SubstMatrix implements Serializable
    {
        public Blosum50() { super(vals); }
        
        static double[][] vals = new double[][] {
          {  3.4217, -0.3440, -1.2002, -0.6534, -1.7445,  0.2784, -1.2822, -0.9673, -0.8392, -1.2123, -0.4572, -0.9498, -0.9689, -0.5729, -1.1626,  0.8948, -0.0958, -0.0535, -2.2183, -1.4913},
          { -0.3440,  8.3766, -2.4505, -2.2656, -1.6474, -1.8628, -1.9012, -1.5362, -2.2414, -1.2426, -1.1088, -1.4691, -2.6218, -2.1100, -2.4498, -0.5381, -0.5800, -0.6128, -3.3595, -1.7849},
          { -1.2002, -2.4505,  5.2230,  1.4586, -3.1360, -0.8514, -0.3564, -2.8677, -0.3351, -2.8461, -2.4018,  1.2792, -0.9897, -0.0968, -1.1190, -0.2359, -0.5544, -2.4237, -3.4069, -1.9743},
          { -0.6534, -2.2656,  1.4586,  4.2994, -2.2654, -1.7784, -0.2674, -2.5338,  0.8252, -2.1273, -1.4495, -0.2413, -0.5072,  1.4820, -0.0690, -0.3620, -0.5760, -2.1233, -1.7353, -1.2355},
          { -1.7445, -1.6474, -3.1360, -2.2654,  5.4607, -2.5449, -0.8120, -0.0298, -2.3364,  0.6708,  0.1451, -2.3774, -2.8244, -2.5240, -2.2294, -1.7094, -1.6202, -0.4364,  0.8618,  2.5553},
          {  0.2784, -1.8628, -0.8514, -1.7784, -2.5449,  5.0682, -1.4678, -2.8659, -1.2214, -2.8155, -1.8900, -0.0908, -1.3785, -1.2813, -1.9036, -0.0022, -1.3046, -2.6243, -1.9719, -2.1953},
          { -1.2822, -1.9012, -0.3564, -0.2674, -0.8120, -1.4678,  6.7706, -2.5618, -0.1616, -1.7384, -0.7914,  0.9586, -1.5623,  0.5479, -0.0536, -0.5637, -1.2313, -2.5472, -2.1480,  1.3015},
          { -0.9673, -1.5362, -2.8677, -2.5338, -0.0298, -2.8659, -2.5618,  3.5403, -2.1890,  1.5256,  1.0483, -2.2985, -1.9360, -1.9847, -2.3997, -1.7479, -0.4328,  2.4192, -1.8730, -0.5512},
          { -0.8392, -2.2414, -0.3351,  0.8252, -2.3364, -1.2214, -0.1616, -2.1890,  3.9128, -2.1226, -1.0907,  0.1693, -0.7779,  1.0088,  2.0918, -0.3266, -0.4866, -1.8713, -2.2306, -1.1652},
          { -1.2123, -1.2426, -2.8461, -2.1273,  0.6708, -2.8155, -1.7384,  1.5256, -2.1226,  3.4694,  1.6798, -2.5540, -2.3449, -1.6579, -1.7039, -1.9394, -0.7950,  0.8104, -1.4666, -0.4406},
          { -0.4572, -1.1088, -2.4018, -1.4495,  0.1451, -1.8900, -0.7914,  1.0483, -1.0907,  1.6798,  4.5354, -1.4191, -1.7523, -0.1180, -1.1212, -1.1031, -0.3613,  0.5579, -0.7861, -0.2474},
          { -0.9498, -1.4691,  1.2792, -0.2413, -2.3774, -0.0908,  0.9586, -2.2985,  0.1693, -2.5540, -1.4191,  4.8040, -1.4535,  0.0991, -0.4561,  0.5840,  0.2981, -2.2203, -2.8247, -1.3957},
          { -0.9689, -2.6218, -0.9897, -0.5072, -2.8244, -1.3785, -1.5623, -1.9360, -0.7779, -2.3449, -1.7523, -1.4535,  6.7020, -0.8282, -1.8944, -0.7940, -0.8356, -1.8623, -2.5025, -2.1898},
          { -0.5729, -2.1100, -0.0968,  1.4820, -2.5240, -1.2813,  0.5479, -1.9847,  1.0088, -1.6579, -0.1180,  0.0991, -0.8282,  4.4636,  0.8816,  0.1918, -0.6069, -1.6865, -0.9692, -0.8607},
          { -1.1626, -2.4498, -1.1190, -0.0690, -2.2294, -1.9036, -0.0536, -2.3997,  2.0918, -1.7039, -1.1212, -0.4561, -1.8944,  0.8816,  4.8541, -0.6277, -0.8824, -1.9288, -1.8383, -0.9194},
          {  0.8948, -0.5381, -0.2359, -0.3620, -1.7094, -0.0022, -0.5637, -1.7479, -0.3266, -1.9394, -1.1031,  0.5840, -0.7940,  0.1918, -0.6277,  3.3043,  1.1647, -1.1185, -2.6998, -1.2385},
          { -0.0958, -0.5800, -0.5544, -0.5760, -1.6202, -1.3046, -1.2313, -0.4328, -0.4866, -0.7950, -0.3613,  0.2981, -0.8356, -0.6069, -0.8824,  1.1647,  3.6581,  0.1679, -1.9846, -1.0676},
          { -0.0535, -0.6128, -2.4237, -2.1233, -0.4364, -2.6243, -2.5472,  2.4192, -1.8713,  0.8104,  0.5579, -2.2203, -1.8623, -1.6865, -1.9288, -1.1185,  0.1679,  3.2807, -2.0889, -0.9186},
          { -2.2183, -3.3595, -3.4069, -1.7353,  0.8618, -1.9719, -2.1480, -1.8730, -2.2306, -1.4666, -0.7861, -2.8247, -2.5025, -0.9692, -1.8383, -2.6998, -1.9846, -2.0889,  9.9418,  1.6396},
          { -1.4913, -1.7849, -1.9743, -1.2355,  2.5553, -2.1953,  1.3015, -0.5512, -1.1652, -0.4406, -0.2474, -1.3957, -2.1898, -0.8607, -0.9194, -1.2385, -1.0676, -0.9186,  1.6396,  5.5702}
        };
    }
            
    static class Blosum62 extends SubstMatrix implements Serializable
    {
        public Blosum62() { super(vals); }

        static double[][] vals = new double[][] {
          {  3.9291, -0.4085, -1.7534, -0.8639, -2.2101,  0.1596, -1.6251, -1.3218, -0.7340, -1.4646, -0.9353, -1.5307, -0.8143, -0.8040, -1.4135,  1.1158, -0.0454, -0.1894, -2.5269, -1.7640},
          { -0.4085,  8.5821, -3.4600, -3.6125, -2.3755, -2.5004, -2.9878, -1.2277, -3.0363, -1.2775, -1.4198, -2.6598, -2.7952, -2.9019, -3.3892, -0.8750, -0.8667, -0.8077, -2.3041, -2.4071},
          { -1.7534, -3.4600,  5.7742,  1.5103, -3.4839, -1.3135, -1.1189, -3.1212, -0.7018, -3.6057, -3.0585,  1.2717, -1.4801, -0.3134, -1.6058, -0.2610, -1.0507, -3.1426, -4.2143, -3.0650},
          { -0.8639, -3.6125,  1.5103,  4.9028, -3.1924, -2.1102, -0.1177, -3.1944,  0.7753, -2.8465, -1.9980, -0.2680, -1.1162,  1.8546, -0.1154, -0.1469, -0.8633, -2.4423, -2.8354, -2.0205},
          { -2.2101, -2.3755, -3.4839, -3.1924,  6.0461, -3.1074, -1.2342, -0.1609, -3.0787,  0.4148,  0.0126, -2.9940, -3.5973, -3.1644, -2.7863, -2.3690, -2.1076, -0.8490,  0.9176,  2.9391},
          {  0.1596, -2.5004, -1.3135, -2.1102, -3.1074,  5.5633, -2.0409, -3.7249, -1.5280, -3.6270, -2.6766, -0.4228, -2.1335, -1.7852, -2.3041, -0.2925, -1.5754, -3.1387, -2.4915, -3.0398},
          { -1.6251, -2.9878, -1.1189, -0.1177, -1.2342, -2.0409,  7.5111, -3.2316, -0.7210, -2.7867, -1.5513,  0.5785, -2.1609,  0.4480, -0.2499, -0.8816, -1.6859, -3.1175, -2.3422,  1.6926},
          { -1.3218, -1.2277, -3.1212, -3.1944, -0.1609, -3.7249, -3.2316,  3.9985, -2.6701,  1.5216,  1.1268, -3.2170, -2.7567, -2.7696, -2.9902, -2.3482, -0.7176,  2.5470, -2.5805, -1.3314},
          { -0.7340, -3.0363, -0.7018,  0.7753, -3.0787, -1.5280, -0.7210, -2.6701,  4.5046, -2.4468, -1.3547, -0.1790, -1.0136,  1.2726,  2.1087, -0.2034, -0.6696, -2.2624, -2.9564, -1.8200},
          { -1.4646, -1.2775, -3.6057, -2.8465,  0.4148, -3.6270, -2.7867,  1.5216, -2.4468,  3.8494,  1.9918, -3.3789, -2.8601, -2.1339, -2.1546, -2.4426, -1.1975,  0.7884, -1.6319, -1.0621},
          { -0.9353, -1.4198, -3.0585, -1.9980,  0.0126, -2.6766, -1.5513,  1.1268, -1.3547,  1.9918,  5.3926, -2.1509, -2.4764, -0.4210, -1.3671, -1.4809, -0.6663,  0.6872, -1.4248, -0.9949},
          { -1.5307, -2.6598,  1.2717, -0.2680, -2.9940, -0.4228,  0.5785, -3.2170, -0.1790, -3.3789, -2.1509,  5.6532, -2.0004,  0.0017, -0.4398,  0.6009, -0.0461, -2.8763, -3.6959, -2.0818},
          { -0.8143, -2.7952, -1.4801, -1.1162, -3.5973, -2.1335, -2.1609, -2.7567, -1.0136, -2.8601, -2.4764, -2.0004,  7.3646, -1.2819, -2.1086, -0.8090, -1.0753, -2.3487, -3.6542, -2.9198},
          { -0.8040, -2.9019, -0.3134,  1.8546, -3.1644, -1.7852,  0.4480, -2.7696,  1.2726, -2.1339, -0.4210,  0.0017, -1.2819,  5.2851,  0.9828, -0.1011, -0.6753, -2.1984, -1.9465, -1.4211},
          { -1.4135, -3.3892, -1.6058, -0.1154, -2.7863, -2.3041, -0.2499, -2.9902,  2.1087, -2.1546, -1.3671, -0.4398, -2.1086,  0.9828,  5.4735, -0.7648, -1.1223, -2.5026, -2.6794, -1.6939},
          {  1.1158, -0.8750, -0.2610, -0.1469, -2.3690, -0.2925, -0.8816, -2.3482, -0.2034, -2.4426, -1.4809,  0.6009, -0.8090, -0.1011, -0.7648,  3.8844,  1.3811, -1.6462, -2.7519, -1.6858},
          { -0.0454, -0.8667, -1.0507, -0.8633, -2.1076, -1.5754, -1.6859, -0.7176, -0.6696, -1.1975, -0.6663, -0.0461, -1.0753, -0.6753, -1.1223,  1.3811,  4.5453, -0.0555, -2.4289, -1.6060},
          { -0.1894, -0.8077, -3.1426, -2.4423, -0.8490, -3.1387, -3.1175,  2.5470, -2.2624,  0.7884,  0.6872, -2.8763, -2.3487, -2.1984, -2.5026, -1.6462, -0.0555,  3.7689, -2.8343, -1.2075},
          { -2.5269, -2.3041, -4.2143, -2.8354,  0.9176, -2.4915, -2.3422, -2.5805, -2.9564, -1.6319, -1.4248, -3.6959, -3.6542, -1.9465, -2.6794, -2.7519, -2.4289, -2.8343, 10.5040,  2.1542},
          { -1.7640, -2.4071, -3.0650, -2.0205,  2.9391, -3.0398,  1.6926, -1.3314, -1.8200, -1.0621, -0.9949, -2.0818, -2.9198, -1.4211, -1.6939, -1.6858, -1.6060, -1.2075,  2.1542,  6.5950}
        };
    }

    static class Blosum80 extends SubstMatrix implements Serializable
    {
        public Blosum80() { super(vals); }

        static double[][] vals = new double[][] {
          {  4.5099, -0.9010, -2.2982, -1.0162, -2.6670, -0.1265, -1.9179, -1.7603, -0.9363, -1.9737, -1.3538, -1.9422, -0.7500, -1.0467, -1.6972,  1.2358, -0.0594, -0.4155, -3.3860, -2.3937},
          { -0.9010,  8.7434, -4.4500, -4.9539, -2.6819, -3.7561, -4.3516, -1.5659, -4.1040, -2.0427, -2.0034, -3.4702, -3.7868, -3.5213, -4.2097, -1.5911, -1.4656, -1.3128, -3.4526, -3.3986},
          { -2.2982, -4.4500,  6.3737,  1.4195, -4.1855, -1.7744, -1.5011, -4.4480, -1.1261, -4.6850, -3.9776,  1.3280, -2.2886, -0.7794, -2.1381, -0.7390, -1.4221, -4.0619, -5.5762, -4.0560},
          { -1.0162, -4.9539,  1.4195,  5.6125, -4.0119, -2.6485, -0.3000, -3.8400,  0.5139, -3.7195, -2.4448, -0.6040, -1.5647,  1.8618, -0.5312, -0.4858, -1.0894, -2.8736, -4.1081, -3.1739},
          { -2.6670, -2.6819, -4.1855, -4.0119,  6.4914, -4.0139, -1.6126, -0.5011, -3.6451,  0.3123, -0.3255, -3.7494, -4.1521, -3.6203, -3.5991, -2.8758, -2.3350, -1.2479,  0.2471,  2.9504},
          { -0.1265, -3.7561, -1.7744, -2.6485, -4.0139,  5.9573, -2.7381, -4.8773, -2.0972, -4.4956, -3.6111, -0.7883, -3.0515, -2.4717, -2.8113, -0.7010, -2.0451, -3.9890, -3.8385, -4.2467},
          { -1.9179, -4.3516, -1.5011, -0.3000, -1.6126, -2.7381,  8.0125, -3.9099, -0.8684, -3.3396, -2.4213,  0.3379, -2.5054,  0.7913, -0.2233, -1.1928, -1.7802, -3.5824, -2.7152,  1.7268},
          { -1.7603, -1.5659, -4.4480, -3.8400, -0.5011, -4.8773, -3.9099,  4.5664, -3.3483,  1.4710,  1.1938, -3.9114, -3.6139, -3.3880, -3.4802, -2.7995, -1.0270,  2.6391, -3.0862, -1.7816},
          { -0.9363, -4.1040, -1.1261,  0.5139, -3.6451, -2.0972, -0.8684, -3.3483,  5.3224, -2.9732, -1.8080, -0.1835, -1.4899,  1.2158,  2.2647, -0.5719, -0.8855, -2.8673, -4.1007, -2.5832},
          { -1.9737, -2.0427, -4.6850, -3.7195,  0.3123, -4.4956, -3.3396,  1.4710, -2.9732,  4.3161,  2.1719, -4.0036, -3.4443, -2.5945, -2.9253, -2.8807, -1.6687,  0.5749, -2.3768, -1.5692},
          { -1.3538, -2.0034, -3.9776, -2.4448, -0.3255, -3.6111, -2.4213,  1.1938, -1.8080,  2.1719,  6.3022, -2.7797, -2.9325, -0.3474, -1.9663, -2.0090, -0.7996,  0.5835, -1.6696, -1.7253},
          { -1.9422, -3.4702,  1.3280, -0.6040, -3.7494, -0.7883,  0.3379, -3.9114, -0.1835, -4.0036, -2.7797,  6.3279, -2.6593, -0.1233, -0.7428,  0.4415, -0.2784, -3.5027, -4.3518, -2.7578},
          { -0.7500, -3.7868, -2.2886, -1.5647, -4.1521, -3.0515, -2.5054, -3.6139, -1.4899, -3.4443, -2.9325, -2.6593,  7.8434, -1.7887, -2.3324, -1.2354, -1.6712, -2.8670, -4.9796, -3.9142},
          { -1.0467, -3.5213, -0.7794,  1.8618, -3.6203, -2.4717,  0.7913, -3.3880,  1.2158, -2.5945, -0.3474, -0.1233, -1.7887,  6.1201,  0.9590, -0.4386, -0.9304, -2.5660, -2.5874, -2.2290},
          { -1.6972, -4.2097, -2.1381, -0.5312, -3.5991, -2.8113, -0.2233, -3.4802,  2.2647, -2.9253, -1.9663, -0.7428, -2.3324,  0.9590,  6.0869, -1.0519, -1.4817, -2.9986, -3.5298, -2.5184},
          {  1.2358, -1.5911, -0.7390, -0.4858, -2.8758, -0.7010, -1.1928, -2.7995, -0.5719, -2.8807, -2.0090,  0.4415, -1.2354, -0.4386, -1.0519,  4.7043,  1.4669, -2.0367, -3.7697, -2.2280},
          { -0.0594, -1.4656, -1.4221, -1.0894, -2.3350, -2.0451, -1.7802, -1.0270, -0.8855, -1.6687, -0.7996, -0.2784, -1.6712, -0.9304, -1.4817,  1.4669,  5.2671, -0.3314, -3.6211, -2.1513},
          { -0.4155, -1.3128, -4.0619, -2.8736, -1.2479, -3.9890, -3.5824,  2.6391, -2.8673,  0.5749,  0.5835, -3.5027, -2.8670, -2.5660, -2.9986, -2.0367, -0.3314,  4.3929, -3.0944, -2.0629},
          { -3.3860, -3.4526, -5.5762, -4.1081,  0.2471, -3.8385, -2.7152, -3.0862, -4.1007, -2.3768, -1.6696, -4.3518, -4.9796, -2.5874, -3.5298, -3.7697, -3.6211, -3.0944, 10.7537,  2.0515},
          { -2.3937, -3.3986, -4.0560, -3.1739,  2.9504, -4.2467,  1.7268, -1.7816, -2.5832, -1.5692, -1.7253, -2.7578, -3.9142, -2.2290, -2.5184, -2.2280, -2.1513, -2.0629,  2.0515,  7.2162}
        };
    }
}
