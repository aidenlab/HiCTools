/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2020-2022 Rice University, Baylor College of Medicine, Aiden Lab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package hic.tools.utils.norm;

import hic.tools.utils.bigarray.BigContactList;
import hic.tools.utils.bigarray.BigGWContactArrayCreator;
import hic.tools.utils.original.ExpectedValueCalculation;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class GWNorms {
    public static List<NormalizationType> getGWNorms(List<NormalizationType> allNorms,
                                                     Map<NormalizationType, Integer> resolutionsToBuildTo, HiCZoom zoom) {
        List<NormalizationType> norms = new ArrayList<>();
        for (NormalizationType norm : allNorms) {
            if (NormalizationHandler.isGenomeWideNorm(norm) && NormalizationHandler.isGenomeWideNormIntra(norm)) {
                if (zoom.getBinSize() >= resolutionsToBuildTo.get(norm)) {
                    norms.add(norm);
                }
            }
        }
        norms.sort(Comparator.comparing(NormalizationType::getLabel));
        return norms;
    }

    public static List<NormalizationType> getInterNorms(List<NormalizationType> allNorms,
                                                        Map<NormalizationType, Integer> resolutionsToBuildTo, HiCZoom zoom) {
        List<NormalizationType> norms = new ArrayList<>();
        for (NormalizationType norm : allNorms) {
            if (NormalizationHandler.isGenomeWideNorm(norm) && !NormalizationHandler.isGenomeWideNormIntra(norm)) {
                if (zoom.getBinSize() >= resolutionsToBuildTo.get(norm)) {
                    norms.add(norm);
                }
            }
        }
        norms.sort(Comparator.comparing(NormalizationType::getLabel));
        return norms;
    }

    public static void getGWNormMaps(Dataset ds, HiCZoom zoom, int resCutoffForRAM, NormVectorsContainer container) {

        if (container.hasNoGenomewideNorms()) return;

        final ChromosomeHandler handler = ds.getChromosomeHandler();
        final BigContactList ba;
        if (zoom.getBinSize() < 25 * resCutoffForRAM) {
            ba = BigGWContactArrayCreator.createLocalVersionWholeGenome(ds, handler, zoom, container.useGWIntra());
        } else {
            ba = BigGWContactArrayCreator.createForWholeGenome(ds, handler, zoom, container.useGWIntra());
        }

        for (NormalizationType normType : container.getGenomewideNorms()) {
            Map<Chromosome, FloatNormVector> wgVectors = getWGVectors(zoom, normType, ba, handler, "GW");
            if (wgVectors != null) {
                container.put(normType, wgVectors);
            }
        }

        ba.clearIntraAndShiftInter();

        for (NormalizationType normType : container.getGenomewideInterNorms()) {
            Map<Chromosome, FloatNormVector> wgVectors = getWGVectors(zoom, normType, ba, handler, "INTER");
            if (wgVectors != null) {
                container.put(normType, wgVectors);
            }
        }

        ba.clear();
    }

    private static Map<Chromosome, FloatNormVector> getWGVectors(HiCZoom zoom, NormalizationType norm,
                                                                 BigContactList ba, ChromosomeHandler handler,
                                                                 String stem) {
        final int resolution = zoom.getBinSize();
        NormalizationCalculations calculations = new NormalizationCalculations(ba, resolution);
        ListOfFloatArrays vector = calculations.getNormWithFix(norm, stem + "_NORM_" + zoom.getBinSize());
        if (vector == null) {
            return null;
        }
        return NormalizationTools.parCreateNormVectorMap(handler, resolution, vector, norm, zoom);
    }

    public static void populateGWExpecteds(List<ExpectedValueCalculation> expectedValueCalculations,
                                           Map<NormalizationType, ExpectedValueCalculation> expectedMap,
                                           List<NormalizationType> norms) {
        for (NormalizationType norm : norms) {
            if (expectedMap.containsKey(norm)) {
                ExpectedValueCalculation expected = expectedMap.get(norm);
                if (expected.hasData()) {
                    expectedValueCalculations.add(expected);
                }
            }
        }
    }
}
