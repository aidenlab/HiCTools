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

import hic.tools.utils.bigarray.BigContactArray;
import hic.tools.utils.bigarray.BigContactList;
import hic.tools.utils.bigarray.BigGWContactArray;
import hic.tools.utils.bigarray.BigGWContactArrayCreator;
import hic.tools.utils.largelists.BigListOfByteWriters;
import hic.tools.utils.original.ExpectedValueCalculation;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.io.IOException;
import java.util.*;

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

    public static Map<NormalizationType, Map<Chromosome, NormalizationVector>> getGWNormMaps(List<NormalizationType> gwNorms,
                                                                                             List<NormalizationType> interNorms,
                                                                                             Dataset ds, HiCZoom zoom) {

        if (gwNorms.isEmpty() && interNorms.isEmpty()) return new HashMap<>();

        final ChromosomeHandler handler = ds.getChromosomeHandler();
        final BigGWContactArray ba;
        if (gwNorms.isEmpty()) {
            // only INTER_ norms getting used
            ba = BigGWContactArrayCreator.createForWholeGenome(ds, handler, zoom, false);
        } else if (interNorms.isEmpty()) {
            // only GW_ norms getting used
            ba = BigGWContactArrayCreator.createForWholeGenome(ds, handler, zoom, true);
        } else {
            // using both
            ba = BigGWContactArrayCreator.createForWholeGenome(ds, handler, zoom, true);
        }

        Map<NormalizationType, Map<Chromosome, NormalizationVector>> result = new HashMap<>();
        for (NormalizationType normType : gwNorms) {
            Map<Chromosome, NormalizationVector> wgVectors = getWGVectors(zoom, normType, ba, handler, "GW");
            if (wgVectors != null) {
                result.put(normType, wgVectors);
            }
        }

        ba.clearIntraAndShiftInter();

        for (NormalizationType normType : interNorms) {
            Map<Chromosome, NormalizationVector> wgVectors = getWGVectors(zoom, normType, ba, handler, "INTER");
            if (wgVectors != null) {
                result.put(normType, wgVectors);
            }
        }

        ba.clearAll();
        return result;
    }

    private static Map<Chromosome, NormalizationVector> getWGVectors(HiCZoom zoom, NormalizationType norm,
                                                                     BigContactArray ba, ChromosomeHandler handler,
                                                                     String stem) {
        final int resolution = zoom.getBinSize();
        NormalizationCalculations calculations = new NormalizationCalculations(ba, resolution);
        ListOfFloatArrays vector = calculations.getNorm(norm, stem + "_NORM_" + zoom.getBinSize());
        if (vector == null) {
            return null;
        }
        return NormalizationTools.parCreateNormVectorMap(handler, resolution, vector, norm, zoom);
    }

    public static void addGWNormsToBuffer(List<NormalizationType> norms,
                                          Map<NormalizationType, Map<Chromosome, NormalizationVector>> normMap,
                                          Chromosome chromosome, List<NormalizationVectorIndexEntry> normVectorIndices,
                                          BigListOfByteWriters normVectorBuffers, HiCZoom zoom,
                                          Map<NormalizationType, ExpectedValueCalculation> expectedMap,
                                          BigContactList ba) throws IOException {
        for (NormalizationType norm : norms) {
            if (normMap.containsKey(norm)) {
                Map<Chromosome, NormalizationVector> map = normMap.get(norm);
                if (map.containsKey(chromosome)) {
                    ListOfFloatArrays vector = map.get(chromosome).getData().convertToFloats();
                    NormVectorUpdater.updateNormVectorIndexWithVector(normVectorIndices, normVectorBuffers,
                            vector, chromosome.getIndex(),
                            norm, zoom);

                    if (expectedMap.containsKey(norm)) {
                        ExpectedValueCalculation exp = expectedMap.get(norm);
                        final int chrIdx = chromosome.getIndex();
                        ba.updateGenomeWideExpected(chrIdx, vector, exp);
                    }
                }
            }
        }
    }

    public static Map<NormalizationType, ExpectedValueCalculation> createdExpectedMap(
            List<NormalizationType> gwNorms, List<NormalizationType> interNorms,
            ChromosomeHandler chromosomeHandler, int resolution) {

        Map<NormalizationType, ExpectedValueCalculation> expMap = new HashMap<>();

        for (NormalizationType norm : gwNorms) {
            expMap.put(norm, new ExpectedValueCalculation(chromosomeHandler, resolution, norm));
        }

        for (NormalizationType norm : interNorms) {
            expMap.put(norm, new ExpectedValueCalculation(chromosomeHandler, resolution, norm));
        }

        return expMap;
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
