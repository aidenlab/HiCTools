/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
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

package juicebox.tools.utils.norm;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import juicebox.tools.utils.bigarray.BigContactArray;
import juicebox.tools.utils.bigarray.BigContactArrayCreator;
import juicebox.tools.utils.original.ExpectedValueCalculation;
import org.broad.igv.tdf.BufferedByteWriter;

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

    public static Map<NormalizationType, Map<Chromosome, NormalizationVector>> getGWNormMaps(List<NormalizationType> norms,
                                                                                             Dataset ds, HiCZoom zoom,
                                                                                             boolean includeIntra) {

        if (norms.isEmpty()) return new HashMap<>();
        final ChromosomeHandler handler = ds.getChromosomeHandler();
        final BigContactArray ba = BigContactArrayCreator.createForWholeGenome(ds, handler, zoom, includeIntra);

        Map<NormalizationType, Map<Chromosome, NormalizationVector>> result = new HashMap<>();
        for (NormalizationType normType : norms) {
            Map<Chromosome, NormalizationVector> wgVectors = getWGVectors(zoom, normType, ba, handler);
            if (wgVectors != null) {
                result.put(normType, wgVectors);
            }
        }
        ba.clear();
        return result;
    }

    private static Map<Chromosome, NormalizationVector> getWGVectors(HiCZoom zoom, NormalizationType norm,
                                                                     BigContactArray ba, ChromosomeHandler handler) {
        final int resolution = zoom.getBinSize();
        NormalizationCalculations calculations = new NormalizationCalculations(ba, resolution);
        ListOfFloatArrays vector = calculations.getNorm(norm);
        if (vector == null) {
            return null;
        }
        return NormalizationTools.parCreateNormVectorMap(handler, resolution, vector, norm, zoom);
    }

    public static void addGWNormsToBuffer(List<NormalizationType> norms,
                                          Map<NormalizationType, Map<Chromosome, NormalizationVector>> normMap,
                                          Chromosome chromosome, List<NormalizationVectorIndexEntry> normVectorIndices,
                                          List<BufferedByteWriter> normVectorBuffers, HiCZoom zoom,
                                          Map<NormalizationType, ExpectedValueCalculation> expectedMap,
                                          BigContactArray ba) throws IOException {
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
            expMap.put(norm, new ExpectedValueCalculation(chromosomeHandler, resolution, null, norm));
        }

        for (NormalizationType norm : interNorms) {
            expMap.put(norm, new ExpectedValueCalculation(chromosomeHandler, resolution, null, norm));
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
