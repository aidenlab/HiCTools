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
import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import juicebox.HiCGlobals;
import juicebox.tools.utils.bigarray.BigContactArray;
import juicebox.tools.utils.bigarray.BigContactArrayCreator;
import juicebox.tools.utils.original.ExpectedValueCalculation;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.util.Pair;

import java.io.IOException;
import java.util.*;

public class GenomeWideNormalizationVectorUpdater extends NormVectorUpdater {
    
    public static void updateHicFileForGWfromPreAddNormOnly(Dataset ds, HiCZoom zoom, List<NormalizationType> normalizationsToBuild,
                                                            Map<NormalizationType, Integer> resolutionsToBuildTo, List<NormalizationVectorIndexEntry> normVectorIndices,
                                                            List<BufferedByteWriter> normVectorBuffers, List<ExpectedValueCalculation> expectedValueCalculations) throws IOException {
        for (NormalizationType normType : normalizationsToBuild) {
            if (NormalizationHandler.isGenomeWideNorm(normType)) {
                if (zoom.getBinSize() >= resolutionsToBuildTo.get(normType)) {

                    if (HiCGlobals.printVerboseComments) {
                        System.out.println("Now Doing " + normType.getLabel());
                    }
                    long currentTime = System.currentTimeMillis();
                    Pair<Map<Chromosome, NormalizationVector>, ExpectedValueCalculation> wgVectors = getWGVectors(ds, zoom, normType);
                    if (HiCGlobals.printVerboseComments) {
                        System.out.println("\n" + normType.getLabel() + " normalization genome wide at " + zoom + " took " + (System.currentTimeMillis() - currentTime) + " milliseconds");
                    }

                    if (wgVectors != null) {
                        Map<Chromosome, NormalizationVector> nvMap = wgVectors.getFirst();
                        List<Chromosome> chromosomes = new ArrayList<>(nvMap.keySet());
                        chromosomes.sort(Comparator.comparingInt(Chromosome::getIndex));
                        for (Chromosome chromosome : chromosomes) {
                            updateNormVectorIndexWithVector(normVectorIndices, normVectorBuffers,
                                    nvMap.get(chromosome).getData().convertToFloats(), chromosome.getIndex(),
                                    normType, zoom);
                        }

                        expectedValueCalculations.add(wgVectors.getSecond());
                    }
                }
            }
        }
    }

    /**
     * Compute the whole-genome normalization and expected value vectors and return as a pair (normalization vector first)
     */

    private static Pair<Map<Chromosome, NormalizationVector>, ExpectedValueCalculation> getWGVectors(Dataset dataset,
                                                                                                     HiCZoom zoom,
                                                                                                     NormalizationType norm) {
        boolean includeIntraData = NormalizationHandler.isGenomeWideNormIntra(norm); // default INTER type
        final ChromosomeHandler chromosomeHandler = dataset.getChromosomeHandler();
        final int resolution = zoom.getBinSize();
        final BigContactArray ba = BigContactArrayCreator.createForWholeGenome(dataset, chromosomeHandler, zoom,
                includeIntraData);

        NormalizationCalculations calculations = new NormalizationCalculations(ba, resolution);
        ListOfFloatArrays vector = calculations.getNorm(norm);
        if (vector == null) {
            return null;
        }
        ba.clear();

        ExpectedValueCalculation expectedValueCalculation = new ExpectedValueCalculation(chromosomeHandler, resolution, null, norm);
        int addY = 0;
        // Loop through chromosomes
        for (Chromosome chr : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {
            MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chr, chr, zoom);
            if (zd == null) continue;
            final int chrIdx = chr.getIndex();

            Iterator<ContactRecord> iterator = zd.getDirectIterator();
            while (iterator.hasNext()) {
                ContactRecord cr = iterator.next();
                int x = cr.getBinX();
                int y = cr.getBinY();
                final float vx = vector.get(x + addY);
                final float vy = vector.get(y + addY);
                if (isValidNormValue(vx) && isValidNormValue(vy)) {
                    double value = cr.getCounts() / (vx * vy);
                    expectedValueCalculation.addDistance(chrIdx, x, y, value);
                }
            }
            zd.clearCache();
            addY += chr.getLength() / resolution + 1;
        }

        // Split normalization vector by chromosome
        Map<Chromosome, NormalizationVector> normVectorMap =
                NormalizationTools.parCreateNormVectorMap(chromosomeHandler, resolution, vector, norm, zoom);

        return new Pair<>(normVectorMap, expectedValueCalculation);
    }
}
