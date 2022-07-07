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

import hic.HiCGlobals;
import hic.tools.utils.bigarray.BigContactList;
import hic.tools.utils.largelists.BigListOfByteWriters;
import hic.tools.utils.original.ExpectedValueCalculation;
import javastraw.reader.Dataset;
import javastraw.reader.DatasetReaderV2;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.io.IOException;
import java.util.*;

/**
 * Update an existing hic file with new normalization vectors (included expected value vectors)
 *
 * @author jrobinso
 * @since 2/8/13
 */
public class NormalizationVectorUpdater extends NormVectorUpdater {

    // Keep track of chromosomes that fail to converge, so we don't try them at higher resolutions.
    protected final Set<Chromosome> scaleBPFailChroms = new HashSet<>();

    // norms to build; gets overwritten
    protected boolean weShouldBuildVC = true;
    protected boolean weShouldBuildVCSqrt = true;
    protected boolean weShouldBuildScale = true;

    protected void reEvaluateWhichIntraNormsToBuild(List<NormalizationType> normalizationsToBuild) {
        weShouldBuildVC = normalizationsToBuild.contains(NormalizationHandler.VC);
        weShouldBuildVCSqrt = normalizationsToBuild.contains(NormalizationHandler.VC_SQRT);
        weShouldBuildScale = normalizationsToBuild.contains(NormalizationHandler.SCALE);
    }

    public void updateHicFile(String path, List<NormalizationType> normalizationsToBuild,
                              Map<NormalizationType, Integer> resolutionsToBuildTo,
                              int resolutionCutoffToSaveRAM) throws IOException {

        reEvaluateWhichIntraNormsToBuild(normalizationsToBuild);

        Map<Integer, NormVectorsContainer> containers = getAllTheNormVectors(path, normalizationsToBuild,
                resolutionsToBuildTo, resolutionCutoffToSaveRAM);

        System.gc();

        calculateExpectedsAndWriteToFile(path, containers, resolutionCutoffToSaveRAM);
    }


    public Map<Integer, NormVectorsContainer> getAllTheNormVectors(String path, List<NormalizationType> normalizationsToBuild,
                                                                   Map<NormalizationType, Integer> resolutionsToBuildTo,
                                                                   int resolutionCutoffToSaveRAM) throws IOException {

        //System.out.println("test: using old norm code");
        int minResolution = Integer.MAX_VALUE;
        for (Map.Entry<NormalizationType, Integer> entry : resolutionsToBuildTo.entrySet()) {
            if (entry.getValue() < minResolution) {
                minResolution = entry.getValue();
            }
        }

        Map<Integer, NormVectorsContainer> containers = new HashMap<>();

        DatasetReaderV2 reader = new DatasetReaderV2(path, false, false);
        Dataset ds = reader.read();
        HiCGlobals.verifySupportedHiCFileWritingVersion(reader.getVersion());

        List<HiCZoom> resolutions = ds.getAllPossibleResolutions();

        boolean interDataAvailable = NormalizationTools.checkIfInterDataAvailable(resolutions.get(0), ds);
        if (!interDataAvailable) {
            System.out.println("Skipping GW_* and INTER_* normalizations because these regions " +
                    "have no data in this .hic file");
        }

        for (HiCZoom zoom : resolutions) {
            if (zoom.getBinSize() < minResolution) {
                System.out.println("Skipping zoom" + zoom);
                continue;
            }
            if (zoom.getUnit() == HiCZoom.HiCUnit.FRAG) continue;

            System.out.println();
            System.out.print("Calculating norms for zoom " + zoom);

            NormVectorsContainer container = new NormVectorsContainer(normalizationsToBuild, resolutionsToBuildTo, zoom);

            if (interDataAvailable) {
                GWNorms.getGWNormMaps(ds, zoom, resolutionCutoffToSaveRAM, container);
            }
            ds.clearCache(true);

            IntraNorms.getAllTheNorms(ds, zoom, resolutionCutoffToSaveRAM, container,
                    weShouldBuildVC, weShouldBuildVCSqrt, weShouldBuildScale, resolutionsToBuildTo, scaleBPFailChroms);

        }
        ds.clearCache(false);
        ds = null;
        System.out.println("Balancing calculations completed");
        return containers;
    }

    public void calculateExpectedsAndWriteToFile(String path, Map<Integer, NormVectorsContainer> containers,
                                                 int resolutionCutoffToSaveRAM) throws IOException {

        DatasetReaderV2 reader = new DatasetReaderV2(path, false, false);
        Dataset ds = reader.read();
        ds.clearCache(true);

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        List<HiCZoom> resolutions = ds.getAllPossibleResolutions();

        final BigListOfByteWriters finalNormVectorBuffers = new BigListOfByteWriters();
        final List<NormalizationVectorIndexEntry> finalNormVectorIndices = new ArrayList<>();
        final List<ExpectedValueCalculation> finalExpectedValueCalculations = new ArrayList<>();
        final List<NormalizationType> sortedNorms = NormVectorsContainer.sortedNorms();

        for (HiCZoom zoom : resolutions) {
            finalNormVectorBuffers.expandBuffer();
            if (zoom.getUnit() == HiCZoom.HiCUnit.FRAG) continue;

            int resolution = zoom.getBinSize();
            if (containers.containsKey(resolution)) {
                NormVectorsContainer container = containers.get(resolution);
                if (container == null) continue;

                // Loop through chromosomes
                for (Chromosome chrom : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {

                    BigContactList ba = IntraNorms.getBigArrayFromAndClearCache(ds, chrom, zoom, resolutionCutoffToSaveRAM);
                    if (ba == null) continue;

                    for (NormalizationType norm : sortedNorms) {
                        if (container.containsNorm(norm)) {
                            Map<Chromosome, FloatNormVector> map = container.get(norm);
                            if (map.containsKey(chrom)) {
                                ListOfFloatArrays vector = map.get(chrom).getData();
                                NormVectorUpdater.updateNormVectorIndexWithVector(finalNormVectorIndices,
                                        finalNormVectorBuffers,
                                        vector, chrom.getIndex(), norm, zoom);

                                ExpectedValueCalculation expected = new ExpectedValueCalculation(chromosomeHandler, zoom.getBinSize(), norm);
                                ba.updateGenomeWideExpected(chrom.getIndex(), vector, expected);
                                finalExpectedValueCalculations.add(expected);
                            }
                        }
                    }

                    ba.clear();
                }
                container.clear();
                containers.remove(resolution);
            }
        }

        ds.clearCache(false);
        writeNormsToUpdateFile(reader, path, true, finalExpectedValueCalculations,
                null, finalNormVectorIndices,
                finalNormVectorBuffers, "Finished writing norms");
    }
}
