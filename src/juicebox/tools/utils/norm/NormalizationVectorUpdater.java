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
import javastraw.reader.DatasetReaderV2;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
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

import java.io.IOException;
import java.util.*;

/**
 * Update an existing hic file with new normalization vectors (included expected value vectors)
 *
 * @author jrobinso
 * @since 2/8/13
 */
public class NormalizationVectorUpdater extends NormVectorUpdater {

    protected List<BufferedByteWriter> normVectorBuffers = new ArrayList<>();
    protected List<NormalizationVectorIndexEntry> normVectorIndices = new ArrayList<>();
    protected List<ExpectedValueCalculation> expectedValueCalculations = new ArrayList<>();

    // Keep track of chromosomes that fail to converge, so we don't try them at higher resolutions.
    protected Set<Chromosome> scaleBPFailChroms = new HashSet<>();
    protected Set<Chromosome> scaleFragFailChroms = new HashSet<>();

    // norms to build; gets overwritten
    protected boolean weShouldBuildVC = true;
    protected boolean weShouldBuildVCSqrt = true;
    protected boolean weShouldBuildScale = true;

    protected static void printNormTiming(String norm, Chromosome chr, HiCZoom zoom, long currentTime) {
        if (HiCGlobals.printVerboseComments) {
            System.out.println(norm + " normalization of " + chr + " at " + zoom + " took " + (System.currentTimeMillis() - currentTime) + " milliseconds");
        }
    }

    protected static void updateExpectedValueCalculationForChr(final int chrIdx, NormalizationCalculations nc, ListOfFloatArrays vec, NormalizationType type, HiCZoom zoom, MatrixZoomData zd,
                                                               ExpectedValueCalculation ev, List<BufferedByteWriter> normVectorBuffers, List<NormalizationVectorIndexEntry> normVectorIndex) throws IOException {
        double factor = nc.getSumFactor(vec);
        vec.multiplyEverythingBy(factor);
        updateNormVectorIndexWithVector(normVectorIndex, normVectorBuffers, vec, chrIdx, type, zoom);
        ev.addDistancesFromIterator(chrIdx, zd.getDirectIterator(), vec);
    }

    protected void reEvaluateWhichIntraNormsToBuild(List<NormalizationType> normalizationsToBuild) {
        weShouldBuildVC = normalizationsToBuild.contains(NormalizationHandler.VC);
        weShouldBuildVCSqrt = normalizationsToBuild.contains(NormalizationHandler.VC_SQRT);
        weShouldBuildScale = normalizationsToBuild.contains(NormalizationHandler.SCALE);
    }

    public void updateHicFile(String path, List<NormalizationType> normalizationsToBuild,
                              Map<NormalizationType, Integer> resolutionsToBuildTo,
                              int genomeWideLowestResolutionAllowed, boolean noFrag) throws IOException {

        //System.out.println("test: using old norm code");
        int minResolution = Integer.MAX_VALUE;
        for (Map.Entry<NormalizationType, Integer> entry : resolutionsToBuildTo.entrySet()) {
            if (entry.getValue() < minResolution) {
                minResolution = entry.getValue();
            }
        }

        DatasetReaderV2 reader = new DatasetReaderV2(path, false);
        Dataset ds = reader.read();
        HiCGlobals.verifySupportedHiCFileWritingVersion(reader.getVersion());

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        Map<String, Integer> fragCountMap = ds.getFragmentCounts();
        List<HiCZoom> resolutions = ds.getAllPossibleResolutions();

        reEvaluateWhichIntraNormsToBuild(normalizationsToBuild);

        normVectorBuffers.add(new BufferedByteWriter());
        for (HiCZoom zoom : resolutions) {
            if (zoom.getBinSize() < minResolution) {
                System.out.println("skipping zoom" + zoom);
                continue;
            }
            if (noFrag && zoom.getUnit() == HiCZoom.HiCUnit.FRAG) continue;

            System.out.println();
            System.out.print("Calculating norms for zoom " + zoom);

            List<NormalizationType> gwNormalizations = GWNorms.getGWNorms(normalizationsToBuild, resolutionsToBuildTo, zoom);
            List<NormalizationType> interNormalizations = GWNorms.getInterNorms(normalizationsToBuild, resolutionsToBuildTo, zoom);

            Map<NormalizationType, Map<Chromosome, NormalizationVector>>
                    gwNormMaps = GWNorms.getGWNormMaps(gwNormalizations, ds, zoom, true);
            Map<NormalizationType, Map<Chromosome, NormalizationVector>>
                    interNormMaps = GWNorms.getGWNormMaps(interNormalizations, ds, zoom, false);

            Map<NormalizationType, ExpectedValueCalculation> gwMapExpected = GWNorms.createdExpectedMap(gwNormalizations,
                    interNormalizations, chromosomeHandler, zoom.getBinSize());

            ds.clearCache(true);

            Map<String, Integer> fcm = zoom.getUnit() == HiCZoom.HiCUnit.FRAG ? fragCountMap : null;

            ExpectedValueCalculation evVC = new ExpectedValueCalculation(chromosomeHandler, zoom.getBinSize(), fcm, NormalizationHandler.VC);
            ExpectedValueCalculation evVCSqrt = new ExpectedValueCalculation(chromosomeHandler, zoom.getBinSize(), fcm, NormalizationHandler.VC_SQRT);
            ExpectedValueCalculation evSCALE = new ExpectedValueCalculation(chromosomeHandler, zoom.getBinSize(), fcm, NormalizationHandler.SCALE);

            // Loop through chromosomes
            for (Chromosome chrom : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {

                MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, zoom);
                if (zd == null) continue;

                if (HiCGlobals.printVerboseComments) {
                    System.out.println("Now Doing " + chrom.getName());
                }

                BigContactArray ba = BigContactArrayCreator.createFromZD(zd);

                GWNorms.addGWNormsToBuffer(gwNormalizations, gwNormMaps, chrom, normVectorIndices,
                        normVectorBuffers, zoom, gwMapExpected, ba);
                GWNorms.addGWNormsToBuffer(interNormalizations, interNormMaps, chrom, normVectorIndices,
                        normVectorBuffers, zoom, gwMapExpected, ba);

                NormalizationCalculations nc = new NormalizationCalculations(ba, zd.getBinSize());
                zd.clearCache();

                if (weShouldBuildVC || weShouldBuildVCSqrt || weShouldBuildScale) {
                    boolean saveVC = weShouldBuildVC && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.VC);
                    boolean saveVCSqrt = weShouldBuildVCSqrt && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.VC_SQRT);
                    boolean saveScale = weShouldBuildScale && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.SCALE);

                    if (saveVC || saveVCSqrt || saveScale) {
                        buildTheNorms(saveVC, saveVCSqrt, saveScale, chrom, nc, zoom, zd, evVC, evVCSqrt, evSCALE);
                    }
                }

                zd.clearCache();
                ba.clear();
            }

            GWNorms.populateGWExpecteds(expectedValueCalculations, gwMapExpected, gwNormalizations);
            GWNorms.populateGWExpecteds(expectedValueCalculations, gwMapExpected, interNormalizations);

            if (evVC.hasData()) {
                expectedValueCalculations.add(evVC);
            }
            if (weShouldBuildVCSqrt && evVCSqrt.hasData() && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.VC_SQRT)) {
                expectedValueCalculations.add(evVCSqrt);
            }
            if (weShouldBuildScale && evSCALE.hasData() && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.SCALE)) {
                expectedValueCalculations.add(evSCALE);
            }

            ds.clearCache(false);
        }
        writeNormsToUpdateFile(reader, path, true, expectedValueCalculations, null, normVectorIndices,
                normVectorBuffers, "Finished writing norms");
    }

    private void buildTheNorms(boolean saveVC, boolean saveVCSqrt, boolean saveScale, Chromosome chr,
                               NormalizationCalculations nc, HiCZoom zoom, MatrixZoomData zd,
                               ExpectedValueCalculation evVC, ExpectedValueCalculation evVCSqrt, ExpectedValueCalculation evSCALE) throws IOException {

        final int chrIdx = chr.getIndex();
        ListOfFloatArrays vc = nc.computeVC();

        if (saveScale) {
            Set<Chromosome> failureSet = zoom.getUnit() == HiCZoom.HiCUnit.FRAG ? scaleFragFailChroms : scaleBPFailChroms;

            if (!failureSet.contains(chr)) {
                ListOfFloatArrays scale = nc.computeSCALE(vc);
                if (scale == null) {
                    failureSet.add(chr);
                } else {
                    updateExpectedValueCalculationForChr(chrIdx, nc, scale, NormalizationHandler.SCALE, zoom, zd, evSCALE, normVectorBuffers, normVectorIndices);
                }
            }
        }
        if (saveVC) {
            updateExpectedValueCalculationForChr(chrIdx, nc, vc, NormalizationHandler.VC, zoom, zd, evVC, normVectorBuffers, normVectorIndices);
        }
        if (saveVCSqrt) {
            ListOfFloatArrays vcSqrt = new ListOfFloatArrays(vc.getLength());
            for (int i = 0; i < vc.getLength(); i++) {
                vcSqrt.set(i, (float) Math.sqrt(vc.get(i)));
            }
            updateExpectedValueCalculationForChr(chrIdx, nc, vcSqrt, NormalizationHandler.VC_SQRT, zoom, zd, evVCSqrt, normVectorBuffers, normVectorIndices);
        }
    }
}
