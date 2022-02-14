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
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import juicebox.HiCGlobals;
import juicebox.tools.utils.original.ExpectedValueCalculation;
import org.broad.igv.tdf.BufferedByteWriter;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class CustomNormVectorFileHandler extends NormVectorUpdater {

    public static void updateHicFile(String path, String vectorPath) throws IOException {
        DatasetReaderV2 reader = new DatasetReaderV2(path, false);
        Dataset ds = reader.read();
        HiCGlobals.verifySupportedHiCFileVersion(reader.getVersion());

        String[] vectorPaths = vectorPath.split(",");
        NormVectorInfo normVectorInfo = completeCalculationsNecessaryForUpdatingCustomNormalizations(ds, vectorPaths, true);
        writeNormsToUpdateFile(reader, path, false, null, normVectorInfo.getExpectedValueFunctionMap(),
                normVectorInfo.getNormVectorIndices(), normVectorInfo.getNormVectorBuffers(), "Finished adding another normalization.");

        System.out.println("all custom norms added");
    }

    private static NormVectorInfo completeCalculationsNecessaryForUpdatingCustomNormalizations(
            final Dataset ds, String[] filePaths, boolean overwriteHicFileFooter) throws IOException {

        Map<NormalizationType, Map<String, NormalizationVector>> normalizationVectorMap = readVectorFile(filePaths,
                ds.getChromosomeHandler(), ds.getNormalizationHandler());

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        Map<String, Integer> fragCountMap = ds.getFragmentCounts();
        List<HiCZoom> resolutions = ds.getAllPossibleResolutions();

        List<BufferedByteWriter> normVectorBuffers = new ArrayList<>();
        List<NormalizationVectorIndexEntry> normVectorIndices = new ArrayList<>();
        Map<String, ExpectedValueFunction> expectedValueFunctionMap = ds.getExpectedValueFunctionMap();

        expectedValueFunctionMap.entrySet().removeIf(entry -> entry.getKey().contains("NONE"));

        // Get existing norm vectors so we don't lose them
        if (overwriteHicFileFooter) {
            for (HiCZoom zoom : resolutions) {
                for (NormalizationType type : NormalizationHandler.getAllNormTypes()) {
                    for (Chromosome chr : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {
                        NormalizationVector existingNorm = ds.getNormalizationVector(chr.getIndex(), zoom, type);
                        if (existingNorm != null) {
                            long position = 0;
                            for (int i=0; i < normVectorBuffers.size(); i++) {
                                position += normVectorBuffers.get(i).bytesWritten();
                            }
                            // todo @suhas
                            putFloatArraysIntoBufferList(normVectorBuffers, existingNorm.getData().convertToFloats().getValues());

                            long newPos = 0;
                            for (int i=0; i < normVectorBuffers.size(); i++) {
                                newPos += normVectorBuffers.get(i).bytesWritten();
                            }
                            int sizeInBytes = (int) (newPos - position);
                            normVectorIndices.add(new NormalizationVectorIndexEntry(
                                    type.toString(), chr.getIndex(), zoom.getUnit().toString(), zoom.getBinSize(), position, sizeInBytes));
                        }
                    }
                }
            }
        }

        //ExecutorService executor = HiCGlobals.newFixedThreadPool();
        for (NormalizationType customNormType : normalizationVectorMap.keySet()) {
            final Map<String, NormalizationVector> normVectorsByChrAndZoom = normalizationVectorMap.get(customNormType);
            final Set<String> keySet = new HashSet<>(normVectorsByChrAndZoom.keySet());
            final Map<Integer, Integer> chrAndResolutionWhichFailed = new HashMap<>();

            for (final String key : keySet) {
                final NormalizationVector nv = normVectorsByChrAndZoom.get(key);
                if (chrAndResolutionWhichFailed.containsKey(nv.getChrIdx()) && nv.getResolution() < chrAndResolutionWhichFailed.get(nv.getChrIdx())) {
                    normVectorsByChrAndZoom.remove(key);
                    continue;
                }
            }
        }

        /*
        executor.shutdown();
        // Wait until all threads finish
        while (!executor.isTerminated()) {
        }
        */

        for (HiCZoom zoom : resolutions) {
            Map<String, Integer> fcm = zoom.getUnit() == HiCZoom.HiCUnit.FRAG ? fragCountMap : null;

            for (NormalizationType customNormType : normalizationVectorMap.keySet()) {

                ExpectedValueCalculation evLoaded = new ExpectedValueCalculation(chromosomeHandler, zoom.getBinSize(), fcm, customNormType);
                String key = ExpectedValueFunction.getKey(zoom, customNormType, false, 0);

                // Loop through chromosomes
                for (Chromosome chr : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {
                    MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr, chr, zoom);
                    if (zd == null) continue;

                    handleLoadedVector(customNormType, chr.getIndex(), zoom, normalizationVectorMap.get(customNormType),
                                normVectorBuffers, normVectorIndices, zd, evLoaded);
                }
                expectedValueFunctionMap.put(key, evLoaded.getExpectedValueFunction());
            }
        }

        ds.setExpectedValueFunctionMap(expectedValueFunctionMap);
        return new NormVectorInfo(normalizationVectorMap, normVectorBuffers, normVectorIndices, expectedValueFunctionMap);
    }

    private static void handleLoadedVector(NormalizationType customNormType, final int chrIndx, HiCZoom zoom, Map<String, NormalizationVector> normVectors,
                                           List<BufferedByteWriter> normVectorBuffers, List<NormalizationVectorIndexEntry> normVectorIndex,
                                           MatrixZoomData zd, ExpectedValueCalculation evLoaded) throws IOException {

        String key = NormalizationVector.getKey(customNormType, chrIndx, zoom.getUnit().toString(), zoom.getBinSize());
        if (normVectors.containsKey(key)) {
            NormalizationVector vector = normVectors.get(key);
            if (vector == null || vector.getData() == null) return;
            // Write custom norm
            long position = 0;
            for (int i=0; i < normVectorBuffers.size(); i++) {
                position += normVectorBuffers.get(i).bytesWritten();
            }
            // todo @suhas
            putFloatArraysIntoBufferList(normVectorBuffers, vector.getData().convertToFloats().getValues());

            long newPos = 0;
            for (int i=0; i < normVectorBuffers.size(); i++) {
                newPos += normVectorBuffers.get(i).bytesWritten();
            }

            int sizeInBytes = (int) (newPos - position);
            normVectorIndex.add(new NormalizationVectorIndexEntry(
                    customNormType.toString(), chrIndx, zoom.getUnit().toString(), zoom.getBinSize(), position, sizeInBytes));
    
            evLoaded.addDistancesFromIterator(chrIndx, zd.getDirectIterator(), vector.getData().convertToFloats());
        }
    }

    private static Map<NormalizationType, Map<String, NormalizationVector>> readVectorFile(String[] fnames, ChromosomeHandler chromosomeHandler, NormalizationHandler normalizationHandler) throws IOException {

        Map<NormalizationType, Map<String, NormalizationVector>> normVectors = new HashMap<>();

        for (String fname : fnames) {
            BufferedReader vectorReader;
            if (fname.endsWith(".gz")) {
                InputStream fileStream = new FileInputStream(fname);
                InputStream gzipStream = new GZIPInputStream(fileStream);
                Reader decoder = new InputStreamReader(gzipStream, StandardCharsets.UTF_8);
                vectorReader = new BufferedReader(decoder, 4194304);
            } else {
                //this.reader = org.broad.igv.util.ParsingUtils.openBufferedReader(path);
                vectorReader = new BufferedReader(new InputStreamReader(new FileInputStream(fname)), HiCGlobals.bufferSize);
            }

            Chromosome chr = null;
            int resolution = -1;
            HiCZoom.HiCUnit unit = null;
            NormalizationType customNormType = null;
            boolean needsToBeScaledTo = false;

            String nextLine = vectorReader.readLine();
            while (nextLine != null) {
                // Header: vector  type  chr1    2048000 BP
                if (nextLine.startsWith("vector")) {
                    String[] tokens = nextLine.split("\\s+");
                    chr = chromosomeHandler.getChromosomeFromName(tokens[2]);
                    if (chr == null) {
                        System.err.println("Skipping " + tokens[2] + " which isn't in dataset");
                        nextLine = skipLinesUntilTextEncountered(vectorReader, "vector");
                        continue;
                    }

                    customNormType = normalizationHandler.getNormTypeFromString(tokens[1]);
                    resolution = Integer.parseInt(tokens[3]);
                    unit = HiCZoom.HiCUnit.valueOf(tokens[4]);
                    needsToBeScaledTo = tokens[0].toLowerCase().contains("scale");
                }
                if (chr != null && customNormType != null) {
                    if (HiCGlobals.printVerboseComments) {
                        System.out.println("Adding norm " + customNormType + " for chr " + chr.getName() + " at " + resolution + " " + unit + " resolution.");
                    }
    
                    // Now do work on loaded norm vector
                    // Create the new vector by looping through the loaded vector file line by line
                    // assume custom norm vectors aren't for indices requiring long
                    long size = (chr.getLength() / resolution + 1);
                    ListOfDoubleArrays data = new ListOfDoubleArrays(size);
                    int i = 0;
                    nextLine = vectorReader.readLine();
                    // List<Double> data = new ArrayList<Double>();
                    while (nextLine != null && !(nextLine.startsWith("vector"))) {
                        if (nextLine.equalsIgnoreCase("nan") || nextLine.equals(".")) {
                            data.set(i, Double.NaN);
                        } else {
                            data.set(i, Double.parseDouble(nextLine));
                        }
                        i++;
                        if (i > size) {
                            throw new IOException("More values than resolution would indicate");
                        }
                        nextLine = vectorReader.readLine();
                    }

                    if (!normVectors.containsKey(customNormType)) {
                        normVectors.put(customNormType, new HashMap<>());
                    }
                    NormalizationVector vector = new NormalizationVector(customNormType, chr.getIndex(), unit,
                            resolution, data);
                    normVectors.get(customNormType).put(vector.getKey(), vector);

                } else {
                    System.err.println("Chromosome vector null"); // this shouldn't happen due to continue above
                }
            }
        }

        return normVectors;
    }

    private static String skipLinesUntilTextEncountered(BufferedReader vectorReader, String string) throws IOException {
        String nextLine = vectorReader.readLine();
        while (nextLine != null && !(nextLine.startsWith(string))) {
            nextLine = vectorReader.readLine();
        }
        return nextLine;
    }
}
