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

package hic.tools.utils.original;

import hic.tools.utils.cleaner.ContactCleaner;
import hic.tools.utils.iterators.mnd.AlignmentPair;
import hic.tools.utils.iterators.mnd.PairIterator;
import htsjdk.tribble.util.LittleEndianOutputStream;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class Preprocessor extends HiCFileBuilder {


    public Preprocessor(File outputFile, String genomeId, double hicFileScalingFactor, String tmp) {
        super(outputFile, genomeId, hicFileScalingFactor, tmp);
    }

    public void preprocess(final String inputFile, final String headerFile, final String footerFile,
                           Map<Integer, List<Chunk>> mndIndex) throws IOException {
        File file = new File(inputFile);
        if (!file.exists() || file.length() == 0) {
            System.err.println(inputFile + " does not exist or does not contain any reads.");
            System.exit(57);
        }

        try {
            LittleEndianOutputStream[] losFooter = initializeLosArrays(headerFile, footerFile);
            writeHeader();
            writeBody(inputFile, mndIndex);
            writeFooter(losFooter);
            closeLosArray(losFooter);
        } finally {
            closeLosArray(losArray);
        }

        updateMasterIndex(headerFile);
        System.out.println("\nFinished preprocess");
    }

    /**
     * @param file List of files to read
     * @return Matrix with counts in each bin
     */
    private MatrixPP computeWholeGenomeMatrix(String file) throws IOException {

        MatrixPP matrix = getInitialGenomeWideMatrixPP(chromosomeHandler);

        PairIterator iter = null;

        // int belowMapq = 0;
        // int intraFrag = 0;
        // int totalRead = 0;
        // int contig = 0;
        // int hicContact = 0;

        // Create an index the first time through
        try {
            iter = PairIterator.getIterator(file, chromosomeIndexes, chromosomeHandler);
            //ContactFilter filter = getContactFilter();

            while (iter.hasNext()) {
                // totalRead++;
                AlignmentPair pair = iter.next();
                if (pair.isNotContigPair()) {
                    int bp1 = pair.getPos1();
                    int bp2 = pair.getPos2();
                    int chr1 = pair.getChr1();
                    int chr2 = pair.getChr2();

                    int pos1, pos2;
                    if (shouldSkipContact(pair)) continue;
                    pos1 = ContactCleaner.getWholeGenomePosition(chr1, bp1, chromosomeHandler);
                    pos2 = ContactCleaner.getWholeGenomePosition(chr2, bp2, chromosomeHandler);
                    matrix.incrementCount(pos1, pos2, pair.getScore(), expectedValueCalculations, tmpDir);
                    // hicContact++;
                }
            }
        } finally {
            if (iter != null) iter.close();
        }
        matrix.parsingComplete();
        return matrix;
    }



    protected void writeBody(String inputFile, Map<Integer, List<Chunk>> mndIndex) throws IOException {
        System.out.println("Writing body");
        MatrixPP wholeGenomeMatrix = computeWholeGenomeMatrix(inputFile);
        writeMatrix(wholeGenomeMatrix, losArray, compressor, matrixPositions,
                -1, false, outputFile);

        PairIterator iter = PairIterator.getIterator(inputFile, chromosomeIndexes, chromosomeHandler);

        Set<String> writtenMatrices = Collections.synchronizedSet(new HashSet<>());

        int currentChr1 = -1;
        int currentChr2 = -1;
        MatrixPP currentMatrix = null;
        String currentMatrixKey = null;
        ContactCleaner cleaner = new ContactCleaner(chromosomeHandler);

        while (iter.hasNext()) {
            AlignmentPair pair = iter.next();
            // skip pairs that mapped to contigs
            if (pair.isNotContigPair()) {
                if (shouldSkipContact(pair)) continue;
                // Flip pair if needed so chr1 < chr2
                cleaner.updateLatestContact(pair);

                // Randomize position within fragment site

                // only increment if not intraFragment and passes the mapq threshold
                if (cleaner.doesntMatchCurrentBlock(currentChr1, currentChr2)) {
                    // Starting a new matrix
                    if (currentMatrix != null) {
                        currentMatrix.parsingComplete();
                        writeMatrix(currentMatrix, losArray, compressor, matrixPositions,
                                -1, false, outputFile);
                        writtenMatrices.add(currentMatrixKey);
                        currentMatrix = null;
                        System.gc();
                        //System.out.println("Available memory: " + RuntimeUtils.getAvailableMemory());
                    }

                    // Start the next matrix
                    currentChr1 = cleaner.getChr1();
                    currentChr2 = cleaner.getChr2();
                    currentMatrixKey = currentChr1 + "_" + currentChr2;

                    if (writtenMatrices.contains(currentMatrixKey)) {
                        System.err.println("Error: the chromosome combination " + currentMatrixKey + " appears in multiple blocks");
                        if (outputFile != null) outputFile.deleteOnExit();
                        System.exit(58);
                    }
                    currentMatrix = new MatrixPP(currentChr1, currentChr2, chromosomeHandler, bpBinSizes,
                            countThreshold, v9DepthBase, BLOCK_CAPACITY);
                }
                cleaner.incrementCount(currentMatrix, expectedValueCalculations, tmpDir);
            }
        }

        if (currentMatrix != null) {
            currentMatrix.parsingComplete();
            writeMatrix(currentMatrix, losArray, compressor, matrixPositions,
                    -1, false, outputFile);
        }

        iter.close();

        masterIndexPosition = losArray[0].getWrittenCount();
    }

    protected boolean shouldSkipContact(AlignmentPair pair) { // static
        int chr1 = pair.getChr1();
        int chr2 = pair.getChr2();
        if (intraChromosomalOnly && chr1 != chr2) return true;
        if (onlyNearDiagonalContacts && tooFarFromDiagonal(pair.getPos1(), pair.getPos2())) return true;
        if (includedChromosomes != null && chr1 != 0) {
            String c1Name = chromosomeHandler.getChromosomeFromIndex(chr1).getName();
            String c2Name = chromosomeHandler.getChromosomeFromIndex(chr2).getName();
            if (!includedChromosomes.contains(c1Name) || !includedChromosomes.contains(c2Name)) {
                return true;
            }
        }
        if (filter != null) {
            if (!filter.pairTypesAreEqual(pair)) {
                return true;
            }
        }
        int mapq = Math.min(pair.getMapq1(), pair.getMapq2());
        return mapq < mapqThreshold;
    }
}
