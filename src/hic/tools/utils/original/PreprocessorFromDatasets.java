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

package hic.tools.utils.original;

import hic.tools.utils.iterators.contacts.Contact;
import hic.tools.utils.iterators.contacts.ContactIterator;
import hic.tools.utils.merge.HiCMergeTools;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.zip.Deflater;


public class PreprocessorFromDatasets extends HiCFileBuilder {

    public PreprocessorFromDatasets(File outputFile, String genomeId, ChromosomeHandler chromosomeHandler,
                                    double hicFileScalingFactor) {
        super(outputFile, genomeId, chromosomeHandler, hicFileScalingFactor);
    }

    public void preprocess(List<Dataset> inputDS, final String headerFile, final String footerFile) throws IOException {

        HiCMergeTools.mergeStatsAndGraphs(inputDS, tmpDir, this);

        try {
            LittleEndianOutputStream[] losFooter = initializeLosArrays(headerFile, footerFile);
            writeHeader();
            writeBody(inputDS);
            writeFooter(losFooter);
            closeLosArray(losFooter);
        } finally {
            closeLosArray(losArray);
        }

        updateMasterIndex(headerFile);
        System.out.println("\nFinished preprocess");
    }

    private MatrixPP computeWholeGenomeMatrix(List<Dataset> inputDS) throws IOException {
        MatrixPP matrix = getInitialGenomeWideMatrixPP(chromosomeHandler);
        ContactIterator iter = null;
        try {
            iter = ContactIterator.getAllByAllIterator(inputDS);
            while (iter.hasNext()) {
                Contact contact = iter.next();
                matrix.incrementCount(contact.getPos1(), contact.getPos2(), contact.getScore(), expectedValueCalculations, tmpDir);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(98);
        }
        matrix.parsingComplete();
        return matrix;
    }

    private void writeBody(List<Dataset> inputDS) throws IOException {
        System.out.println("Writing body");
        writeWholeGenomeMatrix(inputDS, losArray, compressor, matrixPositions);

        Chromosome[] chromosomes = chromosomeHandler.getChromosomeArrayWithoutAllByAll();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                int currentChr1 = chromosomes[i].getIndex();
                int currentChr2 = chromosomes[j].getIndex();
                writeChromosomeRegionMatrix(currentChr1, currentChr2, inputDS);
            }
        }

        masterIndexPosition = losArray[0].getWrittenCount();
    }

    private void writeChromosomeRegionMatrix(int chrom1Index, int chrom2Index, List<Dataset> inputDS) throws IOException {
        int newBlockCapacity = (int) (Math.sqrt(inputDS.size()) * BLOCK_CAPACITY);
        MatrixPP mergedMatrix = new MatrixPP(chrom1Index, chrom2Index, chromosomeHandler, bpBinSizes,
                countThreshold, v9DepthBase, newBlockCapacity);

/*
        PairIterator iter = PairIterator.getDatasetIterator(inputDS, chromosomeHandler);
        while (iter.hasNext()) {
            AlignmentPair pair = iter.next();
            mergedMatrix.incrementCount(bp1, bp2, score, expectedValueCalculations, tmpDir);
        }
        iter.close();



        mergedMatrix.mergeMatrices();

 */

        mergedMatrix.parsingComplete();
        writeMatrix(mergedMatrix, losArray, compressor, matrixPositions, -1, false);
        mergedMatrix = null;
    }

    private void writeWholeGenomeMatrix(List<Dataset> inputDS, LittleEndianOutputStream[] losArray, Deflater compressor,
                                        Map<String, IndexEntry> matrixPositions) throws IOException {
        MatrixPP wholeGenomeMatrix = computeWholeGenomeMatrix(inputDS);
        writeMatrix(wholeGenomeMatrix, losArray, compressor, matrixPositions, -1, false);
        wholeGenomeMatrix = null;
    }
}
