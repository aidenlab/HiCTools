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

package hic.tools.clt.old;

import hic.HiCGlobals;
import hic.tools.clt.CommandLineParser;
import hic.tools.clt.JuiceboxCLT;
import hic.tools.utils.ShellCommandRunner;
import hic.tools.utils.original.MultithreadedPreprocessor;
import hic.tools.utils.original.Preprocessor;
import javastraw.reader.type.NormalizationType;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class PreProcessing extends JuiceboxCLT {


    private String inputFile;
    private String outputFile;
    private Preprocessor preprocessor;
    private boolean noNorm = false;
    private boolean noFragNorm = false;
    private int genomeWide;
    private String shell = "sh";
    private final List<NormalizationType> normalizationTypes = new ArrayList<>();

    public PreProcessing() {
        super(getBasicUsage() + "\n"
                + "           : --intra only calculate intra chromosomal maps [false]\n"
                + "           : --near-diagonal only retain reads within 10MB of diagonal [false]\n"
                + "           : --block-size <int> set block capacity [1000]\n"
                + "           : -m <int> only write cells with count above threshold m [0]\n"
                + "           : -q <int> filter by MAPQ score greater than or equal to q [not set]\n"
                + "           : -c <chromosome ID> only calculate map on specific chromosome [not set]\n"
                + "           : -r <comma-separated list of resolutions> Only calculate specific resolutions [not set]\n"
                + "           : -t <tmpDir> Set a temporary directory for writing\n"
                + "           : -s <statistics file> Add the text statistics file to the Hi-C file header\n"
                + "           : -g <graphs file> Add the text graphs file to the Hi-C file header\n"
                + "           : -n Don't normalize the matrices\n"
                + "           : -z <double> scale factor for hic file\n"
                + "           : -a <1, 2, 3, 4, 5> filter based on inner, outer, left-left, right-right, tandem pairs respectively\n"
                + "           : --random_seed <long> for seeding random number generator\n"
                + "           : -k normalizations to include\n"
                + "           : -j number of CPU threads to use\n"
                + "           : --threads <int> number of threads \n"
                + "           : --mndindex <filepath> to mnd chr block indices\n"
                + "           : --conserve-ram will minimize RAM usage\n"
                + "           : --check-ram-usage will check ram requirements prior to running\n"
                + "           : --shell how to execute shell (sh, bash, zsh, etc); default: sh"
        );
    }

    public static String getBasicUsage() {
        return "pre [options] <infile> <outfile> <genomeID>";
    }

    @Override
    public void readArguments(String[] args, CommandLineParser parser) {

        String genomeId = "";
        try {
            genomeId = args[3];
        } catch (ArrayIndexOutOfBoundsException e) {
            System.err.println("No genome ID given");
            printUsageAndExit();
        }

        inputFile = args[1];
        outputFile = args[2];
        String tmpDir = parser.getTmpdirOption();
        double hicFileScalingFactor = parser.getScalingOption();

        HiCGlobals.primaryThreads = updateNumberOfCPUThreads(parser, 1);
        HiCGlobals.normThreads = updateSecondaryNumberOfCPUThreads(parser, 10);

        if (HiCGlobals.primaryThreads < 2) {
            preprocessor = new Preprocessor(new File(outputFile), genomeId, hicFileScalingFactor, tmpDir);
            usingMultiThreadedVersion = false;
        } else {
            try {
                preprocessor = new MultithreadedPreprocessor(new File(outputFile), genomeId,
                        hicFileScalingFactor, HiCGlobals.primaryThreads, parser.getMndIndexOption(), tmpDir);
                usingMultiThreadedVersion = true;
            } catch (Exception e) {
                System.err.println(e.getLocalizedMessage() + "\nUsing single threaded preprocessor");
                preprocessor = new Preprocessor(new File(outputFile), genomeId, hicFileScalingFactor, tmpDir);
                usingMultiThreadedVersion = false;
            }
        }

        preprocessor.setIncludedChromosomes(parser.getChromosomeSetOption());
        preprocessor.setCountThreshold(parser.getCountThresholdOption());
        preprocessor.setV9DepthBase(parser.getV9DepthBase());
        preprocessor.setMapqThreshold(parser.getMapqThresholdOption());
        preprocessor.setIntraChromosomalOnly(parser.getDiagonalsOption());
        preprocessor.setStatisticsFile(parser.getStatsOption());
        preprocessor.setGraphFile(parser.getGraphOption());
        preprocessor.setGenome(parser.getGenomeOption());
        preprocessor.setResolutions(parser.getResolutionOption());
        preprocessor.setFilter(parser.getAlignmentOption());
        int blockCapacity = parser.getBlockCapacityOption();
        if (blockCapacity > 10) {
            Preprocessor.BLOCK_CAPACITY = blockCapacity;
        }

        String customShell = parser.getShellOption();
        if (customShell != null && customShell.length() > 0) {
            shell = customShell;
        }
        noNorm = parser.getNoNormOption();
        genomeWide = parser.getGenomeWideOption();
        noFragNorm = parser.getNoFragNormOption();
        normalizationTypes.addAll(parser.getAllNormalizationTypesOption());
    }

    @Override
    public void run() {
        try {
            long currentTime = System.currentTimeMillis();
            if (usingMultiThreadedVersion) {
                preprocessor.preprocess(inputFile, null, null, null);
                ShellCommandRunner.runShellFile(shell, outputFile + MultithreadedPreprocessor.CAT_SCRIPT);
            } else {
                preprocessor.preprocess(inputFile, outputFile, outputFile, null);
            }

            if (HiCGlobals.printVerboseComments) {
                System.out.println("\nBinning contact matrices took: " + (System.currentTimeMillis() - currentTime) + " milliseconds");
            }

            if (!noNorm) {
                Map<NormalizationType, Integer> resolutionsToBuildTo = AddNorm.defaultHashMapForResToBuildTo(normalizationTypes);
                AddNorm.launch(outputFile, normalizationTypes, 0, resolutionsToBuildTo);
            } else {
                System.out.println("Done creating .hic file. Normalization not calculated due to -n flag.");
                System.out.println("To run normalization, run: java -jar juicer_tools.jar addNorm <hicfile>");
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(56);
        }
    }
}