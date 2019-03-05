/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2019 Broad Institute, Aiden Lab
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
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package juicebox.tools.clt.juicer;

import juicebox.data.ChromosomeHandler;
import juicebox.data.Dataset;
import juicebox.data.HiCFileTools;
import juicebox.tools.clt.CommandLineParserForJuicer;
import juicebox.tools.clt.JuicerCLT;
import juicebox.tools.utils.juicer.hiccups.HiCCUPSConfiguration;
import juicebox.tools.utils.juicer.hiccups.HiCCUPSUtils;
import juicebox.track.feature.Feature2DList;
import juicebox.track.feature.Feature2DParser;
import juicebox.track.feature.Feature2DTools;
import juicebox.windowui.HiCZoom;
import juicebox.windowui.NormalizationType;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * HiCCUPS Diff
 * <p/>
 * Developed by Suhas Rao, ported by Suhas Rao + Neva Durand
 * <p/>
 * -------
 * Takes as input two Hi-C maps and their associated loop calls <br>
 * Outputs two differential loop lists: loops that are in the first but not in the second
 * and loops that are in the second but not in the first<br>
 * Other parameters are used for the two HiCCUPS calls on the alternate loop lists
 * @see juicebox.tools.clt.juicer.HiCCUPS
 * hiccupsdiff [-m matrixSize] [-k normalization (NONE/VC/VC_SQRT/KR)] [-c chromosome(s)] [-r resolution(s)]
 * [-f fdr] [-p peak width] [-i window] [-t thresholds] [-d centroid distances]
 *  <firstHicFile> <secondHicFile> <firstLoopList> <secondLoopList> <outputDirectory>
 *  firstLoopList is the loop list generated by running HiCCUPS on firstHicFile
 *  secondLoopList is the loop list generated by running HiCCUPS on secondHicFile
 */
public class HiCCUPSDiff extends JuicerCLT {

    private final float maxEnrich = 1.3f;
    private HiCCUPS hiccups1 = null;
    private HiCCUPS hiccups2 = null;
    private Feature2DList looplist1;
    private Feature2DList looplist2;
    private File outputDirectory;
    private ChromosomeHandler commonChromosomesHandler;
    private List<HiCCUPSConfiguration> configs;

    public HiCCUPSDiff() {
        // what variables should they be able to send in?
        // need to add maxEnrich
        super("hiccupsdiff [-m matrixSize] [-k normalization (NONE/VC/VC_SQRT/KR)] [-c chromosome(s)] " +
                "[-f fdr] [-p peak width] [-i window] [-t thresholds] [-d centroid distances] " +
                "<firstHicFile> <secondHicFile> <firstLoopList> <secondLoopList> <outputDirectory>");
    }

    public static String getBasicUsage() {
        return "hiccupsdiff <firstHicFile> <secondHicFile> <firstLoopList> <secondLoopList> <outputDirectory>";
    }

    @Override
    protected void readJuicerArguments(String[] args, CommandLineParserForJuicer juicerParser) {
        if (args.length != 6) {
            printUsageAndExit();
        }

        outputDirectory = HiCFileTools.createValidDirectory(args[5]);

        Dataset ds1 = HiCFileTools.extractDatasetForCLT(Arrays.asList(args[1].split("\\+")), true);
        Dataset ds2 = HiCFileTools.extractDatasetForCLT(Arrays.asList(args[2].split("\\+")), true);

        if (!(ds1.getGenomeId().equals(ds2.getGenomeId()))) {
            System.err.println("Hi-C maps must be from the same genome");
            System.exit(27);
        }
        // intersecting for the edge case where one of the hic files may not be using all chromosomes
        // e.g. the mbr_19 files for testing
        commonChromosomesHandler = HiCFileTools.getChromosomeSetIntersection(ds1.getChromosomeHandler(),
                ds2.getChromosomeHandler());

        if (givenChromosomes != null && givenChromosomes.size() > 0)
            commonChromosomesHandler = HiCFileTools.stringToChromosomes(givenChromosomes, commonChromosomesHandler);

        List<HiCZoom> availableZooms = new ArrayList<>(HiCFileTools.getZoomSetIntersection(
                ds1.getBpZooms(), ds1.getBpZooms()));

        looplist1 = Feature2DParser.loadFeatures(args[3], commonChromosomesHandler, true, null, false);
        looplist2 = Feature2DParser.loadFeatures(args[4], commonChromosomesHandler, true, null, false);

        configs = HiCCUPSConfiguration.extractConfigurationsFromCommandLine(juicerParser, availableZooms);

        if (configs == null) {
            configs = new ArrayList<>();
            if (Feature2DTools.isResolutionPresent(looplist1, 5000) && Feature2DTools.isResolutionPresent(looplist2, 5000)) {
                configs.add(HiCCUPSConfiguration.getDefaultConfigFor5K());
            }
            if (Feature2DTools.isResolutionPresent(looplist1, 10000) && Feature2DTools.isResolutionPresent(looplist2, 10000)) {
                configs.add(HiCCUPSConfiguration.getDefaultConfigFor10K());
            }
            if (Feature2DTools.isResolutionPresent(looplist1, 25000) && Feature2DTools.isResolutionPresent(looplist2, 25000)) {
                configs.add(HiCCUPSConfiguration.getDefaultConfigFor25K());
            }
            if (configs.size() == 0) {
                System.err.println("The loop lists have no resolutions in common.");
                System.exit(28);
            }
        }

        System.out.println("Running differential HiCCUPs with resolutions:");
        for (HiCCUPSConfiguration config : configs) {
            System.out.println(config);
        }

        boolean processed = true;
        for (HiCCUPSConfiguration config : configs) {
            String fname = outputDirectory + File.separator + "file1" + File.separator + "requested_list_" + config.getResolution() + ".bedpe";
            if (!new File(fname).exists()) processed = false;
            fname = outputDirectory + File.separator + "file2" + File.separator + "requested_list_" + config.getResolution() + ".bedpe";
            if (!new File(fname).exists()) processed = false;
        }

        if (processed) {
            System.out.println("Using already created differential lists in " + outputDirectory + File.separator +
                    "file1 and " + outputDirectory + File.separator + "file2");
        }
        else {

            NormalizationType preferredNorm = juicerParser.getNormalizationTypeOption(ds1.getNormalizationHandler());
            if (preferredNorm != null)
                norm = preferredNorm;

            int matrixSize = juicerParser.getMatrixSizeOption();
            if (matrixSize <= 0) matrixSize = 1024;

            boolean usingCPUVersion = false;
            if (juicerParser.getCPUVersionOfHiCCUPSOptions()) {
                usingCPUVersion = true;
                System.out.println(HiCCUPS.CPU_VERSION_WARNING);
            }

            double[] thresholds = null;
            List<String> t = juicerParser.getThresholdOptions();
            if (t != null && t.size() == 4) {
                thresholds = HiCCUPSUtils.extractDoubleValues(t, 4, Double.NaN);
            }

            int numThreads = juicerParser.getNumThreads();

            System.out.println("Running HiCCUPS with alternate loop lists");
            hiccups1 = new HiCCUPS();
            hiccups2 = new HiCCUPS();
            hiccups1.initializeDirectly(ds1, outputDirectory + File.separator + "file1", args[4],
                    norm, matrixSize, commonChromosomesHandler, configs, thresholds, usingCPUVersion, numThreads);
            hiccups2.initializeDirectly(ds2, outputDirectory + File.separator + "file2", args[3],
                    norm, matrixSize, commonChromosomesHandler, configs, thresholds, usingCPUVersion, numThreads);
        }
    }

    @Override
    public void run() {

        if (hiccups1 != null && hiccups2 != null) {
            hiccups1.run();
            hiccups2.run();
        }

        // for every feature in second loop list, see if there's a reasonably close one in first list
        Feature2DList conservedLoopList2 = Feature2DTools.extractReproducibleCentroids(looplist1, looplist2, 50000, 0.2);
        // for every feature in first loop list, see if there's a reasonably close one in second list
        Feature2DList conservedLoopList1 = Feature2DTools.extractReproducibleCentroids(looplist2, looplist1, 50000, 0.2);

        // get the differences - loops that appear only in looplist1
        Feature2DList diff1 = Feature2DTools.compareLists(conservedLoopList1, looplist1, false);
        // get the differences - loops that appear only in looplist2
        Feature2DList diff2 = Feature2DTools.compareLists(conservedLoopList2, looplist2, false);

        // load all the loops resulting from running HiCCUPs on the first HiC file with the second loop list
        // then filter by max enrichment: observed < maxEnrich*expected BL & donut & V & H
        Feature2DList results1 = new Feature2DList();
        for (HiCCUPSConfiguration config : configs) {
            String fname = outputDirectory + File.separator + "file1" + File.separator + "requested_list_" + config.getResolution() + ".bedpe";
            Feature2DList requestedList = Feature2DParser.loadFeatures(fname, commonChromosomesHandler, true, null, false);
            HiCCUPSUtils.filterOutFeaturesByEnrichment(requestedList, maxEnrich);
            results1.add(requestedList);
        }

        // load all the loops resulting from running HiCCUPs on the second HiC file with the first loop list
        // then filter by max enrichment: observed < maxEnrich*expected BL & donut & V & H
        Feature2DList results2 = new Feature2DList();
        for (HiCCUPSConfiguration config : configs) {
            String fname = outputDirectory + File.separator + "file2" + File.separator + "requested_list_" + config.getResolution() + ".bedpe";
            Feature2DList requestedList = Feature2DParser.loadFeatures(fname, commonChromosomesHandler, true, null, false);
            HiCCUPSUtils.filterOutFeaturesByEnrichment(requestedList, maxEnrich);
            results2.add(requestedList);
        }
        // differential loop list 1 is loops that appeared in list1 that are not enriched in Hi-C file 2
        Feature2DList differentialList1 = Feature2DList.getIntersection(diff1, results2);
        // differential loop list 2 is loops that appeared in list2 that are not enriched in Hi-C file 1
        Feature2DList differentialList2 = Feature2DList.getIntersection(diff2, results1);

        differentialList1.exportFeatureList(new File(outputDirectory, "differential_loops1.bedpe"), true, Feature2DList.ListFormat.FINAL);
        differentialList2.exportFeatureList(new File(outputDirectory, "differential_loops2.bedpe"), true, Feature2DList.ListFormat.FINAL);
    }
}
