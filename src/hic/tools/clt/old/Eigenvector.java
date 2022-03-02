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

package hic.tools.clt.old;

import hic.HiCGlobals;
import hic.tools.clt.CommandLineParser;
import hic.tools.clt.JuiceboxCLT;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Class for calculating Eigenvector (separated out from Dump)
 * @author Neva Durand
 */
public class Eigenvector extends JuiceboxCLT {

    private HiCZoom.HiCUnit unit = null;
    private int binSize = 0;
    private Chromosome chromosome1;
    private PrintWriter pw;

    public Eigenvector() {
        super(getUsage() + "\n\t-p, --pearsons_all_resolutions: calculate eigenvector at all resolutions");
    }

    public static String getUsage(){
        return "eigenvector -p <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]";
    }

    @Override
    public void readArguments(String[] args, CommandLineParser parser) {
        if (args.length != 7 && args.length != 6) {
            printUsageAndExit();
        }

        setDatasetAndNorm(args[2], args[1], false);

        ChromosomeHandler chromosomeHandler = dataset.getChromosomeHandler();

        if (chromosomeHandler.doesNotContainChromosome(args[3])) {
            System.err.println("Unknown chromosome: " + args[3]);
            System.exit(18);
        }
        chromosome1 = chromosomeHandler.getChromosomeFromName(args[3]);

        try {
            unit = HiCZoom.valueOfUnit(args[4]);
        } catch (IllegalArgumentException error) {
            System.err.println("Unit must be in BP or FRAG.");
            System.exit(20);
        }

        String binSizeSt = args[5];

        try {
            binSize = Integer.parseInt(binSizeSt);
        } catch (NumberFormatException e) {
            System.err.println("Integer expected for bin size.  Found: " + binSizeSt + ".");
            System.exit(21);
        }

        if ((unit == HiCZoom.HiCUnit.BP && binSize < HiCGlobals.MAX_EIGENVECTOR_ZOOM) ||
                (unit == HiCZoom.HiCUnit.FRAG && binSize < HiCGlobals.MAX_EIGENVECTOR_ZOOM / 1000)) {
            System.out.println("WARNING: Eigenvector calculation at high resolution can take a long time");
        }


        if (args.length == 7) {
            try {
                pw = new PrintWriter(new FileOutputStream(args[6]));
            } catch (IOException error) {
                System.err.println("Cannot write to " + args[6]);
                error.printStackTrace();
                System.exit(22);
            }
        }
        else pw = new PrintWriter(System.out);
    }

    @Override
    public void run() {
        HiCZoom zoom = new HiCZoom(unit, binSize);

        MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chromosome1, chromosome1, zoom);
        if (zd == null) {
            System.err.println("No reads in " + chromosome1);
            System.err.println("Unknown resolution: " + zoom);
            System.err.println("This data set has the following bin sizes (in bp): ");
            for (int zoomIdx = 0; zoomIdx < dataset.getNumberZooms(HiCZoom.HiCUnit.BP); zoomIdx++) {
                System.err.print(dataset.getZoom(HiCZoom.HiCUnit.BP, zoomIdx).getBinSize() + " ");
            }
            System.err.println("\nand the following bin sizes (in frag): ");
            for (int zoomIdx = 0; zoomIdx < dataset.getNumberZooms(HiCZoom.HiCUnit.FRAG); zoomIdx++) {
                System.err.print(dataset.getZoom(HiCZoom.HiCUnit.FRAG, zoomIdx).getBinSize() + " ");
            }
            System.exit(13);
        }
        ExpectedValueFunction df = dataset.getExpectedValuesOrExit(zd.getZoom(), norm, chromosome1, true, true);
        double[] vector = zd.getEigenvector(df, 0);

        // mean center and print
        int count = 0;
        double total = 0;

        for (double element : vector) {
            if (!Double.isNaN(element)) {
                total += element;
                count++;
            }
        }

        double mean = total / count; // sum is now mean

        // print out vector
        for (double element : vector) {
            pw.println(element - mean);
        }
        pw.close();

    }
}