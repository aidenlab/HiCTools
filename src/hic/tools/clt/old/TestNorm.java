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
import hic.tools.utils.bigarray.BigContactList;
import hic.tools.utils.norm.IntraNorms;
import hic.tools.utils.norm.NormalizationCalculations;
import javastraw.reader.Dataset;
import javastraw.reader.DatasetReaderV2;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.type.HiCZoom;

import java.io.*;
import java.nio.charset.StandardCharsets;


public class TestNorm extends JuiceboxCLT {

    private final String norm = "SCALE";
    private int resolution = 1000;
    private String file, name, output;

    public TestNorm() {
        super(getBasicUsage());
    }

    public static String getBasicUsage() {
        return "testNorm [--threads number] <input.hic> <chromosome> <resolution> <output>";
    }

    @Override
    public void readArguments(String[] args, CommandLineParser parser) {
        if (parser.getHelpOption() || args.length != 5) {
            printUsageAndExit();
        }

        HiCGlobals.normThreads = updateNumberOfCPUThreads(parser, 10);

        file = args[1];
        name = args[2];
        resolution = Integer.parseInt(args[3]);
        output = args[4];
    }

    @Override
    public void run() {
        HiCGlobals.allowDynamicBlockIndex = false;

        try {
            DatasetReaderV2 reader = new DatasetReaderV2(file, false, false);
            Dataset ds = reader.read();
            HiCGlobals.verifySupportedHiCFileVersion(reader.getVersion());

            Chromosome chrom = ds.getChromosomeHandler().getChromosomeFromName(name);
            HiCZoom zoom = new HiCZoom(resolution);
            String stem = "NORM_" + chrom.getIndex() + "_" + zoom.getBinSize();

            BigContactList ba = IntraNorms.getBigArrayFromAndClearCache(ds, chrom, zoom, 0);
            System.out.println("Contacts loaded to RAM");
            if (ba != null) {
                NormalizationCalculations nc = new NormalizationCalculations(ba, zoom.getBinSize());
                ListOfFloatArrays vc = nc.computeVC();
                long s0 = System.nanoTime();
                ListOfFloatArrays scale = nc.computeSCALE(vc, stem);
                long s1 = System.nanoTime();
                System.out.println("Overall time taken " + ((s1 - s0) * 1e-9) + " seconds");
                export(scale, output);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void export(ListOfFloatArrays scale, String output) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(new OutputStreamWriter(
                    new BufferedOutputStream(new FileOutputStream(output)), StandardCharsets.UTF_8));
            for (float[] row : scale.getValues()) {
                for (float val : row) {
                    out.println(val);
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } finally {
            if (out != null) {
                out.flush();
                out.close();
            }
        }
    }

    public static void export(ListOfDoubleArrays scale, String output) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(new OutputStreamWriter(
                    new BufferedOutputStream(new FileOutputStream(output)), StandardCharsets.UTF_8));
            for (double[] row : scale.getValues()) {
                for (double val : row) {
                    out.println(val);
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } finally {
            if (out != null) {
                out.flush();
                out.close();
            }
        }
    }
}