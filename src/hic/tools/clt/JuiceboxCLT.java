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

package hic.tools.clt;

import javastraw.reader.Dataset;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

/**
 * All command line tools should extend from this class
 */
public abstract class JuiceboxCLT {

    private static String usage;
    protected Dataset dataset = null;
    protected NormalizationType norm = null;
    protected boolean usingMultiThreadedVersion = false;

    protected JuiceboxCLT(String usage) {
        setUsage(usage);
    }

    public static String[] splitToList(String nextLine) {
        return nextLine.trim().split("\\s+");
    }

    public abstract void readArguments(String[] args, CommandLineParser parser);

    public abstract void run();

    private void setUsage(String newUsage) {
        usage = newUsage;
    }

    public void printUsageAndExit() {
        System.out.println("Usage:   juicer_tools " + usage);
        System.exit(0);
    }

    public void printUsageAndExit(int exitcode) {
        System.out.println("Usage:   juicer_tools " + usage);
        System.exit(exitcode);
    }

    protected void setDatasetAndNorm(String file, String normType, boolean allowPrinting) {
        dataset = HiCFileTools.extractDatasetForCLT(file, allowPrinting, false, false);
        norm = dataset.getNormalizationHandler().getNormTypeFromString(normType);
        if (norm == null) {
            System.err.println("Normalization type " + norm + " unrecognized.  Normalization type must be one of \n" +
                    "\"NONE\", \"VC\", \"VC_SQRT\", \"KR\", \"GW_KR\"," +
                    " \"GW_VC\", \"INTER_KR\", \"INTER_VC\", or a custom added normalization.");
            System.exit(16);
        }
    }

    public static int getAppropriateNumberOfThreads(int numThreads, int defaultNum) {
        if (numThreads > 0) {
            return numThreads;
        } else if (numThreads < 0) {
            return Math.abs(numThreads) * Runtime.getRuntime().availableProcessors();
        } else {
            return defaultNum;
        }
    }

    protected int updateNumberOfCPUThreads(CommandLineParser parser, int numDefaultThreads) {
        int numCPUThreads = getAppropriateNumberOfThreads(parser.getNumThreads(), numDefaultThreads);
        System.out.println("Using " + numCPUThreads + " CPU thread(s) for primary task");
        return numCPUThreads;
    }

    protected int updateSecondaryNumberOfCPUThreads(CommandLineParser parser, int numDefaultThreads) {
        int numCPUThreadsForSecondTask = getAppropriateNumberOfThreads(parser.getNumMatrixOperationThreads(), numDefaultThreads);
        System.out.println("Using " + numCPUThreadsForSecondTask + " CPU thread(s) for secondary task");
        return numCPUThreadsForSecondTask;
    }
}

