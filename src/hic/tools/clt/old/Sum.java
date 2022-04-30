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

import hic.tools.clt.CommandLineParser;
import hic.tools.clt.JuiceboxCLT;
import hic.tools.utils.original.Preprocessor;
import hic.tools.utils.original.PreprocessorFromDatasets;
import javastraw.reader.Dataset;
import javastraw.tools.HiCFileTools;

import java.io.File;

public class Sum extends JuiceboxCLT {

    private String outputFile;
    private PreprocessorFromDatasets pfd;

    public Sum() {
        super(getBasicUsage() + "\n"
                + " : [--intra]             only calculate intra chromosomal maps       [default: false]\n"
                + " : [--near-diagonal]     only retain reads within 10MB of diagonal   [default: false]\n"
                + " : [-t <string>]         set a temporary directory for writing       [default: temp_folder]\n"
                + " : [-r <int>]            set the highest resolution to build to      [default: highest available]\n"
                + " : [--block-size <int>]  scale factor to increase block capacity by  [default: 1]");
    }

    public static String getBasicUsage() {
        return "sum [options] <outfile.hic> <infile1.hic> <infile2.hic> ... <infileN.hic";
    }

    @Override
    public void readArguments(String[] args, CommandLineParser parser) {

        outputFile = args[1];
        Dataset[] datasets = new Dataset[args.length - 2];
        for (int z = 2; z < args.length; z++) {
            datasets[z - 2] = HiCFileTools.extractDatasetForCLT(args[z], false, false, false);
        }

        String tmpDir = parser.getTmpdirOption();
        double hicFileScalingFactor = parser.getScalingOption();
        updateNumberOfCPUThreads(parser, 10);

        pfd = new PreprocessorFromDatasets(new File(outputFile), datasets, hicFileScalingFactor);
        pfd.setTmpdir(tmpDir);
        pfd.setHighestResolution(parser.getResolutionOption());
        boolean getOnlyNearDiagonal = parser.getOnlyNearDiagonalOption();
        pfd.setIntraChromosomalOnly(parser.getDiagonalsOption() || getOnlyNearDiagonal);
        pfd.setOnlyNearDiagonalsOnly(getOnlyNearDiagonal);
        int blockCapacity = parser.getBlockCapacityOption();
        if (blockCapacity > 1) {
            Preprocessor.BLOCK_CAPACITY *= blockCapacity;
        }
    }

    @Override
    public void run() {
        try {
            pfd.preprocess();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(56);
        }
    }
}