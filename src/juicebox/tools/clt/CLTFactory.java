/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
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

package juicebox.tools.clt;

import juicebox.HiCGlobals;
import juicebox.tools.clt.old.*;


/**
 * Factory for command line tools to call different functions
 *
 * @author Muhammad Shamim
 * @since 1/30/2015
 */
public class CLTFactory {

    // Commenting some out because we're not going to release all these when we release CLT
    private final static String[] commandLineToolUsages = {
            PreProcessing.getBasicUsage(),
            AddNorm.getBasicUsage(),
            Pearsons.getBasicUsage(),
            Eigenvector.getUsage(),
    };

    public static void generalUsage() {
        System.out.println("Juicer Tools Version " + HiCGlobals.versionNum);
        System.out.println("Usage:");
        for (String usage : commandLineToolUsages) {
            System.out.println("\t" + usage);
        }
        System.out.println("\t" + "-h, --help print help");
        System.out.println("\t" + "-v, --verbose verbose mode");
        System.out.println("\t" + "-V, --version print version");
        System.out.println("Type juicer_tools <commandName> for more detailed usage instructions");
    }

    public static JuiceboxCLT getCLTCommand(String cmd) {

        cmd = cmd.toLowerCase();
        if (cmd.equals("pre")) {
            return new PreProcessing();
        } else if (cmd.equals("addnorm")) {
            return new AddNorm();
        } else if (cmd.equals("bigwig")) {
            return new BigWig();
        } else if (cmd.equals("bintopairs")) {
            return new BinToPairs();
        } else if (cmd.equals("bptofrag")) {
            return new BPToFragment();
        } else if (cmd.equals("fragmenttobed")) {
            return new FragmentToBed();
        } else if (cmd.equals("pairstobin")) {
            return new PairsToBin();
        } else if (cmd.equals("pearsons")) {
            return new Pearsons();
        } else if (cmd.equals("eigenvector")) {
            return new Eigenvector();
        } else if (cmd.contains("librarycomplexity")) {
            return new LibraryComplexity(cmd);
        } else if (cmd.equals("statistics")) {
            return new Statistics();
        }

        return null;
    }
}
