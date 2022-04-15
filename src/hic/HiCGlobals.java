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

package hic;

import javastraw.reader.mzd.MatrixZoomData;

/**
 * @author Muhammad Shamim
 * @since 11/25/14
 */
public class HiCGlobals {
    public static final String versionNum = "3.14.01";
    public static final int writingVersion = 9;
    public static final int bufferSize = 2097152;
    public static int MAX_PEARSON_ZOOM = 50000;
    public static final int MAX_EIGENVECTOR_ZOOM = 250000;
    public static boolean useCache = true;
    public static boolean allowDynamicBlockIndex = true;
    public static boolean printVerboseComments = false;
    public static boolean USE_ITERATOR_NOT_ALL_IN_RAM = false;
    public static boolean CHECK_RAM_USAGE = false;
    public static int numCPUMatrixThreads = 1;

    public static void verifySupportedHiCFileVersion(int version) throws RuntimeException {
        if (version != writingVersion) {
            System.err.println("This file is version " + version +
                    ". Only version " + writingVersion + " files are supported with this jar.");
            System.exit(3);
        }
    }

    public static void verifySupportedHiCFileWritingVersion(int version) throws RuntimeException {
        if (version != writingVersion) {
            System.err.println("This file is version " + version +
                    ". Only version " + writingVersion + " files can be edited using this jar.");
            System.exit(2);
        }
    }

    public static void setMatrixZoomDataRAMUsage() {
        MatrixZoomData.useIteratorDontPutAllInRAM = USE_ITERATOR_NOT_ALL_IN_RAM;
        MatrixZoomData.shouldCheckRAMUsage = CHECK_RAM_USAGE;
    }
}
