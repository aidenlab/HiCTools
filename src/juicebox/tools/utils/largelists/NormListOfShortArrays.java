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

package juicebox.tools.utils.largelists;

import javastraw.tools.ParallelizationTools;
import juicebox.HiCGlobals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * can't use <T> because we need to instantiate the array, otherwise that would have been nice
 */
public class NormListOfShortArrays {

	public static final int DEFAULT_LENGTH = 10000000;
	final long overallLength;
	final List<short[]> internalList = new ArrayList<>();

	public NormListOfShortArrays(long length) {
		this.overallLength = length;
		long tempLength = length;
		while (tempLength > 0) {
			if (tempLength < DEFAULT_LENGTH) {
				internalList.add(new short[(int) tempLength]);
				break;
			} else {
				internalList.add(new short[DEFAULT_LENGTH]);
				tempLength -= DEFAULT_LENGTH;
			}
		}
	}

	public NormListOfShortArrays(long totSize, short defaultValue) {
		this(totSize);
		for (short[] array : internalList) {
			Arrays.fill(array, defaultValue);
		}
	}

	public void clear() {
		internalList.clear();
	}

	public int get(long index) {
		if (index < overallLength) {
			int pseudoRow = (int) (index / DEFAULT_LENGTH);
			int pseudoCol = (int) (index % DEFAULT_LENGTH);
			return internalList.get(pseudoRow)[pseudoCol];
		} else {
			System.err.println("long index exceeds max size of list of int arrays while getting");
			return -Integer.MAX_VALUE;
		}
	}

	public void set(long index, short value) {
		if (index < overallLength) {
			int pseudoRow = (int) (index / DEFAULT_LENGTH);
			int pseudoCol = (int) (index % DEFAULT_LENGTH);
			internalList.get(pseudoRow)[pseudoCol] = value;
		} else {
			System.err.println("long index exceeds max size of list of arrays while setting");
			return;
		}
		//System.err.println("unusual - long index exceeds max size of list of arrays while setting");
		return;
	}

	public long getLength() {
		return overallLength;
	}

	public NormListOfFloatArrays deepConvertedClone() {
		NormListOfFloatArrays clone = new NormListOfFloatArrays(overallLength);
		for (int k = 0; k < internalList.size(); k++) {
			float[] dest = clone.internalList.get(k);
			short[] src = internalList.get(k);
			for (int q = 0; q < dest.length; q++) {
				dest[q] = src[q];
			}
		}
		return clone;
	}

	public List<short[]> getValues() {
		return internalList;
	}

	public void parSetTo(NormListOfFloatArrays srcArrays) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(HiCGlobals.numCPUMatrixThreads, () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				short[] dest = internalList.get(i);
				float[] src = srcArrays.internalList.get(i);
				for (int z = 0; z < dest.length; z++) {
					dest[z] = (short) src[z];
				}
				i = index.getAndIncrement();
			}
		});
	}
}
