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

package hic.tools.utils.largelists;

import hic.HiCGlobals;
import javastraw.tools.ParallelizationTools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * can't use <T> because we need to instantiate the array, otherwise that would have been nice
 */
public class BigIntsArray {

	public static final int DEFAULT_LENGTH = 10000;
	final long overallLength;
	final List<int[]> internalList = new ArrayList<>();

	public BigIntsArray(long length) {
		this.overallLength = length;
		long tempLength = length;
		while (tempLength > 0) {
			if (tempLength < DEFAULT_LENGTH) {
				internalList.add(new int[(int) tempLength]);
				break;
			} else {
				internalList.add(new int[DEFAULT_LENGTH]);
				tempLength -= DEFAULT_LENGTH;
			}
		}
	}

	public BigIntsArray(long totSize, short defaultValue) {
		this(totSize);
		for (int[] array : internalList) {
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
		}
	}

	public long getLength() {
		return overallLength;
	}

	public BigFloatsArray deepConvertedClone() {
		BigFloatsArray clone = new BigFloatsArray(overallLength);
		for (int k = 0; k < internalList.size(); k++) {
			float[] dest = clone.internalList.get(k);
			int[] src = internalList.get(k);
			for (int q = 0; q < dest.length; q++) {
				dest[q] = src[q];
			}
		}
		return clone;
	}

	public List<int[]> getValues() {
		return internalList;
	}

	public void parSetTo(BigFloatsArray srcArrays) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				int[] dest = internalList.get(i);
				float[] src = srcArrays.internalList.get(i);
				for (int z = 0; z < dest.length; z++) {
					dest[z] = (short) src[z];
				}
				i = index.getAndIncrement();
			}
		});
	}

	public BigIntsArray deepClone() {
		BigIntsArray clone = new BigIntsArray(overallLength);
		for (int k = 0; k < internalList.size(); k++) {
			System.arraycopy(internalList.get(k), 0, clone.internalList.get(k), 0, internalList.get(k).length);
		}
		return clone;
	}

	private int getNumThreads() {
		return Math.min(HiCGlobals.normThreads, internalList.size());
	}

	public void setAll(int value) {
		for (int[] row : internalList) {
			Arrays.fill(row, value);
		}
	}

	public void setToOneMinus(BigIntsArray input) {
		for (int i = 0; i < internalList.size(); i++) {
			int[] dest = internalList.get(i);
			int[] in = input.internalList.get(i);

			for (int j = 0; j < dest.length; j++) {
				dest[j] = 1 - in[j];
			}
		}
	}
}
