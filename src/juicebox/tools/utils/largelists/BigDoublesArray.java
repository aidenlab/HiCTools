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

import java.util.ArrayList;
import java.util.List;

/**
 * can't use <T> because we need to instantiate the array, otherwise that would have been nice
 */
public class BigDoublesArray {

	final long DEFAULT_LENGTH = BigShortsArray.DEFAULT_LENGTH;
	final long overallLength;
	final List<double[]> internalList = new ArrayList<>();

	public BigDoublesArray(long length) {
		this.overallLength = length;
		long tempLength = length;
		while (tempLength > 0) {
			if (tempLength < DEFAULT_LENGTH) {
				internalList.add(new double[(int) tempLength]);
				break;
			} else {
				internalList.add(new double[(int) DEFAULT_LENGTH]);
				tempLength -= DEFAULT_LENGTH;
			}
		}
	}

	public void clear() {
		internalList.clear();
	}

	public double get(long index) {
		if (index < overallLength) {
			int pseudoRow = (int) (index / DEFAULT_LENGTH);
			int pseudoCol = (int) (index % DEFAULT_LENGTH);
			return internalList.get(pseudoRow)[pseudoCol];
		} else {
			System.err.println("long index exceeds max size of list of arrays while getting: " + index + " " + overallLength);
			Exception ioe = new Exception();
			ioe.printStackTrace();
			return Double.NaN;
		}
	}

	public void set(long index, double value) {
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

	public void addTo(long index, double value) {
		if (index < overallLength) {
			int pseudoRow = (int) (index / DEFAULT_LENGTH);
			int pseudoCol = (int) (index % DEFAULT_LENGTH);
			try {
				internalList.get(pseudoRow)[pseudoCol] += value;
			} catch (Exception e) {
				System.err.println(index + " " + pseudoCol);
				e.printStackTrace();
			}
		} else {
			System.err.println("long index exceeds max size of list of arrays while adding: " + index + " " + overallLength);
			Exception ioe = new Exception();
			ioe.printStackTrace();
		}
	}

	public void addValuesFrom(BigDoublesArray other) {
		if (overallLength == other.overallLength) {
			for (int i = 0; i < internalList.size(); i++) {
				for (int j = 0; j < internalList.get(i).length; j++) {
					internalList.get(i)[j] += other.internalList.get(i)[j];
				}
			}
		} else {
			System.err.println("Adding objects of different sizes!");
		}
	}

	public List<double[]> getValues() {
		return internalList;
	}

	public BigFloatsArray convertToFloats() {
		BigFloatsArray newList = new BigFloatsArray(overallLength);
		for (int j = 0; j < internalList.size(); j++) {
			float[] dest = newList.getValues().get(j);
			double[] src = internalList.get(j);
			for (int k = 0; k < dest.length; k++) {
				dest[k] = (float) src[k];
			}
		}
		return newList;
	}
}
