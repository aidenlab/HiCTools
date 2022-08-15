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

package hic.tools.utils.norm.RU;

import java.util.Arrays;

public class RUCode() {
    public static int balance(long m, double[] b, double[] report, double[] allIters, float tol, float pppp,
                              int maxiter, double del, double dp, int totIter, int threads, int k) {

        int p;
        int n0;
        double low;
        int[] bad;
        int lind;
        double[] row;
        double[] r0;
        double[] one;
        double[] row0;

        double[] current = new double[k];
        row = new double[k];
        row0 = new double[k];
        r0 = new double[k];
        bad = new int[k];

        one = new double[k];

        double[] space = new double[threads];

        for (int c = 0; c < threads; c++) {
            space[c] = new double[k];
        }

        double perc = pppp;

        for (int g = 0; g < k; g++){
            one[g] = 1.0;
        }

        for (int g = 0; g < k; g++){
            bad[g] = 0;
        }


        double[] nz = new double[k];

        for (int g = 0; g < k; g++){
            nz[g] = 0;
        }

        for (p = 0; p < m; p++) {
            nz[i[p]] += 1.0;
            if (i[p] != j[p]) nz[j[p]] += 1.0;
        }

        n0 = 0;
        for (int g = 0; g < k; g++) if (nz[g] > 0) r0[n0++] = nz[g];
        Arrays.sort(r0);

        lind = (int) (n0 * perc + 0.5);
        if (lind < 0) lind = 0;
        low = r0[lind];
//	find rows which are identically 0 and the perc proportion of rows with the lowest row sums and exclude them
        for (int g = 0; g < k; g++) if (nz[g] < low) bad[g] = 1;
        //	if bad[p] is 1 we remove row p
        utmvMul(i, j, x, m, one, k, row0, threads, space);
        for (int g = 0; g < k; g++){
            row[g] = row0[g];
        }

        for (int g = 0; g < k; g++){
            one[g] = 1.0 - bad[g];
        }
        for (int g = 0; g < k; g++) if (bad[g] == 1) row[g] = 1.0;
        for (int g = 0; g < k; g++){
            b[g] = one[g] / Math.sqrt(row[g]);
        }

        //	start iterations
//	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
        double ber = 10.0 * (1.0 + tol);
        double err = 10.0 * (1.0 + tol);
        int iter = 0;
        for (int g = 0; g < k; g++){
            current[g] = b[g];
        }
        int fail;   // checks whether the convergence rate seems to be good enough; 0 if yes, 1 if no
        int nerr = 0;

        double[] errors = new double[10000];
        int all_iters = 0;

        while ((ber > tol || err > 5.0 * tol) && iter < maxiter) {
            iter++;
            all_iters++;
            fail = 1;

            utmvMul(i, j, x, m, b, k, row, threads, space);
            for (int g = 0; g < k; g++){
                row[g] *= b[g];
            }
            for (int g = 0; g < k; g++) if (bad[g] == 1){
                row[g] = 1.0;
            }
            for (int g = 0; g < k; g++){
                b[g] *= one[g] / Math.sqrt(row[g]);
            }
            //	calculate the current error
            ber = 0;
            for (int g = 0; g < k; g++) {
                if (bad[g] == 1) continue;
                if (Math.abs((b[g] - current[g]) / (b[g] + current[g])) > ber)
                    ber = Math.abs((b[g] - current[g]) / (b[g] + current[g]));
            }
            report[all_iters - 1] = ber;
            allIters[all_iters - 1] = iter;

            //	since calculating the error in row sums requires matrix-vector multiplication we are are doing this every 10
            // 	iterations
            if (iter % 10 == 0) {
                utmvMul(i, j, x, m, b, k, row, threads, space);
                err = 0;
                for (int g = 0; g < k; g++) {
                    if (bad[g] == 1) continue;
                    if (err < Math.abs(row[g] * b[g] - one[g]))
                        err = Math.abs(row[g] * b[g] - one[g]);
                }
                errors[nerr++] = err;
            }
            for (int g = 0; g < k; g++){
                current[g] = b[g];
            }
            //	check whether convergence rate is satisfactory
//	if less than 5 iterations (so less than 5 errors) and less than 2 row sums errors, there is nothing to chek
            if ((ber < tol) && (nerr < 2 || (nerr >= 2 && errors[nerr - 1] < 0.5 * errors[nerr - 2]))) continue;
//	otherwise check
            if (iter > 5) {
                for (int g = 1; g <= 5; g++) if (report[all_iters - g] * (1.0 + del) < report[all_iters - g - 1]) fail = 0;
                if (nerr >= 2 && errors[nerr - 1] > 0.75 * errors[nerr - 2]) fail = 1;
                if (iter >= maxiter) fail = 1;

                //	the scaling vector does not seen to converge well enough or the row sums erros do not decrease enough
                // 	increase perc and perc1 and continue
                if (fail == 1) {
                    perc += dp;
                    nerr = 0;
                    lind = (int) (n0 * perc + 0.5);
                    low = r0[lind];
                    for (int g = 0; g < k; g++) {
                        if (nz[g] < low) {
                            bad[g] = 1;
                            one[g] = 0;
                        }
                    }

                    ber = 10.0 * (1.0 + tol);
                    err = 10.0 * (1.0 + tol);
                    //	if the current error is larger than 5 iteration ago start from scratch, otherwise continue from the current
                    // 	position
                    if (report[all_iters - 1] > report[all_iters - 6]) {
                        for (p = 0; p < k; p++){
                            one[p] = 1.0 - bad[p];
                        }
                        for (p = 0; p < k; p++){
                            row0[p] *= one[p];
                        }
                        for (p = 0; p < k; p++){
                            row[p] = row0[p];
                        }
                        for (p = 0; p < k; p++){
                            current[p] = b[p];
                        }
                    } else for (p = 0; p < k; p++){
                        b[p] *= (1.0 - bad[p]);
                    }
                    iter = 0;

                }
            }
            //	if perc or perc1 reached upper bound or the total number of iterationbs is too high, exit
            if (perc > 0.2) break;
            if (all_iters > totIter) break;
        }
        //	find the final error in row sums
        utmvMul(i, j, x, m, b, k, row, threads, space);
        err = 0;
        for (p = 0; p < k; p++) {
            if (bad[p] == 1) continue;
            if (err < Math.abs(row[p] * b[p] - one[p])) err = Math.abs(row[p] * b[p] - one[p]);
        }
        report[all_iters + 1] = ber;
        report[all_iters + 2] = err;
        pppp = (float) perc;
        totIter = all_iters;

        for (p = 0; p < k; p++) if (bad[p] == 1) b[p] = Float.NaN;

        return (iter);
    }
}

