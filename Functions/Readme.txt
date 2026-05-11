Modification by Daniel Kaufmann, 11/05/2026: Files have been modified to run in Octave 8.x. Find the original files on https://karelmertens.com/research/


Code for: `A Robust Test for Weak Instruments with Multiple Endogenous Regressors', by Daniel Lewis and Karel Mertens 
This version 22/09/2024

gweakivtest_Example.m  : Example code using simulated data, start with this file. Includes implementations in both local-to-zero and local-to-rank-reduction-of-1 settings and full-vector and single-coefficient tests
gweakivtest.m: Outputs generalized test statistic and critical values in local-to-zero (instruments are uniformly weak) setting. Also outputs Stock-Yogo statistic and critical value (based on Nagar approximation).
gweakivtest_critical_values.m: Calculates critical values based on both the worst case Nagar bias and the simplified bound in local-to-zero setting. 
gweakivtestLRR1.m: Outputs generalized test statistic and critical values in local-to-rank-reduction-of-1 (first stage is close to rank deficient) setting. Also outputs Stock-Yogo statistic and critical value (based on Nagar approximation).
gweakivtest_critical_valuesLRR1.m: Calculates critical values based on both the worst case Nagar bias and the simplified bound in local-to-rank-reduction-of-1 setting. 

Note: gweakivtest_critical_values.m and gweakivtest_critical_valuesLRR1.m make use of the curvilinear search algorithm for optimization on Stiefel manifold by Z. Wen and W. Yin (A feasible method for optimization with orthogonality constraints, In: Mathematical Programming 142, pp. 397–434. doi: 10.1007/s10107-012-0584-1.20) 
 





