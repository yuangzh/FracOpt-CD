CDDC README File
=================================================
1. Running the Code
To reproduce the experiments in [1], you need to perform three steps: 
(a) Run mexC.m file. Note that this step is optional and only used for the One-Hidden-Layer Networks Problem.
(b) To obtain the accuracy for the 6 applications, run demo_accuracy1_*.m / demo_accuracy2_*.m / demo_accuracy3_*.m / demo_accuracy4_*.m / demo_accuracy5_*.m / demo_accuracy6_*.m 
(c) To obtain the efficiency for the 6 applications, run demo_cpu1_*.m / demo_cpu2_*.m / demo_cpu3_*.m / demo_cpu4_*.m / demo_cpu5_*.m / demo_cpu6_*.m 

2. Directory tree
We demonstrate the directory tree of our code, including the names of some main functions and their descriptions.

------------------- solvers
     |              |
     |              |-------- ICA_FractionalCD.m            : Implementation of Toland-Dual Method for the L1 Norm Generalized Eigenvalue Problem
     |              |-------- ICA_PGSA.m                    : Implementation of SubGradient Method for the L1 Norm Generalized Eigenvalue Problem
     |              |-------- ICA_PowerMethod.m             : Implementation of CD-SCA for the L1 Norm Generalized Eigenvalue Problem
     |              |
     |              |
     |              |-------- SpaseRecovery_FractionalCD.m           :  The objective function used in the L1 Norm Generalized Eigenvalue Problem
     |              |-------- SpaseRecovery_DPA.m           :  The objective function used in the L1 Norm Generalized Eigenvalue Problem
     |              |-------- SpaseRecovery_DPA2.m          : The objective function used in the L1 Norm Generalized Eigenvalue Problem
     |              |-------- SpaseRecovery_PGSA.m          : The objective function used in the L1 Norm Generalized Eigenvalue Problem
     |              |-------- SpaseRecovery_PGSA2.m         : The objective function used in the L1 Norm Generalized Eigenvalue Problem
     |              |-------- SpaseRecovery_ParametricCD.m  : The objective function used in the L1 Norm Generalized Eigenvalue Problem
     |              |-------- SpaseRecovery_QTA.m           : The objective function used in the L1 Norm Generalized Eigenvalue Problem
     |              |
     |              |
     |              |-------- AccerlatedProximalGradient.m  : Implementation of Multi-Stage Convex Relaxation for the L1 Norm Generalized Eigenvalue Problem
     |              |-------- ProximalGradient.m            : Implementation of CD-SNCA for the L1 Norm Generalized Eigenvalue Problem
     |
     |
     |
------------------- util
     |              |
     |              |-------- ComputeCofPoly4.m                          : Accelerated Proximal Gradient Method for solving convex composite optimization problems
     |              |-------- ComputeObjGradNorm4.m        : Compute the auxiliary matrix for the Phase Retrieval Problem
     |              |-------- ComputeObjGradNorm4_sqrt.m                   : Euclidean projection for the binary constraints: min_v ||v - a||_2^2, s.t. v \in {-1,+1}^n
     |              |-------- ComputeObjNorm1.m                      : Proximal operator for the L1 norm: min_{x} 0.5||x-a||^2 + lambda*||x||_1
     |              |-------- SpectralNorm.m                 : Estimating the spectral norm of the matrix G'G using the Lanczos algorithm
     |              |-------- ComputeObjNorm4.m                : Compute the subgradient of the top k norm function (in absolute value)
     |              |-------- GetDataStr.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- GetDataStr2.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- generate_y_DCSparseOpt.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- getdata_ica.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- prox_l1.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- prox_l2.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- quadfrac2.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- solve_4order.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- topksum.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |
------------------- proximal
     |              |
     |              |-------- nonconvex_prox_diff_l1topk.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_diff_l2l4.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_frac_l1.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_frac_l1topk.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_frac_l2.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_frac_l2.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_l1topk2.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_frac_topk.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |              |-------- nonconvex_prox_frac_linf.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
     |
     |
     |
------------------- data
     |              |
     |              |-------- *.mat                          : fourth matlab datasets
     |
     |
     |
------------------- DrawFig
                    |
                    |-------- demo_cpu_ICA.m                 : Generate the figures for the ICA problem
                    |-------- demo_cpu_Sparse.m              : Generate the figures for the sparse optimization problem
	 
	 
3. REFERENCES:
[1] Ganzhao Yuan. Ganzhao Yuan. Coordinate Descent Methods for Fractional Minimization. International Conference on Machine Learning (ICML), 2023.


	 
