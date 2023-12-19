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
     |              |-------- ComputeCofPoly4.m              : Compute the coefficient for the 4-order polynomial function
     |              |-------- ComputeObjGradNorm4.m          : Compute the objective value for some function
     |              |-------- ComputeObjGradNorm4_sqrt.m     : Compute the objective value for some function
     |              |-------- ComputeObjNorm1.m              : Compute the objective value for some function
     |              |-------- SpectralNorm.m                 : Estimating the spectral norm of the matrix G'G using the Lanczos algorithm
     |              |-------- ComputeObjNorm4.m              : Compute the objective value for some function
     |              |-------- GetDataStr.m                   : Get the data strings
     |              |-------- GetDataStr2.m                  : Get the data strings
     |              |-------- generate_y_DCSparseOpt.m       : Generate the measurement vector for the sparse optimization problem
     |              |-------- getdata_ica.m                  : Get the data for the ICA problem
     |              |-------- quadfrac2.m                    : solve the quadratic fractional problem
     |              |-------- solve_4order.m                 : Solve the 4-order polynomial one-dimensional subproblem
     |              |-------- topksum.m                      : Compute the sum of top-k value for a vector
     |
     |
------------------- proximal
     |              |
     |              |-------- nonconvex_prox_diff_l1topk.m   : Compute the proximal operator for l1 and top-k norm function
     |              |-------- nonconvex_prox_diff_l2l4.m     : Compute the proximal operator for l2 and l4 norm norm function
     |              |-------- nonconvex_prox_frac_l1.m       : Compute the proximal operator for l1 norm function
     |              |-------- nonconvex_prox_frac_l1topk.m   : Compute the proximal operator for l1 and top-k norm function
     |              |-------- nonconvex_prox_frac_l2.m       : Compute the proximal operator for l2 norm function
     |              |-------- nonconvex_prox_l1topk2.m       : Compute the proximal operator for l1 top-k norm function
     |              |-------- nonconvex_prox_frac_topk.m     : Compute the proximal operator for top-k norm function
     |              |-------- nonconvex_prox_frac_linf.m     : Compute the proximal operator for linf norm function
     |              |-------- prox_l1.m                      : Compute the proximal operator for l1 norm function
     |              |-------- prox_l2.m                      : Compute the proximal operator for l2 norm function
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


	 
