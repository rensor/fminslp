# fminslp
Matlab based optimizer framework for Sequential Linear Programming (SLP) coupled with a global convergence filter and an adaptive move-limit strategy. The algorithm can handle linear equality, linear in-equality, non-linear in-equality, and non-linear equality constraints. A merit function approach is applied to ensure unconditional feasibility of the linearized sub problem.

The global convergence filter monitors progression of the penalized objective function (the merit function), and the associated non-linear constraints. Based on the response, the filter algorithm adjusts the move-limits to ensure stable convergence. The convergence filter is based on the following papers by Chin CM, Fletcher R (1999) and Fletcher R, Leyffer S, Toint PL (1998).

The adaptive move-limit strategy controls the box-constraints (upper and lower bounds for the design variables) and is based on the work by professor Erik Lund from Aalborg University.

The overall framework is based on an implementation developed during my Ph.d. studies. It was first used for the following paper: 

R Soerensen, E Lund (2015): Thickness filters for gradient based multi-material and thickness optimization of laminated composite structures, Structural and Multidisciplinary Optimization 52 (2), 227-250
https://doi.org/10.1007/s00158-015-1230-3.
 
Please refer to this paper when citing the fminslp algorithm.

References: 
Chin CM, Fletcher R (1999): On the global convergence of an SLP-
filter algorithm that takes EQP steps. Numerical Analysis Report
NA/199, Department of Mathematics, University of Dundee,Scotland, UK

Fletcher R, Leyffer S, Toint PL (1998): On the global convergence
of an SLP-filter algorithm. Numerical Analysis Report NA/183,
Department of Mathematics, University of Dundee, Scotland, UK

Usage: see also examples.m 

myProblem = fminslp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon); % Default options are applied

myProblem = fminslp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options); % Input options structure

myProblem = fminslp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,'MoveLimit',0.01,'MaxIterations',100); % Input/modify individual options

[x,fval,exitflag,output] = myProblem.solve(); % Solve problem

Available options:

% Initialize values to default
options = fminslp.slpoptions();

                    Algorithm: 'merit'
                      Display: 'off'
    InfeasibilityPenalization: 1000
       MaxFunctionEvaluations: 1000
             MaxInfeasibility: Inf
                MaxIterations: 1000
                    MoveLimit: 0.1000
              MoveLimitExpand: 1.1000
              MoveLimitMethod: 'adaptive'
              MoveLimitReduce: 0.5000
          OptimalityTolerance: 1.0000e-06
                       Solver: 'linprog'
                StepTolerance: 1.0000e-10
                
% Modify individual options, the rest are initialized to default values
options = fminslp.slpoptions('MoveLimit',0.01); 
                    Algorithm: 'merit'
                      Display: 'off'
    InfeasibilityPenalization: 1000
       MaxFunctionEvaluations: 1000
             MaxInfeasibility: Inf
                MaxIterations: 1000
                    MoveLimit: 0.0100
              MoveLimitExpand: 1.1000
              MoveLimitMethod: 'adaptive'
              MoveLimitReduce: 0.5000
          OptimalityTolerance: 1.0000e-06
                       Solver: 'linprog'
                StepTolerance: 1.0000e-10
