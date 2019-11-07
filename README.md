# fminslp
Matlab based optimizer framework for Sequential Linear Programming (SLP) coupled with a global convergence filter and an adaptive move-limit strategy. The algorithm can handle linear equality, linear in-equality, non-linear in-equality, and non-linear equality constraints. A merit function approach ensures unconditional feasibility of the linearized sub problem. Here, the user can specify two algorithms/apparoches to ensure feasibility. Default is the Merit approach presented in [1]. Alternativly, the Augmented Lagrange method is also available. This method updates the penalty parameters/Lagrange multipliers in each iteration, see [2] for details. The final Lagrange multipliers are available for the user to analyze. 

The global convergence filter monitors progression of the penalized objective function (the merit function), and the associated non-linear constraints. Based on the response, the filter algorithm adjusts the move-limits to ensure stable convergence. The convergence filter is based [3] and [4]

The adaptive move-limit strategy controls the box-constraints (upper and lower bounds for the design variables) and is based on the work by professor Erik Lund from Aalborg University.

The overall framework is based on an implementation developed during my Ph.d. studies, see [1]
Please refer to [1] when citing the fminslp algorithm.

[1] R Soerensen, E Lund (2015): Thickness filters for gradient based multi-material and thickness optimization of laminated composite structures, Structural and Multidisciplinary Optimization 52 (2), 227-250
https://doi.org/10.1007/s00158-015-1230-3.

[2]  Nocedal J, Wright SJ: Numerical Optimization, second edition, ISBN-10:0-387-303003-0, p 514.
 
[3] Chin CM, Fletcher R (1999): On the global convergence of an SLP-
filter algorithm that takes EQP steps. Numerical Analysis Report
NA/199, Department of Mathematics, University of Dundee,Scotland, UK

[4] Fletcher R, Leyffer S, Toint PL (1998): On the global convergence
of an SLP-filter algorithm. Numerical Analysis Report NA/183,
Department of Mathematics, University of Dundee, Scotland, UK

Usage: see also examples.m 

myProblem = fminslp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon); % Default options are applied

myProblem = fminslp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options); % Input options structure

myProblem = fminslp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,'MoveLimit',0.01,'MaxIterations',100); % Input/modify individual options

[x,fval,exitflag,output] = myProblem.solve(); % Solve problem

Available options:

Initialize values to default

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
	SpecifyConstraintGradient: 0
	 SpecifyObjectiveGradient: 0
	           CheckGradients: 0
	 FiniteDifferenceStepSize: sqrt(eps)
		FiniteDifferenceType : 'forward' or 'backward' or 'central'
                
Modify individual options, the rest are initialized to default values

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
	SpecifyConstraintGradient: 0
	 SpecifyObjectiveGradient: 0
	           CheckGradients: 0
	 FiniteDifferenceStepSize: sqrt(eps)
		FiniteDifferenceType : 'forward'

Change log
* v1.1
  * Changed default behaivor wrt., user supplied gradients. Now, the program assumes no user supplied gradients. Same as with fmincon.If the user has analytical gradients, these can be supplied and activated by using the two options 'SpecifyConstraintGradient' and 'SpecifyObjectiveGradient'
  * Added finite difference schemes to approximate objective and constraint gradients
  * Added option to check user supplied gradients agains finite difference approximations
  * Found bug in examples
