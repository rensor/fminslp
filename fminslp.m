classdef fminslp
  
  properties(Constant)
    name = 'fminslp';
    version = 'v1.4';
  end

  properties
    options = [];
  end
  
  properties(SetAccess = public, Hidden = true)
    
    % Inputs
    fun = [];
    x0 = [];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = [];
    nDV = [];
    
    % Merit function variables
    f0 = []; % Initial, unpeanalized objective function value
    nGnl = []; % number of non-linear inequality constraints
    aFac = []; % Scaling parameter for merit function
    
    % Global convergence filter 
    filter = struct();
    
    % Switch
    initialized = false;
    
    
  end
  
  methods
    
    % Construct
    function this = fminslp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,varargin)
      
      % initialize data structures
      if nargin < 1 || ~isa(fun,'function_handle')
        error([this.name,': fun is required to be a function handle'])
      else
        this.fun = fun;
      end
      
      if nargin < 2 || isempty(x0)
        this.x0 = [];
      else
        if isnumeric(x0)
          this.x0 = x0(:);
        else
          error([this.name,' x0 is required to be of type numeric'])
        end
      end
      
      if nargin < 3 || isempty(A)
        this.A = [];
      else
        if isnumeric(A)
          this.A = A;
        else
          error([this.name,' A is required to be of type numeric'])
        end
      end
      
      if nargin < 4 || isempty(b)
        this.b = [];
      else
        if isnumeric(b)
          this.b = b(:);
        else
          error([this.name,' b is required to be of type numeric'])
        end
      end
      
      if size(this.A,1) ~= size(this.b,1)
        error([this.name,' A and b must contain equal number of rows'])
      end
      
      if nargin < 5 || isempty(Aeq)
        this.Aeq = [];
      else
        if isnumeric(Aeq)
          this.Aeq = Aeq;
        else
          error([this.name,' Aeq is required to be of type numeric'])
        end
      end
      
      if nargin < 6 || isempty(beq)
        this.beq = [];
      else
        if isnumeric(beq)
          this.beq = beq(:);
        else
          error([this.name,' beq is required to be of type numeric'])
        end
      end
      
      if size(this.Aeq,1) ~= size(this.beq,1)
        error([this.name,' Aeq and beq must contain equal number of rows'])
      end
          
      if nargin < 7 || isempty(lb)
        this.lb = [];
      else
        if isnumeric(lb)
          this.lb = lb(:);
        else
          error([this.name,' lb (lower bound) is required to be of type numeric'])
        end
      end
      
      if nargin < 8 || isempty(ub)
        this.ub = [];
      else
        if isnumeric(ub)
          this.ub = ub(:);
        else
          error([this.name,' ub (upper bound) is required to be of type numeric'])
        end
      end
      
      if nargin < 9 || isempty(nonlcon)
        this.nonlcon = [];
      elseif ~isempty(nonlcon)
        if ~isa(nonlcon,'function_handle')
          error([this.name,' nonlcon is required to be a function handle'])
        else
          this.nonlcon = nonlcon;
        end
      end
      
      % check that all sizes match
      if ~isempty(this.x0)
        this.nDV = numel(this.x0);
      end
      
      if ~isempty(this.lb) && ~isempty(this.ub)
        if numel(this.lb) ~= numel(this.ub)
          error([this.name,' lb and ub not equal dimensions'])
        end
      end
      
      if ~isempty(this.lb) && ~isempty(this.nDV)
        if numel(this.lb) ~= this.nDV
          error([this.name,' x0 and lb not equal dimensions'])
        end
      end
      
      if ~isempty(this.ub) && ~isempty(this.nDV)
        if numel(this.ub) ~= this.nDV
          error([this.name,' x0 and ub not equal dimensions'])
        end 
      end
      
      if ~isempty(this.A) && ~isempty(this.nDV)
        if size(this.A,2) ~= this.nDV
          error([this.name,' Columns of A(',num2str(size(this.A,2)),') does not match number of design variables(',num2str(this.nDV),')'])
        end
      elseif ~isempty(this.A) && isempty(this.nDV)
        this.nDV = size(this.A,2);
      end
      
      if ~isempty(this.Aeq) && ~isempty(this.nDV)
        if size(this.Aeq,2) ~= this.nDV
          error([this.name,' Columns of Aeq(',num2str(size(this.A,2)),') does not match number of design variables(',num2str(this.nDV),')'])
        end
      elseif ~isempty(this.Aeq) && isempty(this.nDV)
        this.nDV = size(this.Aeq,2);
      end
      
      % initialize options structure
      this.options = fminslp.slpoptions(varargin);
      
      % Initialize global convergence filter
      this.filter = fminslp.initializeGlobalConvergenceFilter(this.options);
      
      % Check Gradients
      if this.options.CheckGradients
        this.CheckUserSuppliedGradients;
      end
      
      
      % We made it this far
      this.initialized = true;
    end
    
    % Main function
    function [x,fval,exitflag,output] = solve(this)
      % Assume the code failed
      exitflag = -1;
      
      
    % Output strucutre
    output = struct('iterations',[],...
                    'funcCount',[],...
                    'constrviolation',[],...
                    'firstorderopt',[],...
                    'message','',...
                    'iterHistory',struct());
                  
      
      if strcmpi(this.options.Display,'iter')
        s2 = sprintf('fminslp optimizer with global convergence filter');
        s3 = sprintf(' %10s      %10s      %10s      %10s      %10s','f(x)','Max inf', 'Norm dx', 'nFeval','IterNo');
        nAst = length(s3);
        s1 = repmat('*',[1,nAst]);
        nwS = round((nAst-length(s2))/2);
        s2center = sprintf('%s%s',repmat(' ',[1,nwS]),s2);
        fprintf('%s\n%s\n%s\n%s\n',s1,s2center,s1,s3);
        
      end
        
      % Allocate iteration history array
      % Store function values, maximum infeasibility from non-linear
      % constraints, norm of design variable change
      output.iterHistory.f = zeros(this.options.MaxIterations,1);
      output.iterHistory.xnorm = zeros(this.options.MaxIterations,1);
      
      % Set minimum "box" limit around each design variables. 
      minDVBoxLimit = 2*this.options.StepTolerance;
      
      % Load initial x value
      x = this.x0(:);
      
      % Initialize the current box constraints with upper and lower bounds.
      xLcur = this.lb(:);
      xUcur = this.ub(:);
      
      % Initialize "old" design variables
      xold1 = x;
      xold2 = x;
      
      % evaluate initial non-linear in-equality constraints
      if ~isempty(this.nonlcon)
        [gnl, gnleq] = this.nonlcon(x);
        % Determine number of constraints
        this.nGnl = numel(gnl) + numel(gnleq)*2;
        % Initialize slag variables to close the infeasibility gap
        y = max([gnl(:);gnleq(:);-gnleq(:)],0);
        % Set upper and lower bounds for slag variables
        ylb = zeros(this.nGnl,1);
        yub = repmat(max([gnl(:);gnleq(:);-gnleq(:);1])*1e6,[this.nGnl,1]);
        % Allocate iteration history for maximum infeasibility
        output.iterHistory.maxInf = zeros(this.options.MaxIterations,1);
      else
        this.nGnl = 0;
        % Assumed allocated to empty when not active
        y = [];
        ylb = [];
        yub = [];
      end
      
      % Expand linear equality constraints to account for slag variables
      if ~isempty(this.Aeq)
        Ameq = zeros(size(this.Aeq,1),size(this.Aeq,2)+this.nGnl);
        Ameq(1:size(this.Aeq,1),1:size(this.Aeq,2)) = this.Aeq;
      else
        Ameq = [];
      end
      
      if ~isempty(this.A)
        % Assemble linear in-equality constraints, and expand design variables/columns to
        % and account for slag variables from the non-linear in-equality constraints
        Am = zeros(size(this.A,1)+this.nGnl,size(this.A,2)+this.nGnl); % Allocate
        Am(1:size(this.A,1),1:size(this.A,2)) = this.A; % Load linear
        bm = zeros(size(this.b,1)+this.nGnl,1); % Allocate
        bm(1:size(this.b,1)) = this.b; % Load linear
      end
      
      if strcmpi(this.options.Algorithm,'AL')
        % Initialize lagrange multipliers to 1, 
        % this is nesseary such that the gradient for the y-variable
        % is at least 1 in the first iteration.
        lambda = ones(this.nGnl,1);
        % Update lagrange multipliers based on the initial infeasibilities
        lambda = this.updateLambda(y,lambda);
      else
        lambda = [];
      end
      
      % evaluate objective function at initial point
      this.f0  = this.fun(x);
      % Get scaling factor for merit function
      this.aFac = max([abs(this.f0),1]);
      
      % Evaluate merit function
      fmerit = this.getMeritObj(x,y,lambda);
      % Store initial value for filter
      this.filter.initF = fmerit;
      % Store "old" value for convergence check
      fOld = fmerit;
      
      % Set counters and switches
      nFeval = 1;
      iterNo = 0;
      optimize = true;
      % Main loop
      while optimize
        
        % update iteration counter
        iterNo = iterNo + 1;
        
        % evaluate gradients
        [~,~,dfmerit] = this.getMeritObj(x,y,lambda);
        [gmerit,~,dgmerit] = this.getMeritConstraints(x,y);
        
        % Setup in-equality constraints to include non-linear parts
        if ~isempty(this.nonlcon) && ~isempty(this.A)
          Am(size(this.A,1)+1:end,:) = dgmerit; % Load non-linear
          bm(size(this.b,1)+1:end) = dgmerit*[x;y]-gmerit; % Load non-linear
        elseif ~isempty(this.nonlcon) && isempty(this.A)
          Am = dgmerit;
          bm = dgmerit*[x;y]-gmerit;
        else
          Am = [];
          bm = [];
        end

        % update move-limits
        reduceSwitch = false;
        [xLcur, xUcur] = this.AdaptiveMoveLimit(x, xLcur, xUcur, this.lb, this.ub, this.options.MoveLimit ,this.options.MoveLimitReduce, this.options.MoveLimitExpand, xold1, xold2, reduceSwitch,minDVBoxLimit);
        
        % update old values
        xold2 = xold1;
        xold1 = x;
        
        % Set switches
        backtrack = true;
        
        % Inner loop
        while backtrack
          
          % Set lower and upper bounds for the lp problem
          xL = [xLcur(:); ylb];
          xU = [xUcur(:); yub];
        
          % Call optimizer
          [xlpNew,exitflag,lpMessage] = this.lpSolver(dfmerit,Am,bm,Ameq,this.beq,xL,xU,[x;y]);
          
          % Terminate if solution failed
          if exitflag~=1
            optimize = false;
            backtrack = false;
            output.message = sprintf([output.message,'\n ',lpMessage],'s');
          end
          
          % Split design variables
          xNew = xlpNew(1:this.nDV);
          yNew = xlpNew(this.nDV+1:end);
          
          % Determine design deltas
          deltaxy = xlpNew-[x;y];
          deltax = xNew - x;
          deltanorm = norm(deltax);
          
          % Determine optimality Norm for convergence check
          optimalityNorm = norm(dfmerit(1:end-this.nGnl));
          
          % evaluate constraints
          [gmerit,greal] = this.getMeritConstraints(xNew,yNew);
          
          % Evaluate objective function at new point
          nFeval = nFeval + 1;
          
          [fmerit,freal] = this.getMeritObj(xNew,yNew,lambda);
          
          % Determine delta for merit function
          deltaf = fOld - fmerit;
          
          % Determine maximum infeasibility for current design
          this.filter.h = max([gmerit;0]); 
          % Store objective function value for current design
          this.filter.f = fmerit/this.filter.initF;
          
          % Evaluate current (h,f) point against the convergence filter
          [this.filter] = this.EvaluateCurrentDesignPointToFilter(this.filter);
          
          % Assume that we don't want to update the convergence filter
          AddToFilter = false;
          if (this.filter.PointAcceptedByFilter)
            % Determine Delta values for filter checks
            deltaL = -dfmerit'*deltaxy;
            % Check if we should add accept the current point
            if ( (deltaf<this.filter.sigma*deltaL) && (deltaL>0.0) )
              % Not accepted
              reduceSwitch = true;
            else
              % Accepted
              reduceSwitch = false;
              backtrack = false;
              % Check if we should but the new point (h,f) into the filter
              if(this.filter.h > 0.0)
                AddToFilter = true;
              end
            end
          else
            reduceSwitch = true;
          end
          
          if reduceSwitch
            [xLcur, xUcur] = this.AdaptiveMoveLimit(x, xLcur, xUcur,this.lb, this.ub, this.options.MoveLimit ,this.options.MoveLimitReduce, this.options.MoveLimitExpand, xold1, xold2, reduceSwitch,minDVBoxLimit);
          end
          
          % check for convergence
          if (optimalityNorm <= this.options.OptimalityTolerance)
            optimize = false;
            backtrack = false;
            exitflag = 1;
            output.message = sprintf([output.message,'Sucessfully solve to Optimality Tolerance: <= %0.5e'],this.options.OptimalityTolerance);
          elseif abs(deltaf)<=this.options.FunctionTolerance 
            optimize = false;
            backtrack = false;
            exitflag = 2;
            output.message = sprintf([output.message,'Sucessfully solve to Function Tolerance: <= %0.5e'],this.options.FunctionTolerance);
          elseif (deltanorm <=this.options.StepTolerance) 
            optimize = false;
            backtrack = false;
            exitflag = 3;
            output.message = sprintf([output.message,'Sucessfully solve to Step Tolerance: <= %0.5e'],this.options.StepTolerance);
          elseif (fmerit <=this.options.ObjectiveLimit) && iterNo > 1
            optimize = false;
            backtrack = false;
            exitflag = 4;
            output.message = sprintf([output.message,'Sucessfully solve to Objective Limit: <= %0.5e'],this.options.ObjectiveLimit);            
          elseif (iterNo >= this.options.MaxIterations) 
            optimize = false;
            backtrack = false;
            exitflag = 0;
            output.message = sprintf([output.message,'Number of iterations exceeded the limit: %i'],this.options.MaxIterations);
          elseif (nFeval >= this.options.MaxFunctionEvaluations)
            optimize = false;
            backtrack = false;
            exitflag = 0;
            output.message = sprintf([output.message,'Number of function evaluations exceeded the limit: %i'],this.options.MaxFunctionEvaluations);
          end
          
        end % Inner loop
        
        % Does the new point(h,f) qualify to be added to the filter?
        if (AddToFilter)
          [ this.filter ] = this.UpdateFilter(this.filter, this.filter.h, this.filter.f);
        end
        
        % Update design variables
        x = xNew;
        y = yNew;
        
        if strcmpi(this.options.Algorithm,'AL')
          % Update lagrange multipliers based on the initial infeasibilities
          lambda = this.updateLambda(y,lambda);
        end
        
        % Update "old" design
        fOld = fmerit;
        % Store history
        maxInf = max([greal;0]);
         
        output.iterHistory.f(iterNo) = freal;
        output.iterHistory.xnorm(iterNo) = deltanorm;
        output.iterHistory.constrviolation(iterNo) = maxInf;
        
        
        if strcmpi(this.options.Display,'iter')
            fprintf(' %6.4e      %6.4e      %6.4e      %10i      %10i \n' ,freal, maxInf, deltanorm, nFeval ,iterNo);
        end
        
      end % Main loop
      
      % Prepare output
      fval = freal;
      
      output.iterHistory.f(iterNo+1:end)=[];
      output.iterHistory.xnorm(iterNo+1:end)=[];
      output.iterHistory.constrviolation(iterNo+1:end)=[];
      
      output.constrviolation = maxInf;
      output.iterations = iterNo;
      output.funcCount = nFeval;
      output.firstorderopt = optimalityNorm;
      
      if strcmpi(this.options.Algorithm,'AL')
        output.lambda = lambda./this.aFac;
      end
      
    end % Solve function
    
    function postprocess(this)
      % Save current "default" window style
      defaultWindowStyle=get(0,'DefaultFigureWindowStyle');
      % Set new window style to docked
      set(0,'DefaultFigureWindowStyle','docked')
      
      % Make iteration vector
      ivec = 1:this.output.iterHistorynIter;
      
      f1=figure();
      plot(ivec,this.output.iterHistoryf)
      title('Objective')
      xlabel('Iteration Number')
      ylabel('Objective value')
      
      figure();
      plot(ivec,this.output.iterHistoryxnorm)
      title('Design change norm')
      xlabel('Iteration Number')
      yl=ylabel('Norm dx');
      set(yl,'Interpreter','none')
      
      figure();
      plot(ivec,this.output.iterHistorymaxInf)
      title('Maximum infeasibility')
      xlabel('Iteration Number')
      ylabel('-')
      
      % Jump back to figure 1
      figure(f1)
      % Restore default window style
      set(0,'DefaultFigureWindowStyle',defaultWindowStyle)
    end
    
  end % methods
  
  
  methods (Hidden = true)
    
    function [fmerit,fval,dfmerit] = getMeritObj(this,x,y,lambda)
      
      if nargout > 2
        if this.options.SpecifyObjectiveGradient
          [fval,df] = this.fun(x);
        else
          fval = this.fun(x);
          df = this.getFunDSA(x,fval);
        end
        
        switch this.options.Algorithm
          case 'merit'
            dy = this.aFac*(this.options.InfeasibilityPenalization+y);
          case 'al'
            dy = this.aFac*(lambda+this.options.InfeasibilityPenalization*y);
        end
        dfmerit = [df;dy];
      else
        fval = this.fun(x);
      end

      switch this.options.Algorithm
        case 'merit'
          fmerit = fval + this.aFac*sum(y*this.options.InfeasibilityPenalization+0.5*y.^2);
        case 'al'
          fmerit = fval + this.aFac*sum(y.*lambda+this.options.InfeasibilityPenalization*0.5*y.^2);
      end

    end
    
    function [gmerit,greal,dgmerit] = getMeritConstraints(this,x,y)
      gmerit = [];
      dgmerit = [];
      greal = [];
      if ~isempty(this.nonlcon)
        if nargout > 2
          if this.options.SpecifyConstraintGradient
            [gn,gneq,dgnl,dgneq] = this.nonlcon(x);
          else
            [gn, gneq] = this.nonlcon(x);
            [dgnl,dgneq] = this.getNonlconDSA(x);
          end
          dgmerit = zeros(this.nGnl,this.nDV+this.nGnl);
          dgmerit(:,1:this.nDV) = [dgnl';dgneq';-dgneq'];
          dgmerit(:,this.nDV+1:end) = -1.0*eye(this.nGnl);
        else
          [gn, gneq] = this.nonlcon(x);
        end
          greal = [gn(:);gneq(:);-gneq(:)];
          gmerit = greal-y;
      end
    end
    
    function lambda = updateLambda(this,y,lambda)
      lambda = max(lambda + this.options.InfeasibilityPenalization*y,1);
    end
    
    function [x,exitflag,message] = lpSolver(this,df,A,b,Aeq,beq,lb,ub,x)
      % Here you can add your own solvers
      switch this.options.Solver
        
        case 'linprog' % MathWorks solver
          linprogOptions = optimoptions('linprog','Algorithm','dual-simplex','Display','off','ConstraintTolerance',this.options.StepTolerance/2);
          [x,~, exitflag,lpOutput]= linprog(df,A,b,Aeq,beq,lb,ub,x,linprogOptions);
          message = lpOutput.message;
        case 'clp'
          % Select open source solver from OPTI tool box (see ToolSet\BladeLib\Ext\OptiToolbox\Solvers)
          rhsLower = [-inf(size(A,1),1);beq];
          rhsUpper = [b;beq];
          [x,~,exitflag,lpOutput] = opti_clp([],df,[A;Aeq],rhsLower,rhsUpper,lb,ub);
          message = lpOutput.Status;
        case 'glpk' % Octave lp solver
          
          % Define constraints
          neq = size(Aeq,1);
          nleq = size(A,1);
          ctype = [repmat('U',[nleq,1]);repmat('S',[neq,1])];
          gA = [A;Aeq];
          gb = [b;beq];
          
          % Define variable type
          ndv = numel(df);
          vartype = repmat('C',[ndv,1]);
          % Change default solver settings
          gParam = struct('dual',2); % Use two-phase dual simplex, and if it fails, switch
                                     % to the primal simplex.
          [x, ~, exitflagOut,lpOutput] = glpk(df, gA, gb, lb, ub, ctype, vartype,1,gParam);
          % exitflag = 1 is a success
          % exitflag < 0 is a fail
          if exitflagOut == 0
            exitflag = 1;
          else
            % All error codes are positive from glpk
            exitflag = -exitflagOut;
          end
          % Convert output status to a message
          switch lpOutput.status
            case 1
              message = 'Solution status is undefined';
            case 2
              message = 'Solution is feasible';
            case 3
              message = 'Solution is infeasible';
            case 4
              message = 'Problem has no feasible solution';
            case 5
              message = 'Solution is optimal';
            case 6
              message = 'Problem has no unbounded solution';
          end
            
        otherwise
          error([this.name,' Unknown LP solver specified: ',this.options.Solver])
      end
    end
    
    function [df] = getFunDSA(this,x,fin)
      df = zeros(this.nDV,1);
      h = this.options.FiniteDifferenceStepSize;
      switch this.options.FiniteDifferenceType
        case 'forward'
          if nargin < 3 || ~isempty(fin)
            fin = this.fun(x);
          end
          for dvNo = 1:this.nDV
            xin = x(dvNo);
            x(dvNo) = xin+h;
            p1 = this.fun(x);
            df(dvNo) = (p1-fin)/(h);
            x(dvNo) = xin;
          end
        case 'backward'
          if nargin < 3 || ~isempty(fin)
            fin = this.fun(x);
          end
          for dvNo = 1:this.nDV
            xin = x(dvNo);
            x(dvNo) = xin-h;
            m1 = this.fun(x);
            df(dvNo) = (fin-m1)/(h);
            x(dvNo) = xin;
          end
        case 'central'
          for dvNo = 1:this.nDV
            xin = x(dvNo);
            x(dvNo) = xin+h;
            p1 = this.fun(x);
            x(dvNo) = xin-h;
            m1 = this.fun(x);
            df(dvNo) = (p1-m1)/(2*h);
            x(dvNo) = xin;
          end
        otherwise
          error([this.name,': Unknown FiniteDifferenceType'])
      end
      
      
    end
    
    function [dg,dgeq] = getNonlconDSA(this,x)
      
      h = this.options.FiniteDifferenceStepSize;
      dg = [];
      dgeq = [];
      switch this.options.FiniteDifferenceType
        case 'forward'
          [gnl0, gnleq0] = this.nonlcon(x);
          ngnl = numel(gnl0);
          ngnleq = numel(gnleq0);
          if (ngnl>0)
            dg = zeros(this.nDV,ngnl);
          end
          
          if (ngnleq>0)
            dgeq = zeros(this.nDV,ngnleq);
          end
          if ngnl > 0 || ngnleq > 0
            for dvNo = 1:this.nDV
              xin = x(dvNo);
              x(dvNo) = xin+h;
              [gp1, geqp1] = this.nonlcon(x);
              for gNo = 1:ngnl
                dg(dvNo,gNo) = (gp1(gNo)-gnl0(gNo))/(h);
              end
              for gNo = 1:ngnleq
                dgeq(dvNo,gNo) = (geqp1(gNo)-gnleq0(gNo))/(h);
              end
              x(dvNo) = xin;
            end
          end
        case 'backward'
          [gnl0, gnleq0] = this.nonlcon(x);
          ngnl = numel(gnl0);
          ngnleq = numel(gnleq0);
          if (ngnl>0)
            dg = zeros(this.nDV,ngnl);
          end
          
          if (ngnleq>0)
            dgeq = zeros(this.nDV,ngnleq);
          end
          if ngnl > 0 || ngnleq > 0
            for dvNo = 1:this.nDV
              xin = x(dvNo);
              x(dvNo) = xin-h;
              [gm1, geqm1] = this.nonlcon(x);
              for gNo = 1:ngnl
                dg(dvNo,gNo) = (gnl0(gNo)-gm1(gNo))/(h);
              end
              for gNo = 1:ngnleq
                dgeq(dvNo,gNo) = (gnleq0(gNo)-geqm1(gNo))/(h);
              end
              x(dvNo) = xin;
            end
          end
        case 'central'
          firstTime = true;
          for dvNo = 1:this.nDV
            xin = x(dvNo);
            x(dvNo) = xin+h;
            [gp1, geqp1] = this.nonlcon(x);
            x(dvNo) = xin-h;
            [gm1, geqm1] = this.nonlcon(x);
            
            if firstTime
              firstTime = false;
              ngnl = numel(gp1);
              ngnleq = numel(geqp1);
              if (ngnl>0)
                dg = zeros(this.nDV,ngnl);
              end
              
              if (ngnleq>0)
                dgeq = zeros(this.nDV,ngnleq);
              end
            end
            if ngnl > 0 || ngnleq > 0
              for gNo = 1:ngnl
                dg(dvNo,gNo) = (gp1(gNo)-gm1(gNo))/(2*h);
              end
              for gNo = 1:ngnleq
                dgeq(dvNo,gNo) = (geqp1(gNo)-geqm1(gNo))/(2*h);
              end
            end
            x(dvNo) = xin;
          end
        otherwise
          error([this.name,': Unknown FiniteDifferenceType'])
      end
      
    end
    
    function CheckUserSuppliedGradients(this)
          
          this.options.FiniteDifferenceType = 'central';
          if this.options.SpecifyObjectiveGradient
            [~,dfUser] = this.fun(this.x0);
            dfFiniteDiff = this.getFunDSA(this.x0);
            dsaDiff = dfUser-dfFiniteDiff;
            maxDiff = max(abs(dsaDiff));
            formatSpec = '\n \t Derivative Check Information\n Objective function derivatives: \n Maximum difference between user-supplied and finite-difference derivatives = %0.5e \n';
            fprintf(formatSpec,maxDiff);
          end
          if this.options.SpecifyConstraintGradient
            [~,~,dgnlUser,dgneqUser] = this.nonlcon(this.x0);
            [dgnlFiniteDiff,dgneqFiniteDiff] = this.getNonlconDSA(this.x0);
            if ~isempty(dgnlUser)
              dsaDiff = abs(dgnlUser-dgnlFiniteDiff);
              [maxDiff,idx]=max(dsaDiff(:));
              [dvNo,gNo]=ind2sub(size(dsaDiff),idx);
              formatSpec = '\n \t Derivative Check Information\n Nonlinear inequality constraint derivatives: \n Maximum difference between user-supplied and finite-difference derivatives = %0.5e \n \t User-supplied constraint derivative element (%i,%i): %0.5e \n \t Finite-difference constraint derivative element (%i,%i): %0.5e \n';
              fprintf(formatSpec,maxDiff,dvNo,gNo,dgnlUser(dvNo,gNo),dvNo,gNo,dgnlFiniteDiff(dvNo,gNo));
            end
            if ~isempty(dgneqUser)
              dsaDiff = abs(dgneqUser-dgneqFiniteDiff);
              [maxDiff,idx]=max(dsaDiff(:));
              [gNo,dvNo]=ind2sub(size(dsaDiff),idx);
              formatSpec = '\n \t Derivative Check Information\n Nonlinear equality constraint derivatives: \n Maximum difference between user-supplied and finite-difference derivatives = %0.5e \n \t User-supplied constraint derivative element (%i,%i): %0.5e \n \t Finite-difference constraint derivative element (%i,%i): %0.5e \n';
              fprintf(formatSpec,maxDiff,dvNo,gNo,dgnlUser(dvNo,gNo),dvNo,gNo,dgnlFiniteDiff(dvNo,gNo));
            end
          end
    end
    
  end
  
  methods (Static = true, Hidden = true)
    
    % options initialization
    function options = slpoptions(input)
      % Here you can add new options if needed
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
      checkNumericPositive = @(x) ((isnumeric(x) && all(x > 0)));
      checkLogicalZeroOne = @(x) (islogical(x) || ((x==1) || (x==0)));
      
      % Set parameters
      p.addParameter('Algorithm','merit',  @(x) checkEmpetyOrChar(x));
      p.addParameter('Solver','linprog',  @(x) checkEmpetyOrChar(x));
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      
      % Convergence parameters
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MaxIterations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('InfeasibilityPenalization',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('OptimalityTolerance',1e-6,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('FunctionTolerance',1e-6,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-8,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('ObjectiveLimit',-1e20,  @(x) isnumeric(x));
      
      % Gradient parameters
      p.addParameter('FiniteDifferenceType','forward',  @(x) ischar(x));
      p.addParameter('FiniteDifferenceStepSize',sqrt(eps),@(x) checkNumericPositive(x));
      p.addParameter('SpecifyConstraintGradient',0,@(x) checkLogicalZeroOne(x));
      p.addParameter('SpecifyObjectiveGradient',0,@(x) checkLogicalZeroOne(x)); 
      p.addParameter('CheckGradients',0,@(x) checkLogicalZeroOne(x)); 
      
      % Move-limit parameters
      p.addParameter('MoveLimitMethod','adaptive',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MoveLimit',0.1,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MoveLimitExpand',1.1,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MoveLimitReduce',0.5,  @(x) checkEmptyOrNumericPositive(x));
      
      % Define global convergene filter parameters
      p.addParameter('MaxInfeasibility',inf,  @(x) checkEmptyOrNumericPositive(x));
      
      % pars input
      if nargin < 1 || isempty(input)
        parse(p);
      else
        parse(p,input{:});
      end
      
      % Output results to options structure
      options = p.Results;
      % Enforce lower case
      options.Algorithm = lower(options.Algorithm);
      options.Solver = lower(options.Solver);
      options.FiniteDifferenceType = lower(options.FiniteDifferenceType);
      options.MoveLimitMethod = lower(options.MoveLimitMethod);
      
    end
    
    function filter = initializeGlobalConvergenceFilter(options)
      % SLP Filter settings
      filter.SmallVal = 1.0e-6; % Used for pertubation from zero and one
      filter.gamma = filter.SmallVal;
      filter.beta = 1.0 - filter.SmallVal;
      filter.sigma = 2*filter.SmallVal;
      filter.delta = filter.SmallVal;
      filter.vals = zeros(options.MaxIterations,2); % Filter is initialized
      filter.nVals = 1; % Number of points in the filter
      filter.PointAcceptedByFilter = false; % Is the current filter point accepted 
      filter.h = 1.0e30; % Current infeasibility point, large initialization
      filter.f = 1.0e30; % Current objective point, large initialization
      filter.initF = 0.0; % Initial objective function value
      
      % Initial values in SLP filter
      if options.MaxInfeasibility < inf
        % Here MaxInfeasibility sets the maximum acceptable infeasibility
        % wrt., the global constraints (merit constrints)
        % Low values e.g., 1e-10 can potentially slow the convergence rate
        % and/or result in the optimizer getting trapped in a local minimum
        
        filter.vals(1,1) = options.MaxInfeasibility; % Maximum infeasibility
        filter.vals(1,2) = -inf;  % Related obj value
      else
        % Here, the convergence filter will steadily add new points to the
        % filter to ensure stable convergence
        
        filter.vals(1,1) = inf;  % Maximum infeasibility
        filter.vals(1,2) = inf;  % Related obj value
      end
    end
    
    % SLP Global Convergence filter function
    function [filter] = EvaluateCurrentDesignPointToFilter(filter)
    % Extract current point
    h = filter.h;
    f = filter.f;
    % Loop for all values in filter
      for ii = 1:filter.nVals
        hi = filter.vals(ii,1);
        fi = filter.vals(ii,2);
        if ( (h <= hi*filter.beta) || (f+filter.gamma*h) <= fi )
          filter.PointAcceptedByFilter = true;
        else
          filter.PointAcceptedByFilter = false;
          break
        end % Accepted or not
      end % Loop for filter values
    end
    
    % SLP Global Convergence filter function
    function [newFilter] = UpdateFilter(filter, hk, fk)
    % This function adds the pair (hk,fk) to the filter.
    % In this process, it determines wheater the new point dominates 
    % some of the points already in the filter. If so, it removes these values
    % and adds the new point to the filter.
      newFilter = filter;
      newFilter.nVals = 0;
      newFilter.vals(:,:) = 0;
      Update = true;
      ii = 0;
      if (filter.nVals >= 1)
        while (Update)
            ii = ii + 1;
            hi = filter.vals(ii,1);
            fi = filter.vals(ii,2);
          % Dominated or not?
          if ( (hk <= hi) && (fk <= (fi)) )
            
          else % Copy old data to new filter
            newFilter.nVals = newFilter.nVals + 1;
            newFilter.vals(newFilter.nVals,1) = hi;
            newFilter.vals(newFilter.nVals,2) = fi;
          end
          
          if (ii >= filter.nVals)
            Update = false;
          end
        end 
      end
      
      % Add new values to filter
      newFilter.nVals = newFilter.nVals + 1;
      newFilter.vals(newFilter.nVals,1) = hk;
      newFilter.vals(newFilter.nVals,2) = fk;
    end
    
    % Adaptive move-limit algorithm
    function [xLcur, xUcur] = AdaptiveMoveLimit(x, xLcur, xUcur, xLorg, xUorg, moveLimit ,reduceFac, expandFac, xold1, xold2, reduceSwitch, minDVBoxLimit)
      
      if (reduceSwitch)
        Expand = reduceFac;
        Reduction = reduceFac;
      else
        Reduction = reduceFac;
        Expand = expandFac;
      end
      
      if nargin < 12 || isempty(minDVBoxLimit)
        minDVBoxLimit = 0;
      end
      
      nDv = numel(x);
      for dvNo = 1:nDv
        delta = (xUcur(dvNo)-xLcur(dvNo))/2; % This was the previous allowable change
        % Use the iteration history to determine whether we have oscillations
        % in the design variables
        if (abs(x(dvNo)-xold1(dvNo)) > 1.e-10)
          s1 = (xold1(dvNo)-xold2(dvNo)) / (x(dvNo)-xold1(dvNo));
          if (s1 < 0.0)
            delta = delta*Reduction;      % oscillation, slow increase
          else
            delta = delta*Expand;       % Stable, allow more rapid increase
          end
        else
          delta = delta*moveLimit;
        end
        dmax = (xUorg(dvNo)-xLorg(dvNo))*moveLimit;
        if (delta > dmax) 
          delta = dmax;
        end
        % Initial extimate of lower and upper bound on x(i)
        xLcur(dvNo) = x(dvNo) - delta;
        xUcur(dvNo) = x(dvNo) + delta;
        % Make sure we are within the feasible domain
        xLcur(dvNo) = max(xLcur(dvNo),xLorg(dvNo));
        xUcur(dvNo) = min(xUcur(dvNo),xUorg(dvNo));
        
        % Take care of extremely small design changes where the bounds may be interchanged 
        if (xLcur(dvNo) > xUcur(dvNo))
          if minDVBoxLimit > 0
            xLcur(dvNo) = xUcur(dvNo)-minDVBoxLimit/2;
          else
            xLcur(dvNo) = (1-1e-6)*xUcur(dvNo);
          end
        end
        if (xUcur(dvNo) < xLcur(dvNo))
          if minDVBoxLimit > 0
            xUcur(dvNo) = xLcur(dvNo)-minDVBoxLimit/2;
          else
            xUcur(dvNo) = (1+1e-6)*xLcur(dvNo);
          end
        end
        
        % Enforce the minimum box limit size
        if (xUcur(dvNo)-xLcur(dvNo)) < minDVBoxLimit
          if (xUcur(dvNo) <= xUorg(dvNo)-minDVBoxLimit/2) && (xLcur(dvNo) <= xLorg(dvNo)-minDVBoxLimit/2)
            xUcur(dvNo) = xUcur(dvNo) + minDVBoxLimit/2;
            xLcur(dvNo) = xLcur(dvNo) - minDVBoxLimit/2;
          elseif (xUcur(dvNo) <= xUorg(dvNo)-minDVBoxLimit)  
            xUcur(dvNo) = xUcur(dvNo) + minDVBoxLimit;
          elseif (xLcur(dvNo) >= xLorg(dvNo)-minDVBoxLimit)
            xLcur(dvNo) = xLcur(dvNo) - minDVBoxLimit;
          end
        end
      end
    end
    
  end
  
end