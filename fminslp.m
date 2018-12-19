classdef fminslp < handle
  
  properties(Constant)
    name = 'fminslp';
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
    
    % Iteration history
    history = struct();
    
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
      
      % We made it this far
      this.initialized = true;
      
      % initialize options structure
      this.options = fminslp.slpoptions(varargin);
      
      % Initialize global convergence filter
      this.filter = fminslp.initializeGlobalConvergenceFilter(this.options);
      
    end
    
    % Main function
    function [x,fval,exitflag,output] = solve(this)
      % Assume the code failed
      exitflag = -1;
      
      % Allocate iteration history array
      % Store function values, maximum infeasibility from non-linear
      % constraints, norm of design variable change
      this.history.f = zeros(this.options.MaxIterations,1);
      this.history.xnorm = zeros(this.options.MaxIterations,1);
      
      
      x = this.x0(:);
      
      % Initialize the current box constraints with upper and lower bounds.
      xLcur = this.lb(:);
      xUcur = this.ub(:);
      
      % Initialize "old" design variables
      xold1 = x;
      xold2 = x;
      
      % evaluate objective function at initial point
      this.f0  = this.fun(x);
      
      % Get scaling factor for merit function
      this.aFac = max([abs(this.f0),1]);
      
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
        this.history.maxInf = zeros(this.options.MaxIterations,1);
      else
        this.nGnl = 0;
        % Assumed allocated to empty when not active
        y = [];
        ylb = [];
        yub = [];
        this.history.maxInf = [];
      end
      
      % Evaluate merit function
      fmerit = this.getMeritObj(x,y);
      % Store initial value for filter
      this.filter.initF = fmerit;
      % Store "old" value for convergence check
      fOld = fmerit;
      
      % Set counters and switches
      nFval = 1;
      iterNo = 0;
      optimize = true;
      % Main loop
      while optimize
        
        % update iteration counter
        iterNo = iterNo + 1;
        
        % evaluate gradients
        [~,dfmerit] = this.getMeritObj(x,y);
        [gmerit] = this.getMeritConstraints(x,y);
        [~,dgmerit] = this.getMeritConstraints(x,y);

        if ~isempty(this.A)
          % Assemble linear and non-linear in-equality constraints
          Am = zeros(size(this.A,1)+this.nGnl,size(this.A,2)+this.nGnl); % Allocate
          Am(1:size(this.A,1),1:size(this.A,2)) = this.A; % Load linear
          bm = zeros(size(this.b,1)+this.nGnl,1); % Allocate
          bm(1:size(this.b,1)) = this.b; % Load linear
          
          if ~isempty(this.nonlcon)
            Am(size(this.A,1)+1:end,:) = dgmerit; % Load non-linear
            bm(size(this.b,1)+1:end) = dgmerit*[x;y]-gmerit; % Load non-linear
          end
        else % Only non-linear
          if ~isempty(this.nonlcon)
            Am = dgmerit;
            bm = dgmerit*[x;y]-gmerit;
          else
            Am = [];
            bm = [];
          end
        end

        % Expand linear equality constraints to account for slag variables
        if ~isempty(this.Aeq)
          Ameq = zeros(size(this.Aeq,1),size(this.Aeq)+this.nGnl);
          Ameq(1:size(this.Aeq,1),1:size(this.Aeq,2)) = this.Aeq;
        else
          Ameq = [];
        end
        
        % update move-limits
        reduceSwitch = false;
        [xLcur, xUcur] = this.AdaptiveMoveLimit(x, xLcur, xUcur, this.lb, this.ub, this.options.MoveLimit ,this.options.MoveLimitReduce, this.options.MoveLimitExpand, xold1, xold2, reduceSwitch);
        
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
          [xlpNew,exitflag] = this.lpSolver(dfmerit,Am,bm,Ameq,this.beq,xL,xU,[x;y]);
          
          % Split design variables
          xNew = xlpNew(1:this.nDV);
          yNew = xlpNew(this.nDV+1:end);
          
          deltaxy = xlpNew-[x;y];
          deltax = xNew - x;
          deltanorm = norm(deltax);
          
          
          optimalityNorm = norm(dfmerit(1:end-this.nGnl)) + norm(yNew);
          
          % evaluate new design
          nFval = nFval + 1;
          fmerit = this.getMeritObj(xNew,yNew);
          gmerit = this.getMeritConstraints(xNew,yNew);
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
            [xLcur, xUcur] = this.AdaptiveMoveLimit(x, xLcur, xUcur,this.lb, this.ub, this.options.MoveLimit ,this.options.MoveLimitReduce, this.options.MoveLimitExpand, xold1, xold2, reduceSwitch);
          end
          
          % check for convergence
          if (optimalityNorm <= this.options.OptimalityTolerance) || (deltanorm <=this.options.StepTolerance) || (iterNo >= this.options.MaxIterations) || (nFval >= this.options.MaxFunctionEvaluations)
            optimize = false;
            backtrack = false;
            exitflag = 1;
            fval = this.fun(x);
          end
          
          
        end % Inner loop
        
        % Does the new point(h,f) qualify to be added to the filter?
        if (AddToFilter)
          [ this.filter ] = this.UpdateFilter(this.filter, this.filter.h, this.filter.f);
        end
        
        % Update design variables
        x = xNew;
        y = yNew;
        % Update "old" design
        fOld = fmerit;
        this.history.f(iterNo) = fmerit;
        this.history.xnorm(iterNo) = deltanorm;
        this.history.maxInf(iterNo) = this.filter.h;
        this.history.nIter = iterNo;
        this.history.nFeval = nFval;
      end % Main loop
      
      this.history.f(iterNo+1:end)=[];
      this.history.xnorm(iterNo+1:end)=[];
      this.history.maxInf(iterNo+1:end)=[];
      output = this.history;
      
    end % Solve function
    
    function postprocess(this)
      % Save current "default" window style
      defaultWindowStyle=get(0,'DefaultFigureWindowStyle');
      % Set new window style to docked
      set(0,'DefaultFigureWindowStyle','docked')
      
      % Make iteration vector
      ivec = 1:this.history.nIter;
      
      figure();
      plot(ivec,this.history.f)
      title('Objective')
      xlabel('Iteration Number')
      ylabel('Objective value')
      
      figure();
      plot(ivec,this.history.xnorm)
      title('Design change norm')
      xlabel('Iteration Number')
      yl=ylabel('|x_i-x_(i-1)| value');
      set(yl,'Interpreter','none')
      
      figure();
      plot(ivec,this.history.maxInf)
      title('Maximum infeasibility')
      xlabel('Iteration Number')
      ylabel('-')
      
      % Restore default window style
      set(0,'DefaultFigureWindowStyle',defaultWindowStyle)
    end
    
  end % methods
  
  
  methods (Hidden = true)
    
    function [fmerit,dfmerit] = getMeritObj(this,x,y)
      fmerit = [];
      dfmerit = [];
      switch nargout
        case 1
          fval = this.fun(x);
          fmerit = fval + this.aFac*sum(y*this.options.InfeasibilityPenalization+0.5*y.^2);
        case 2
          [~,df] = this.fun(x);
          dy = this.aFac*(this.options.InfeasibilityPenalization+y);
          dfmerit = [df;dy];
      end
    end
    
    function [gmerit,dgmerit] = getMeritConstraints(this,x,y)
      gmerit = [];
      dgmerit = [];
      if ~isempty(this.nonlcon)
        switch nargout
          case 1
            [gn, gneq] = this.nonlcon(x);
            gmerit = [gn(:);gneq(:);-gneq(:)]-y;
          case 2
            [~,~,dgnl,dgneq] = this.nonlcon(x);
            dgmerit = zeros(this.nGnl,this.nDV+this.nGnl);
            dgmerit(:,1:this.nDV) = [dgnl';dgneq';-dgneq'];
            dgmerit(:,this.nDV+1:end) = -1.0*eye(this.nGnl);
        end
      end
    end
    
    function [x,exitflag] = lpSolver(this,df,A,b,Aeq,beq,lb,ub,x)
      switch this.options.Solver
        
        case 'linprog'
          linprogOptions = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
          [x,~, exitflag]= linprog(df,A,b,Aeq,beq,lb,ub,x,linprogOptions);
          
         case 'glpk' % Octave lp solver
          
          % Define constraints
          if ~isempty(A) && ~isempty(Aeq)
            nleq = size(A,1);
            neq = size(Aeq,1);
            ctype = char(nleq+neq,1);
            ctype(1:nleq) = 'U';
            ctype(nleq+1:end) = 'S';
            A = [A;Aeq];
            b = [b;beq];
          elseif ~isempty(A)
            nleq = size(A,1);
            ctype = char(nleq,1);
            ctype(1:nleq) = 'U';
          elseif ~isempty(Aeq)
            neq = size(Aeq,1);
            ctype = char(neq,1);
            ctype(1:neq) = 'S';
            A = Aeq;
            b = beq;
          else
            ctype = [];
          end
          
          % Define variable type
          ndv = numel(df);
          vartype = char(ndv,1);
          vartype(1:ndv) = 'C';
          
          [x, ~, exitflag] = glpk (df, A, b, lb, ub, ctype, vartype);
        otherwise
          error([this.name,' Unknown LP solver specified: ',this.options.Solver])
      end
    end
    
  end
  
  methods(Static = true, Hidden = true)
    
    % options initialization
    function options = slpoptions(input)
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));

      
      % Define offset methods
      p.addParameter('Algorithm','merit',  @(x) checkEmpetyOrChar(x));
      p.addParameter('Solver','linprog',  @(x) checkEmpetyOrChar(x));
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MaxIterations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('InfeasibilityPenalization',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('OptimalityTolerance',1e-6,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-10,  @(x) checkEmptyOrNumericPositive(x));
      
      % Move-limit parameters
      p.addParameter('MoveLimitMethod','adaptive',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MoveLimit',0.1,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MoveLimitExpand',1.1,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MoveLimitReduce',0.5,  @(x) checkEmptyOrNumericPositive(x));
      
      
      % Define global convergene filter parameters
      p.addParameter('MaxInfeasibility',inf,  @(x) checkEmptyOrNumericPositive(x));
      
      % pars input
      if isempty(input)
        parse(p);
      else
        parse(p,input{:});
      end
      
      options = p.Results;
      
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
      filter.vals(1,1) = inf;  % Maximum infeasibility
      filter.vals(1,2) = inf;  % Related obj value
    end
    
    % SLP Global Convergence filter function
    function [filter] = EvaluateCurrentDesignPointToFilter(filter)

    h = filter.h;
    f = filter.f;

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
    function [xLcur, xUcur] = AdaptiveMoveLimit(x, xLcur, xUcur, xLorg, xUorg, moveLimit ,reduceFac, expandFac, xold1, xold2, reduceSwitch)

      if (reduceSwitch)
        Expand = reduceFac;
        Reduction = reduceFac;
      else
        Reduction = reduceFac;
        Expand = expandFac;
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
        if (xLcur(dvNo) >= xUcur(dvNo)); xLcur(dvNo) = 0.9999999*xUcur(dvNo); end;
        if (xUcur(dvNo) <= xLcur(dvNo)); xUcur(dvNo) = 1.0000001*xLcur(dvNo); end;
      end
    end
    
  end
  
end