function [x,fval,exitflag,output,ex] = examples(No)
  
  switch No
    case 1
      % exmple 1
      %	Minimize	f(x1,x2) = -2x1 -x2
      %
      %	S.t.  g1(x1,x2) = x1^2 + x2^2 <= 25   (nonlinear inequality constraint)
      %         g2(x1,x2) = x1^2 - x2^2 <= 7    (nonlinear inequality constraint)
      %         0 <= x1 <= 10;  0 <= x2 <= 10   (box constraints)
      %
      %   Starting point: (x1,x2) = (1,1) (a feasible point)
      x0 = [1;1];
      lb = [0;0];
      ub = [7;10];
      ex = fminslp(@f1,x0,[],[],[],[],lb,ub,@g1);
      [x,fval,exitflag,output] = ex.solve;
      ex.postprocess;
    
    case 2
      % exmple 2  
      %	Minimize	f(x1,x2) = x1^4 - 2*x1*x1*x2 + x1*x1 + x1*x2*x2 - 2*x1 + 4
      %
      %	S.t. h(x1,x2) = x1^2 + x2^2 = 2            (nonlinear equality constraint)
      %        g(x1,x2) = 1/4*x1^2 +3/4*x2^2 <= 1  (nonlinear inequality constraint)
      %        0 <= x1 <= 4;  0 <= x2 <= 4         (box constraints)
      %
      %   Starting point: (x1,x2) = (sqrt(2),0)    (a feasible point)
      x0 = [sqrt(2);0];
      lb = [0;0];
      ub = [4;4];
      ex = fminslp(@f2,x0,[],[],[],[],lb,ub,@g2);
      [x,fval,exitflag,output] = ex.solve;
      ex.postprocess;

  end

  function [f,df] = f1(x)
    f = [];
    df = [];

    switch nargout
      case 1
        f = -2*x(1) - x(2);
      case 2
        df(1,1) = -2;
        df(2,1) = 1;
    end
  end

  function [g,geq,dg,deg] = g1(x)
    g = [];
    geq = [];
    dg = [];
    deg = [];

    switch nargout
      case {1,2}
        g(1)= x(1)^2 + x(2)^2 - 25;
        g(2) = x(1)^2 - x(2)^2 - 7;
      case {3,4}
        dg(1,1) = 2*x(1);
        dg(2,1) = 2*x(2);
        dg(1,2) = 2*x(1);
        dg(2,2) = -2*x(2);
    end

  end

  function [f,df] = f2(x)
    f = [];
    df = [];

    switch nargout
      case 1
        f = x(1)^4 - 2*x(1)^2*x(2) + x(1)^2 + x(1)*x(2)^2 - 2*x(1) + 4;
      case 2
        df(1,1) = 4*x(1)^3 -4*x(1)*x(2) + 2*x(1) + x(2)^2 -2;
        df(2,1) = -2*x(1)^2 + 2*x(1);
    end
  end

  function [g,geq,dg,deg] = g2(x)
    g = [];
    geq = [];
    dg = [];
    deg = [];

    switch nargout
      case {1,2}
        g(1)= 1/4*x(1)^2 + 3/4*x(2)^2-1;
      case {3,4}
        dg(1,1) = 2/4*x(1);
        dg(2,1) = 6/4*x(2);
    end

  end
end