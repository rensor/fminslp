function [x,fval,exitflag,output,ex] = examples(No)
  
  switch No
    case 1
      % exmple 1
      %	Minimize	f(x1,x2) = -2x1 -x2
      %
      %	S.t.  g1(x1,x2) = x1^2 + x2^2 <= 25       (nonlinear inequality constraint)
      %         g2(x1,x2) = x1^2 - x2^2 <= 7      (nonlinear inequality constraint)
      %         0 <= x1 <= 10;  0 <= x2 <= 10     (box constraints)
      %
      %   Starting point: (x1,x2) = (1,1)         (a feasible point)
      x0 = [1;1];
      lb = [0;0];
      ub = [7;10];
      ex = fminslp(@f1,x0,[],[],[],[],lb,ub,@g1,'display','iter');
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
      ex = fminslp(@f2,x0,[],[],[],[],lb,ub,@g2,'display','iter');
      [x,fval,exitflag,output] = ex.solve;
      ex.postprocess;
    case 3
      % Example 12.3 from Arora (2012) Example 7.3, pp. 285-287 (Arora (2004) pp. 418-420):
      %
      %	Minimize	f(x1,x2) = (x1-10)^3 + (x2-20)^3
      %
      %	S.t.  g1(x1,x2) = -(x1-5)^2 -(x2-5)^2 <= -100     (nonlinear inequality constraint)
      %         g2(x1,x2) = (x1-6)^2 + (x2-5)^2 <= 82.81  (nonlinear inequality constraint)
      %         13 <= x1 <= 100;  0 <= x2 <= 100          (box constraints)
      %
      %   Starting point: (x1,x2) = (20,5)                (a infeasible point) 
      x0 = [13;5];
      lb = [13;0];
      ub = [100;100];
      ex = fminslp(@f3,x0,[],[],[],[],lb,ub,@g3,'display','iter');
      [x,fval,exitflag,output] = ex.solve;
      ex.postprocess;
      
    case 4
      %	Minimize	max (f1(x1,x2) = 7 - 2*x1 + 2*x2; f2(x1,x2) = 5 - x1 + 3*x2)
      %
      %	S.t.  g1(x1,x2) = x1^2 + x2^2 <= 25     (nonlinear inequality constraint)
      %         g2(x1,x2) = x1^2 - x2^2 <= 7    (nonlinear inequality constraint)
      %         0 <= x1 <= 10;  0 <= x2 <= 10   (box constraints)
      %
      %   Starting point: (x1,x2) = (1,1)       (a feasible point) 
      x0 = [1;1;10];
      lb = [0;10;0];
      ub = [0;10;1e6];
      ex = fminslp(@f4,x0,[],[],[],[],lb,ub,@g4,'display','iter');
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

  function [g,geq,dg,dgeq] = g1(x)
    g = [];
    geq = [];
    dg = [];
    dgeq = [];

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

  function [g,geq,dg,dgeq] = g2(x)
    g = [];
    geq = [];
    dg = [];
    dgeq = [];

    switch nargout
      case {1,2}
        g(1)= 1/4*x(1)^2 + 3/4*x(2)^2-1;
        geq(1) = x(1)^2 + x(2)^2 - 2;
      case {3,4}
        dg(1,1) = 2/4*x(1);
        dg(2,1) = 6/4*x(2);
        dgeq(1,1) = 2*x(1);
        dgeq(2,1) = 2*x(2);
    end

  end

  function [f,df] = f3(x)
    f = [];
    df = [];

    switch nargout
      case 1
        f = (x(1)-10)^3 + (x(2)-20)^3;
      case 2
        df(1,1) =3*(x(1)-10)^2;
        df(2,1) = 3*(x(2)-20)^2;
    end
  end

  function [g,geq,dg,dgeq] = g3(x)
    g = [];
    geq = [];
    dg = [];
    dgeq = [];

    switch nargout
      case {1,2}
        g(1)= -(x(1)-5)^2/100 - (x(2)-5)^2/100 + 1;
        g(2) = (x(1)-6)^2/82.81 + (x(2)-5)^2/82.81 - 1;
      case {3,4}
        dg(1,1) = -2*(x(1)-5)/100;
        dg(2,1) = -2*(x(2)-5)/100;
        dg(1,2) = 2*(x(1)-6)/82.81;
        dg(2,2) = -2*(x(2)-5)/82.81;

    end

  end

  function [f,df] = f4(x)
    f = [];
    df = [];

    switch nargout
      case 1
        f = x(3);
      case 2
        df(1,1) = 0;
        df(2,1) = 0;
        df(3,1) = 1;
    end
  end

  function [g,geq,dg,dgeq] = g4(x)
    g = [];
    geq = [];
    dg = [];
    dgeq = [];

    switch nargout
      case {1,2}
        g(1) =  7 - 2*x(1) + 2*x(2) - x(3);
        g(2) =  5 - x(1) + 3*x(2) - x(3);
        g(3)= (x(1)^2 + x(2)^2)/25 - 1;
        g(4) = (x(1)^2 - x(2)^2)/7 - 1;
      case {3,4}
        dg(1,1) = -2;
        dg(2,1) = 2;
        dg(3,1) = -1;
        
        dg(1,2) =  -1;
        dg(2,2) =  3;
        dg(3,2) = -1;
        
        dg(1,3) = 2*(x(1))/25;
        dg(2,3) = 2*(x(2))/25;
        dg(3,3) = 0;
        
        dg(1,4) = 2*(x(1))/7;
        dg(2,4) = -2*(x(2))/7;
        dg(3,4) = 0;
    end

  end
end