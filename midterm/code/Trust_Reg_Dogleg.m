function [vecx, itr, vecgrad] = Trust_Reg_Dogleg(f, grad, hess, x0, delta_max, ...
    delta0, eta, maxitr, tol)

%   Trust region method using Dogleg method
%   Reference: J. Nocedal and S. Wright, Numerical optimization, page 69
%   Algorithm 4.1
%   Input:
%        * f = function
%        * grad = gradient of function
%        * hess = Hessian matrix of function
%        * x0 = starting point
%        * delta_max = maximum value of delta (must be > 0)
%        * delta0 = starting delta in interval (0, delta_max)
%        * eta = parameter in interval [0, 1/4)
%        * maxitr = maximum number of iterations allowed
%        * tol = tolerance threshold

x = x0;
delta = delta0;
itr = 1;

% Vector to save iterates
vecx = zeros(size(x0, 1), maxitr + 1);
vecx(:, itr) = x0;

% Vector to save norm of gradients
vecgrad = zeros(1, maxitr + 1);

for i = 1 : maxitr
    g = grad(x(1), x(2));
    h = hess(x(1), x(2));
    % Vector to save norm of gradients
    vecgrad(itr) = norm(g);

    if norm(g) < tol 
        break;
    end
    
    % Dogleg method
    pB = -(hessian\gradient);

    if norm(pB)<=delta
        p=pB;
        return;
    else
        pV= -((gradient'*gradient)/(gradient'*hessian*gradient))*gradient;
        if (norm(pV)>delta)
            p = (pV/norm(pV))*delta;
        else
            tau= chooseTau(pB,pV,delta);
            p = pV + tau*(pB-pV);
        end
    end

    f_k=f(x(1),x(2));
    m = f_k + p'*g+(1/2)*p'*h*p;

    rho = (f_k-f(x(1)+p(1),x(2)+p(2)))/(f_k-m);

     if rho < 1 / 4
        delta = (1/4)*delta;

    elseif rho > 3 / 4 && norm(p) == delta
        delta = min(2*delta,delta_max);
    end
    
    if rho > eta
        x = x + p;
    end
    
    itr = i + 1;
    
    % Vector to save iterates
    vecx(:, itr) = x;
    
end

vecx = vecx(:, 1:itr);
vecgrad = vecgrad(:, 1:itr);

end