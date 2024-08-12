function [p] = TrCauchy(hessian,gradient,delta)
 % Cauchy point method
    % Input
    % - hessian : hessian matrix (for Newton method) or approximated
    % hessian matrix (Quasi-Newton method)
    % - gradient : value of gradient of the function
    % - delta : radius of the trust region area

    p = (-gradient/norm(gradient))*delta;

    alpha=1;

    if p'*hessian*p >0
        tao = -(p'*gradient)/(p'*hessian*p);
        alpha = min(1,tao);
    end

    p=p*alpha;
end

