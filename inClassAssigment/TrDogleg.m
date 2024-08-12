function [p] = TrDogleg(hessian,gradient,delta)
 % Dogleg method
    % Input
    % - hessian : hessian matrix (for Newton method) or approximated
    % hessian matrix (Quasi-Newton method)
    % - gradient : value of gradient of the function
    % - delta : radius of the trust region area

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
end

