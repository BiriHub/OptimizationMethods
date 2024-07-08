function [alpha] = backTracking(f,f_gradient,x,p,init_alpha)
    %% Backtracking algorithm
    % f is the function on where compute the method
    % f_gradient is a vector containing the partial derivatives of the
    % function 'f'
    % x is the starting point
    % p is the step direction of the iterative process
    % init_alpha initial value of alpha

    %Variables definition
    alpha=init_alpha; % initial alpha
    fka = f(x(1)+alpha*p(1), x(2)+alpha*p(2)); % lhs Wolfe 1st condition
    fk = f(x(1),x(2)); % A part of the rhs of 1st Wolfe condition
    gp = p'*f_gradient(x(1),x(2)); % dot product between gradient and direction, to complete rhs 1st Wolfe condition
    

    i = 0;
    rho = 0.9;
    c1 = 1e-4; % the constant factor generally is quite small
    
    % While loop to guarantee the Wolfe 1° condition
    while(fka > fk + c1*alpha*gp && i < 1000) 
        alpha = rho*alpha; % reduce the alpha
        fka = f(x(1)+alpha*p(1), x(2)+alpha*p(2)); % new 1st Wolfe condition
        i = i + 1;
    end
end

