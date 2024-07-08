function [x,energy_val,gradient_norms,i] = BFGS(x0,H,f,f_gradient,max_iter,tol)
    %% BFGS iterative method
    % x0 is the starting point
    % H is the initial hessian matrix
    % f is the function on which apply the method
    % f_gradient is a vector containing the partial derivatives of the
    % function 'f'
    % max_itr is the number of max iterations for the iterative process
    % tol is the distance for two consecutive points below which we stop

    
    i=1; % counter of iterations
    x=x0;
    I=eye(size(x,2)); %identity matrix
    initial_alpha= 1; % required to the backtracking algorithm

    energy_val = zeros(max_iter, 3); %energy function
    gradient_norms = zeros(max_iter, 1); %norms of gradient

    %Save the information related to the starting point
    energy_val(1,1:3)=[x,f(x(1),x(2))];
    gradient_norms(1)=norm(f_gradient(x(1),x(2)));
     

    while norm(f_gradient(x(1),x(2)))>tol && i<max_iter 
        
        %Computing the gradient of function
        grad_value = f_gradient(x(1),x(2));
        
        % Compute search direction
        step_direction = -H * grad_value;

        % Compute the alpha through the backtracking
        alpha= backTracking(f,f_gradient,x,step_direction,initial_alpha);

        xk= x; %save previous point 

        %Update the point coordinates
        x = x + alpha * step_direction';

        %Update the hessian matrix
        s = x' - xk';
        y = f_gradient(x(1),x(2))-grad_value;
        
        rho= 1/dot(y,s);


        H = (I-rho*s*y')*H*(I-rho*y*s')+rho*(s*s');

        %Save the information of the current point
        energy_val(i+1,:)=[x,f(x(1),x(2))];
        gradient_norms(i+1)=norm(f_gradient(x(1),x(2)));
        i=i+1; 
    end

end
