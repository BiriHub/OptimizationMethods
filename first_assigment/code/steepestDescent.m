function [x,energy_val,gradient_norms,i] = steepestDescent(A,x0,mu,func,gradient,max_itr,tol)
    %% Steepest descent iterative method
    % A is the matrix used in the quadratic form necessary for computing
    % the optimal step length at each iteration of the steepest method
    % x_in is the starting point
    % func is the function on where compute the method
    % gradient is a vector containing the partial derivatives of the
    % function 'func'
    % max_itr is the number of max iterations for the iterative process
    % tol is the distance for two consecutive points below which we stop

    
    i=1; 
    x=x0;
    energy_val = zeros(max_itr, 3);
    gradient_norms = zeros(max_itr, 1);

    energy_val(1,1:2)=x; %save x,y of the point
    energy_val(1,3)=func(x(1),x(2),mu); %save f(x) of the starting point
    

    gradient_norms(1)=norm(gradient(x,mu),2);
    
    % Steepest descent
    while norm(gradient(x,mu))>tol && i<max_itr
        grad_value = gradient(x,mu)';
        
        %compute optimal alpha
        numerator=dot(grad_value,grad_value); 
        alpha = numerator / dot(grad_value*A,grad_value );    % optimal alpha

        x = x - alpha * grad_value;                     

        energy_val(i+1,1:2)=x; %save coordinates x and y of the point
        energy_val(i+1,3)=func(x(1),x(2),mu); %save f(x) of the current point
    
        %Compute the gradient norm
        gradient_norms(i+1)=norm(gradient(x,mu)',2);
    

        i=i+1;
    end
    gradient_norms(i+1)=norm(gradient(x,mu)',2);
    
end

