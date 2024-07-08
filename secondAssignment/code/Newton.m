function [x,energy_val,gradient_norms,i] = Newton(x0,alpha,f,f_gradient,f_hessian,max_itr,tol,useBackTracking)
    %% Newton method 
    % x0 is the starting point
    % alpha is the value of the step size
    % f is the function on which apply the method
    % f_gradient is a vector containing the partial derivatives of the
    % function 'f'
    % f_hessian is a matrix 2x2 containing the second partial derivatives of the function 'f'
    % max_itr is the number of max iterations for the iterative process
    % tol is the distance for two consecutive points below which we stop
    % useBackTracking is a boolean value used to flag the use of the backtracking algorithm for computing the step size value

    i=1; %iteration
    x=x0;

    energy_val = zeros(max_itr, 3); %energy function
    gradient_norms = zeros(max_itr, 1); %norms of gradient

    %Save the information related to the starting point
    energy_val(1,1:3)=[x,f(x(1),x(2))];
    gradient_norms(1)=norm(f_gradient(x(1),x(2)));
     
    % Iterative process
    while norm(f_gradient(x(1),x(2)))>tol && i<max_itr 
        
        grad_values=f_gradient(x(1),x(2));
        step_direction=-(f_hessian(x(1),x(2))\grad_values);
                
        if(useBackTracking)        
            alpha= backTracking(f,f_gradient,x,step_direction,1);
        end
        
        %Update the point coordinates
        x = x + alpha * step_direction';  
    
        %Save the information of the current point
        energy_val(i+1,:)=[x,f(x(1),x(2))];
        gradient_norms(i+1)=norm(f_gradient(x(1),x(2)));
        i=i+1; 
    end    
end

