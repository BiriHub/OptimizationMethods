function [x,points,i] = CG(A,b,x0,max_itr,tol)
% Conjugate Gradient Method

% A is a matrix semi positive definite 
% b is the right-hand side column vector of the system Ax =b
% x0 is the starting point
% max_itr is the number of max iterations for the iterative process
% tol is the distance for two consecutive points below which we stop

%starting point
x = x0;

%row vector of all points computed during the iterative process
points = zeros(max_itr,size(x,1)); 


points(1,:)=x';
%residual of the starting point
r = A*x -b;
p = -r;
i=1;
while norm(r)>tol && i<max_itr 
    
    %Compute the step size 
    alpha = -(dot(r,p) / dot(A*p, p));         
    
    %Update the position
    x = x + alpha * p;                      
    
    %Update the residual
    r = r + alpha * (A*p);                  

    %Update the beta coefficient
    beta = dot(A*r,p) / dot(A*p, p);
    
    %Update the step direction
    p= - r  + beta * p;
    
    %Save the current point
    points(i+1,:)=x';
    
    i=i+1;
end
i = i-1;
end
