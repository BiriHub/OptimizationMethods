%% Exercise n° 1


%% Variables declaration

syms x y;

%Rosenbrock’s function
rosenbrock_f = power((1-x),2)+100*power((y-x.^2),2);
variables = symvar(rosenbrock_f);

%Function
func = matlabFunction(rosenbrock_f);

%Gradient of function
func_grad= matlabFunction(gradient(rosenbrock_f, variables));

%Hessian matrix of function (required to the Newton method)
func_hessian = matlabFunction(hessian(rosenbrock_f, variables));

%Starting point
x0 = [0,0]; 

alpha = 1; %default step size
max_iter = 50000; %Maximum number of iterations
tol = 1e-6; % Tollerance for reaching the convergence

%% Exercise n° 1.2 - Steepest method

%% GD without backtracking application
[x,energy_val,gradient_norms,iter] = GD(x0,alpha,func,func_grad,max_iter,tol,false);

%%Plots
printPlots(func,energy_val,gradient_norms,iter);

%% GD with backtracking application
[x,energy_val,gradient_norms,iter] = GD(x0,alpha,func,func_grad,max_iter,tol,true);

%%Plots
printPlots(func,energy_val,gradient_norms,iter);

%% Exercise 1.3 - Newton method

%% Newton method without backtracking application
[x,energy_val,gradient_norms,iter] = Newton(x0,alpha,func,func_grad,func_hessian,max_iter,tol,false);

%%Plots 
%Note : I sum the first and last value of "grad_norms" and "energy_val"
%with "tol" value given that log function does not accept the 0 value inside its domain.
%To show the gradient norm and the function values I used the tollerance value as a lower bound 


% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape 2D');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= energy_val(1:iter,1);
y_points= energy_val(1:iter,2);

plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

%Plot the logarithm of the gradient norm at each iteration of the iterative
%process

if gradient_norms(1)==0 
    gradient_norms(1)=gradient_norms(1)+tol;
end
if gradient_norms(iter)==0
    gradient_norms(iter)=gradient_norms(iter)+tol;
end

figure;
semilogy(0:iter-1, gradient_norms(1:iter), 'b-');
title('Gradient norms along iterations');
xlabel('number of iterations');
ylabel('||\nabla f||');

hold on;
plot(0,gradient_norms(1),'bo'); % print the first iteration on plot
plot(iter-1,gradient_norms(iter),'bo'); % print the last iteration on plot



% Plot the function value at each iteration of the iterative
%process

if energy_val(1,3)==0 
    energy_val(1,3)=energy_val(1,3)+tol;
end
if energy_val(iter,3)==0
    energy_val(iter,3)=energy_val(iter,3)+tol;
end

figure;
semilogy(0:iter-1,energy_val(1:iter,3), 'b-', 'MarkerSize', 4);
title('Function values along iterations');
xlabel('number of iterations');
ylabel('f(x)');

hold on;
plot(0,energy_val(1,3),'bo'); % print the first iteration on plot
plot(iter-1,energy_val(iter,3),'bo'); % print the last iteration on plot


%% Newton method with backtracking application
[x,energy_val,gradient_norms,iter] = Newton(x0,alpha,func,func_grad,func_hessian,max_iter,tol,true);

%%Plots 
printPlots(func,energy_val,gradient_norms,iter);