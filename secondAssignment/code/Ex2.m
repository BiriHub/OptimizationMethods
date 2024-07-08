%% Exercise n°2 - BFGS

%% Variables declaration

syms x y;

%Rosenbrock’s function
rosenbrock_f = power((1-x),2)+100*power((y-x.^2),2);
variables = symvar(rosenbrock_f);

%Function
func = matlabFunction(rosenbrock_f);

%Gradient of function
func_grad= matlabFunction(gradient(rosenbrock_f, variables));



%Starting point
x0 = [0,0]; 

alpha = 1; %default step size
max_iter = 500; %Maximum number of iterations
tol = 1e-6; % Tollerance for reaching the convergence

H = eye(2); %initial hessian matrix

[x,energy_val,gradient_norms,iter] = BFGS(x0,H,func,func_grad,max_iter,tol);

%%Plots
printPlots(func,energy_val,gradient_norms,iter);
