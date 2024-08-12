%% Exercise 1 - in class assignment

%

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
x0 = [0;-2]; 

delta0 = 0.75;
delta_max = 2;
eta = 0.22; 
max_iter = 10000; %Maximum number of iterations
tol = 1e-6; % Tollerance for reaching the convergence

[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta0, eta, max_iter, tol,false);


%% Cauchy point energy landscape plot

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);


%% Dogleg method 
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta0, eta, max_iter, tol,true);


%% Dogleg energy landscape plot

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Dogleg method');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

%% Comment - 1.3
% Cauchy point method number of iterations : 9651
% Dogleg method number of iterations : 21
% It is evident from the final results that the Dogleg method results to be more efficient in
% reaching the minimum of the Rosembrock's function compared to the Cauchy
% method


%% Bonus - 1.6