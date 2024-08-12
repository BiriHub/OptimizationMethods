%% 1.4 

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



delta0 = 0.75;
delta_max = 2;
eta = 0.22; 
max_iter = 10000; %Maximum number of iterations
tol = 1e-6; % Tollerance for reaching the convergence

%% Cauchy point method

x1=[0;0];
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x1, delta_max, delta0, eta, max_iter, tol,false);


%% Cauchy point energy landscape plot for x1

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method for x1');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

x2=[-1;-1];
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x2, delta_max, delta0, eta, max_iter, tol,false);


%% Cauchy point energy landscape plot for x2

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method for x2');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

x3=[-3;-3];
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x3, delta_max, delta0, eta, max_iter, tol,false);


%% Cauchy point energy landscape plot for x3

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method for x3');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

%% Dogleg point method - 1.4

x1=[0;0];
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x1, delta_max, delta0, eta, max_iter, tol,true);


%% Cauchy point energy landscape plot for x1

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Dogleg method for x_1');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

x2=[-1;-1];
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x2, delta_max, delta0, eta, max_iter, tol,true);


%% Cauchy point energy landscape plot for x2

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Dogleg method for x_2');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

x3=[-3;-3];
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x3, delta_max, delta0, eta, max_iter, tol,true);


%% Dogleg method energy landscape plot for x3

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method for x_3');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);

%% Comment - 1.4
% x1 = [0;0]
% Cauchy number of iterations : 1926
% Dogleg number of iterations : 16
% x2= [-1;-1]
% Cauchy number of iterations : 10001
% Dogleg number of iterations : 22
% x3=[-3;-3]
% Cauchy number of iterations : 725
% Dogleg number of iterations : 25

% It can be notice that the Dogleg method easily reaches the minimum of the function in just
% few steps in each starting point compared to the Cauchy method approach
% which one requires many iterations to converge. In particular , it can be
% notice that the method does not reach the minimun in 10000 iteration when
% the starting point is [-1;-1]

