%% 1.5 

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

delta_max = 2;
eta = 0.22; 
max_iter = 10000; %Maximum number of iterations
tol = 1e-6; % Tollerance for reaching the convergence

% Delta parameters definition
delta1=0.05;
delta2= 0.50;
delta3=1.5;


%% Cauchy point method 

% delta1
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta1, eta, max_iter, tol,false);


%% Cauchy point energy landscape plot for delta1

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method for \Delta_1');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);


% delta2
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta2, eta, max_iter, tol,false);


%% Cauchy point energy landscape plot for delta2

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method for \Delta_2');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);



% delta3
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta3, eta, max_iter, tol,false);


%% Cauchy point energy landscape plot for delta3

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Cauchy method for \Delta_3');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);


%% Dogleg method


% delta1
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta1, eta, max_iter, tol,true);


%% Dogleg energy landscape plot for delta1

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Dogleg method for \Delta_1');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);



% delta2
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta2, eta, max_iter, tol,true);


%% Dogleg energy landscape plot for delta2

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Dogleg method for \Delta_2');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);



% delta3
[vecx, iter, vecgrad] = trustRegion(func, func_grad, func_hessian, x0, delta_max, delta3, eta, max_iter, tol,true);


%% Dogleg energy landscape plot for delta3

% Define the range for the plot
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y));
title('Energy landscape of Dogleg method for \Delta_3');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= vecx(1,1:iter);
y_points= vecx(2,1:iter); 
plot(x_points, y_points, 'ro-', 'MarkerSize', 4);


%% Comment - 1.5
% delta=0.01
% Cauchy number of iterations : 10001
% Dogleg number of iterations : 24
% delta2=0.50
% Cauchy number of iterations : 10001
% Dogleg number of iterations : 21
% delta3=1.5
% Cauchy number of iterations : 10001
% Dogleg number of iterations : 19

% It can be notice that apparently changing the value of delta does not really bring
% to better improvements of the final result for the Cauchy point method given that it does not converge in 10000 iterations. 
% On the other hand the Dogleg method seems to be slightly influenced by
% the changing of the radius of the trust region area , indeed when the
% delta is very far from the maximum delta (2), the method tends to increase the number of iterations (24
% iterations), but when the delta is near to the delta_max value it seems
% to improve the convergence saving few iterations (19 iterations).