%% Exercise n° 3

%% 3.2


% Define the range for the plot
x = linspace(-10, 10, 100);
y = linspace(-10, 10, 100);
[X, Y] = meshgrid(x, y);


% Compute the function values
z_mu_1 = X.^2 + 1 * Y.^2; % For mu = 1
z_mu_10 = X.^2 + 10 * Y.^2; % For mu = 10



% mu =1

% surface plot

figure;
% subplot(1, 2, 1);
surf(X, Y, z_mu_1);
title('Surface Plot for \mu = 1');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');

%contour plot
figure;
contour(X, Y, z_mu_1, 20); % 20 contour levels
title('Contour Plot for \mu = 1');
xlabel('x');
ylabel('y');

% mu =10 

%surface plot

figure;

surf(X, Y, z_mu_10);
title('Surface Plot for \mu = 10');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');

%contour plot
figure;
contour(X, Y, z_mu_10, 20); % 20 contour levels
title('Contour Plot for \mu = 10');
xlabel('x');
ylabel('y');


%% 3.4 and 3.5

%Variables declaration


N = 100; % Max number of iterations
epsilon = 1e-8;

f = @(x, y, mu) x.^2 + mu * y.^2; 
f_gradient = @(x,mu) ([2*x(1); 2*mu*x(2)]);
A=@(mu)[2 0 ; 0 2*mu];


%% Minimize f for mu = 1

mu=1;

% Point (10,0)
x=[10,0];

[~,energy_f,grad_norms,iter]=steepestDescent(A(mu),x,mu,f,f_gradient,N,epsilon);

%Print the plots for Ex. 3.5
printPlots(f,mu,energy_f,grad_norms,iter,epsilon);


% Point (0,10)
x=[0,10];
[~,energy_f,grad_norms,iter]=steepestDescent(A(mu),x,mu,f,f_gradient,N,epsilon);

%Print the plots for Ex. 3.5
printPlots(f,mu,energy_f,grad_norms,iter,epsilon);

% Point (10,10)
x=[10,10];
[~,energy_f,grad_norms,iter]=steepestDescent(A(mu),x,mu,f,f_gradient,N,epsilon);

%Print the plots for Ex. 3.5
printPlots(f,mu,energy_f,grad_norms,iter,epsilon);

%% Minimize f for mu = 10

mu=10;

% Point (10,0)
x=[10,0];

[~,energy_f,grad_norms,iter]=steepestDescent(A(mu),x,mu,f,f_gradient,N,epsilon);

%Print the plots for Ex. 3.5
printPlots(f,mu,energy_f,grad_norms,iter,epsilon);


% Point (0,10)
x=[0,10];
[~,energy_f,grad_norms,iter]=steepestDescent(A(mu),x,mu,f,f_gradient,N,epsilon);

%Print the plots for Ex. 3.5
printPlots(f,mu,energy_f,grad_norms,iter,epsilon);
 
% Point (10,10)
x=[10,10];
[~,energy_f,grad_norms,iter]=steepestDescent(A(mu),x,mu,f,f_gradient,N,epsilon);

%Print the plots for Ex. 3.5
printPlots(f,mu,energy_f,grad_norms,iter,epsilon);












