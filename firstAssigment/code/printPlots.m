function printPlots(func,mu,energy_f,grad_norms,iter,tol)
%printPlot function : print all plots for the exercise 3.5
% func is the function
% mu is the value required for represent the function
% energy_f is the matrix iter x 3 that contains the x,y fun(x,y) of the
% function
% grad_norms is a column vector containing the gradient norm at each
% itearation of the steepest algorithm
% iter is the number of iteration computed by the iterative algorithm
% tol is the tollerance value for makes the steepest method converge

% Define the range for the plot
[X, Y] = meshgrid(linspace(-10, 10, 100), linspace(-10, 10, 100));


% Plot the iterations on the energy landscape in 2D
figure;
contourf(X,Y,func(X,Y,mu));
title('Energy landscape 2D');
xlabel('x');
ylabel('y');

%Show the values of vector x at each iteration
hold on;
x_points= energy_f(:,1);
y_points= energy_f(:,2);

plot(x_points, y_points, 'ro');
plot(x_points, y_points, 'r-');

figure;
% Plot the logarithm of the norm of the gradient
subplot(1, 2, 1);
semilogy(grad_norms(1:iter)+tol,"b-");
hold on;
semilogy(grad_norms(1:iter)+tol,"bo");
%Note : I sum "grad_norms(1:iter)+tol" due to the fact that the minimum of
%the function is x=[0 0] and the log function does not accept the 0 value
%inside its domain in order to show the gradient norm at that point I used
%the tollerance value as a lower bound 

title('Logarithm of Gradient Norm');
xlabel('Iterations');
ylabel('log_{10} ||\nabla f||');
grid on;

% Plot the value of the energy function
subplot(1, 2, 2);
plot(1:iter, energy_f(1:iter,3),"b-");
hold on;
plot(1:iter, energy_f(1:iter,3),"bo");
title('Value of Energy Function');
xlabel('Iterations');
ylabel('f(x,y)');
grid on;

end

