function printPlots(func,energy_val,gradient_norms,iter)
%printPlot function : print all plots the iterative method
% func is the function
% energy_f is the matrix iter x 3 that contains the x,y fun(x,y) of the
% function
% grad_norms is a column vector containing the gradient norm at each
% iteration of the iterative method
% iter is the number of iteration computed by the iterative algorithm

%% Plots 

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
figure;
semilogy(0:iter-1,energy_val(1:iter,3), 'b-', 'MarkerSize', 4);
title('Function values along iterations');
xlabel('number of iterations');
ylabel('f(x)');

hold on;
plot(0,energy_val(1,3),'bo'); % print the first iteration on plot
plot(iter-1,energy_val(iter,3),'bo'); % print the last iteration on plot



end

