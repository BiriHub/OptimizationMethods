function printPlots(points,solution,matrix,iter)
% points is a squared matrix nxn in which are contained the points values
% solution is the solution computed by the Conjugate gradient
% matrix is the matrix required to compute the logarithm energy norm of the error 
% iter is the number of iterations


%Logarithm energy norm of the error 
log_energy_values= log_energy_func(points, solution,matrix,iter);



%Print the logarithm energy norm of the error at each iteration
figure;
plot(0:iter-1, log_energy_values(1:iter), 'b-');
title('Logarithm energy norm of the error');
xlabel('Number of iterations');
ylabel('Energy norm');

hold on;
plot(0,log_energy_values(1),'bo'); % print the first iteration on plot
plot(iter-1,log_energy_values(iter),'bo'); % print the last iteration on plot

end

