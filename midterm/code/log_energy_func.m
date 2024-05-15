function log_energy_values = log_energy_func(points,solution,matrix,iter)
%Function that computes the logarithm energy norm of the error 
% points is a squared matrix nxn in which are contained the points values
% solution is the solution computed by the Conjugate gradient
% matrix is the matrix required to compute the logarithm energy norm of the error
% iter is the number of iterations


log_energy_values = zeros(1, iter);
i = 1;
while i<=iter
    % Extract the column vector x 
    x = points(i,:)';
    % Calculate the energy logarithm of the error for the current vector
    log_energy_values(i) = log((x-solution)' * matrix * (x-solution));
    i=i+1;
end

end

