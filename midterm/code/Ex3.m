%% Exercise n° 3


% 3.1 

% A1
A1 = diag(1 : 10);

%Print A1 eigenvalues
num_eigvals_A1=count_unique_eigvals(A1);

% A2
A2 = diag(ones(1, 10));

%Print A2 eigenvalues
num_eigvals_A2=count_unique_eigvals(A2)


% A3
A3 = diag([1, 1, 1, 3, 4, 5, 5, 5, 10, 10]);

%Print A3 eigenvalues
num_eigvals_A3=count_unique_eigvals(A3)

% A4
A4 = diag([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]);

%Print A4 eigenvalues
num_eigvals_A4=count_unique_eigvals(A4)


%% 3.2 & 3.3

%Set the seed for "rand" to guarantee the same repetition
rng(1);
%Right-hand side column vector
b = rand(10,1);

sizeA1=size(A1,1);
sizeA2=size(A2,1);
sizeA3=size(A3,1);
sizeA4=size(A4,1);

%Parameters

tol = 1e-6;
max_itr=10;

%Conjugate gradient applied to A1
[x_A1,points_A1,iter_A1] = CG(A1,b,ones(1,sizeA1)',max_itr,tol);


%Conjugate gradient applied to A2
[x_A2,points_A2,iter_A2] = CG(A2,b,ones(1,sizeA2)',max_itr,tol);


%Conjugate gradient applied to A3
[x_A3,points_A3,iter_A3] = CG(A3,b,ones(1,sizeA3)',max_itr,tol);


%Conjugate gradient applied to A4
[x_A4,points_A4,iter_A4] = CG(A4,b,ones(1,sizeA4)',max_itr,tol);


%% 3.4

%A1 
printPlots(points_A1,x_A1,A1,iter_A1);

%A2
printPlots(points_A2,x_A2,A2,iter_A2);

log_energy_values= log_energy_func(points_A2,x_A2,A2,iter_A2);

hold on;
plot([0,1],[log_energy_values(1),tol],'-bo'); % print the last iteration on plot

%A3
printPlots(points_A3,x_A3,A3,iter_A3);

%A4
printPlots(points_A4,x_A4,A4,iter_A4);


function num_eigvals = count_unique_eigvals(matrix)
    %Support function use to count the number of unique eigenvalues in a
    %matrix
    [~, eigvals] = eig(matrix);
    matrix_eigs=unique(eigvals)';
    num_eigvals =size(matrix_eigs(2:end),2); 
end
