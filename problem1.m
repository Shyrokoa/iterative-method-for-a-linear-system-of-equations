function problem1
n = 3;              % dimensions of an A-matrix 

A = [ -6 -1  3
       4  6 -3
       1  2  8 ];
   
b = [1;2;32];

x = zeros(n,1);     % is an n-by-1 matrix of zeros.
A = A + eye(n);     % is the n-by-n identity matrix.

[~, res_jacobi] = jacobi(A,b,x);
[~, res_gausa_seidla] = liebman(A,b,x);
[~, res_SOR] = SOR(A,b,x);

figure
hold on
plot(log(res_jacobi), 'o')
plot(log(res_gausa_seidla), 'o')
plot(log(res_SOR), 'o')
legend('Jacobi','Gauss-Seidel','SOR');
title('A logarithmic scale');
grid on;
hold off

figure
hold on
plot((res_jacobi), 'o')
plot((res_gausa_seidla), 'o')
plot((res_SOR), 'o')
legend('Jacobi','Gauss-Seidel','SOR');
title('A linear scale');
grid on;
hold off
end

%% Jacobi method implementation
function [x, R] = jacobi(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1); % Extract lower triangular part.
U = triu(A,1);  % Extract upper triangular part.
for i=1:60
x = -inv(D)*(L+U) * x + D\b;
R(i) = norm(A*x - b)/norm(x);
end
end

%% Gauss–Seidel method (Liebmann method or the method of successive displacement)
function [x, R] = liebman(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
for i=1:60
x = -inv(L+D)*(U) * x + (L+D)\b;
R(i) = norm(A*x - b)/norm(x);
end
end

%% The method of successive over-relaxation (SOR)
function [x, R] = SOR(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
w = 1.15;                            % the relaxation factor
B = (D+w*L)\((1-w)*D - w*U);
c = w*((D+w*L)\b);
for i=1:60
x = B*x+c;
R(i) = norm(A*x - b)/norm(x);
end
end