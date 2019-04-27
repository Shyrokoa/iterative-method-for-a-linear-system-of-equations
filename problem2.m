function problem2
close all;
n = 4;

A =[ 2     0    -1    -1
     0     2    -1    -1
    -1    -1     2     0
    -1    -1     0     2];

b=[1 1 2 3]';

x = zeros(n,1);
A = A + eye(n);

[xj, Rjac] = jacobi(A,b,x);
[xgs, Rleib] = gauss_leib(A,b,x);
[xSOR, RSOR] = SOR(A,b,x);

res_jac  = norm(A*xj - b)/norm(xj);
res_leib = norm(A*xgs - b)/norm(xgs);
res_SOR  = norm(A*xSOR-b)/norm(xSOR);

hold on
plot(log(Rjac), 'o')
plot(log(Rleib), 'o')
plot(log(RSOR), 'o')
legend('Jacobi','Gauss-Seidel','SOR');
title('logarithmic scale');
grid on;
hold off

figure
hold on
plot((Rjac), 'o')
plot((Rleib), 'o')
plot((RSOR), 'o')
legend('Jacobi','Gauss-Seidel','SOR');
title('linear scale');
grid on;
hold off
end

%% Jacobi method implementation
function [x, R] = jacobi(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
B = -inv(D)*(L+U);
sr_j = max(abs(eig(B)))
c = D\b;
for i=1:70
x = B* x + c;
R(i) = norm(A*x - b)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_jac = i;
 it_jac
end

%% Liebmann method implementation
function [x, R] = gauss_leib(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
B = -inv(L+D)*(U);
sr_leib = max(abs(eig(B)))
c = (L+D)\b;
for i=1:70
x = B*x+c;
R(i) = norm(A*x - b)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_leib = i;
 it_leib
end

%% The method of successive over-relaxation
function [x, R] = SOR(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
w = 0.7;
% w=1.3;
% w=1;
B = (D+w*L) \ ((1-w)*D - w*U);
sr_SOR = max(abs(eig(B)))
c = w*((D+w*L)\b);
for i=1:70
x = B*x+c;
R(i) = norm(A*x - b)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
it_SOR = i;
it_SOR
end