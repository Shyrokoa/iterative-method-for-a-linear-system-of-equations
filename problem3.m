function problem3
close all;

n = 4;

A =[ 2     0    -1    -1
     0     2    -1    -1
    -1    -1     2     0
    -1    -1     0     2];

b = [1 45 23 1]';

%% Jacobies preconditioner
Mjac = diag(diag(A));
Ajac = Mjac\A;
bjac = Mjac\b;

%% Liebmann`s preconditioner
Mlib = tril(A,-1) + diag(diag(A));
Alib = Mlib\A;
blib = Mlib\b;

x = zeros(n,1);
A = A + eye(n);
Ajac = Ajac + eye(n);
Alib = Alib + eye(n);

[~, Rj]     = jacobi(A,b,x);
[~, Rj_pj]  = jacobi_pr_jac(Ajac,bjac,x);
[~, Rj_pgs] = jacobi_pr_lieb(Alib,blib,x);

%% Jacobi method in logarithmic scale
figure
subplot(2,3,1)
hold on
plot(log(Rj), 'o')
hold on
plot(log(Rj_pj), 'o')
hold on
plot(log(Rj_pgs), 'o')
legend('without preconditioner','Jacobies preconditioner','Liebmann`s preconditioner');
title('Jacobi method (logarithmic scale)');
grid on;
hold off

%% Jacobi method in linear scale
subplot(2,3,4)
hold on
plot((Rj), 'o')
hold on
plot((Rj_pj), 'o')
hold on
plot((Rj_pgs), 'o')
legend('without preconditioner','Jacobies preconditioner','Liebmann`s preconditioner');
title('Jacobi method (linear scale)');
grid on;
hold off
[xgs, Rgs] = liebman(A,b,x);
[xgs_pj, Rgs_pj]=liebman_pr_jac(Ajac,bjac,x);
[xgs_pgs, Rgs_pgs]=liebman_pr_lieb(Alib,blib,x);

%% Liebmann method in logarithmic scale

subplot(2,3,2)
hold on
plot(log(Rgs), 'o')
hold on
plot(log(Rgs_pj), 'o')
hold on
plot(log(Rgs_pgs), 'o')
legend('without preconditioner','Jacobies preconditioner','Liebmann`s preconditioner');
title('Liebmann method (logarithmic scale)');
grid on;
hold off

%% Liebmann method in linear scale
subplot(2,3,5)
hold on
plot((Rgs), 'o')
hold on
plot((Rgs_pj), 'o')
hold on
plot((Rgs_pgs), 'o')
legend('without preconditioner','Jacobies preconditioner','Liebmann`s preconditioner');
title('Liebmann method (linear scale)');
grid on;
hold off

[xSOR, RSOR] = SOR(A,b,x);
[xSOR_pj, RSOR_pj]=SOR_pr_jac(Ajac,bjac,x);
[xSOR_pgs, RSOR_pgs]=SOR_pr_lieb(Alib,blib,x);

%% SOR in logarithmic scale

subplot(2,3,3)
hold on
plot(log(RSOR), 'o')
hold on
plot(log(RSOR_pj), 'o')
hold on
plot(log(RSOR_pgs), 'o')
legend('without preconditioner','Jacobies preconditioner','Leibmann`s preconditioner');
title('SOR (logarithmic scale)');
grid on;
hold off

%% SOR in linear scale
subplot(2,3,6)
hold on
plot((RSOR), 'o')
hold on
plot((RSOR_pj), 'o')
hold on
plot((RSOR_pgs), 'o')
legend('without preconditioner','Jacobies preconditioner','Leibmann`s preconditioner');
title('SOR (linear scale)');
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
sr_jac = max(abs(eig(B)))               % spectral radius
c = D\b;
for i = 1:70
x = B*x + c;
R(i) = norm(A*x - b)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_jac = i;
 it_jac                                 % number of iterations 
end

%% Jacobi method implementation with Jacobies preconditioner
function [x, R] = jacobi_pr_jac(Ajac,bjac,x)
R = [];
D = diag(diag(Ajac));
L = tril(Ajac,-1);
U = triu(Ajac,1);
B = -inv(D)*(L+U);
sr_jac_pr_jac = max(abs(eig(B)))        % spectral radius
c = D\bjac;
for i=1:70
x = B* x + c;
R(i) = norm(Ajac*x - bjac)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_jac_pr_jac = i;
 it_jac_pr_jac                          % number of iterations
end

%% Jacobi method implementation with Leibman`s preconditioner
function [x, R] = jacobi_pr_lieb(Aleib,bleib,x)
R = [];
D = diag(diag(Aleib));
L = tril(Aleib,-1);
U = triu(Aleib,1);
B = -inv(D)*(L+U);
sr_jac_pr_leib = max(abs(eig(B)))       % spectral radius
c = D\bleib;
for i=1:70
x = B* x + c;
R(i) = norm(Aleib*x - bleib)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_jac_pr_leib = i;
 it_jac_pr_leib                         % number of iterations
end

%% Liebmann method implementation
function [x, R] = liebman(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
B = -inv(L+D)*(U);
sr_leib = max(abs(eig(B)))              % spectral radius
c = (L+D)\b;
for i=1:70
x = B*x+c;
R(i) = norm(A*x - b)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_leib = i;
 it_leib                                % number of iterations
end

%% Liebmann method implementation with Jacobies preconditioner
function [x, R] = liebman_pr_jac(Ajac,bjac,x)
R = [];
D = diag(diag(Ajac));
L = tril(Ajac,-1);
U = triu(Ajac,1);
B = -inv(L+D)*(U);
sr_leib_pr_jac = max(abs(eig(B)))       % spectral radius
c = (L+D)\bjac;
for i=1:70
x = B*x+c;
R(i) = norm(Ajac*x - bjac)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_leib_pr_jac = i;
 it_leib_pr_jac                         % number of iterations
end

%% Liebmann method implementation with Leibman`s preconditioner
function [x, R] = liebman_pr_lieb(Aleib, bleib, x)
R = [];
D = diag(diag(Aleib));
L = tril(Aleib,-1);
U = triu(Aleib,1);
B = -inv(L+D)*(U);
sr_leib_pr_leib = max(abs(eig(B)))      % spectral radius
c = (L+D)\bleib;
for i=1:70
x = B*x+c;
R(i) = norm(Aleib*x - bleib)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
 it_leib_pr_leib = i;
 it_leib_pr_leib                        % number of iterations
end

%% The method of successive over-relaxation
function [x, R] = SOR(A,b,x)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
w = 0.7;                                % the relaxation factor
% w = 1;
% w = 1.3;
B = (D+w*L) \ ((1-w)*D - w*U);
sr_sor = max(abs(eig(B)))               % spectral radius
c = w*((D+w*L)\b);
for i=1:70
x = B*x+c;
R(i) = norm(A*x - b)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
it_sor = i;
it_sor                                  % number of iterations
end

%% The method of successive over-relaxation with Jacobies preconditioner
function [x, R] = SOR_pr_jac(Ajac,bjac,x)
R = [];
D = diag(diag(Ajac));
L = tril(Ajac,-1);
U = triu(Ajac,1);
w = 0.7;                                % the relaxation factor
% w = 1;
% w = 1.3;
B = (D+w*L) \ ((1-w)*D - w*U);
sr_sor_pr_jac = max(abs(eig(B)))        % spectral radius
c = w*((D+w*L)\bjac);
for i=1:70
x = B*x+c;
R(i) = norm(Ajac*x - bjac)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
it_sor_pr_jac = i;
it_sor_pr_jac                           % number of iterations
end

%% The method of successive over-relaxation with Leibman`s preconditioner
function [x, R] = SOR_pr_lieb(Aleib,bleib,x)
R = [];
D = diag(diag(Aleib));
L = tril(Aleib,-1);
U = triu(Aleib,1);
w = 0.7;                                % the relaxation factor
% w = 1;
% w = 1.3;
B = (D+w*L) \ ((1-w)*D - w*U);
rs_sor_pr_leib = max(abs(eig(B)))       % spectral radius
c = w*((D+w*L)\bleib);
for i=1:70
x = B*x+c;
R(i) = norm(Aleib*x - bleib)/norm(x);
 if (R(i) < 1e-4)
 break
 end
end
it_sor_pr_leib = i;
it_sor_pr_leib                          % number of iterations
end