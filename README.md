# iterative-method-for-a-linear-system-of-equations
Numerical methods
Problem 1.
Zaimplementować metody iteracyjne Jacobiego, Gaussa-Seidla i SOR do rozwiązywania układów równań liniowych. Należy zaprezentować wynik oraz metodę weryfikacji poprawności  rozwiązania (implementując reziduum rozwiązania) na dowolnej przykładowej macierzy 3x3 oraz dowolnie dobranego współczynnika relaksacji omega.
Wybrane następujące macierzy:

A=\left|\begin{matrix}-6&-1&3\\4&6&-3\\1&2&8\\\end{matrix}\right|
	

(1)
b=\left[\begin{matrix}1&2&32\\\end{matrix}\right]
	(2)
\left|\begin{matrix}-6&-1&3\\4&6&-3\\1&2&8\\\end{matrix}\right|\bullet x=\left[\begin{matrix}1\\2\\32\\\end{matrix}\right]
	
(3)
 
Współczynnik relaksacji odpowiednio ω = 1.15.
Przebiegi rezidua (rys.1) rozwiązania w kolejnych iteracjach w skali liniowej wyglądają następująco:
 
1.Przebiegi residua rozwiązań w kolejnych iteracjach s wskali liniowej.
W skali logarytmicznej (rys.2):
 
2.Przebiegi residua rozwiązań w kolejnych iteracjach s wskali logarytmicznej.

Listing do problemu 1:
function problem1
n = 3;              		% dimensions of an A-matrix 
 
A = [ -6 -1  3
       4  6 -3
       1  2  8 ];
   
b = [1;2;32];
 
x = zeros(n,1);     	% is an n-by-1 matrix of zeros.
A = A + eye(n);     	% is the n-by-n identity matrix.
 
[~, res_jacobi] = jacobi(A,b,x);
[~, res_gausa_seidla] = gauss_seidel(A,b,x);
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
L = tril(A,-1); 		% Extract lower triangular part.
U = triu(A,1); 	 	% Extract upper triangular part.
for i=1:60
x = -inv(D)*(L+U) * x + D\b;
R(i) = norm(A*x - b)/norm(x);
end
end
 
%% Gauss–Seidel method (Liebman method or the method of successive displacement)
function [x, R] = gauss_seidel(A,b,x)
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




Wnioski.
Jak widać z rysunków, każda z badanych metod pozwoliła na rozwiązanie naszego układu równań. Jednak trzeba zaznaczyć, że najlepszy wynik prezentuję metoda Gausa-Seidela, w przypadku którego residium najszybciej spada do wartości 0. W przypadku metody SOR (Successive over-relaxation) szybkość zbieżnośći residium do zera zależy przede wszystkim od wartości współczynnika relaksacji.
Problem 2.
Sprawdzić zbieżność macierzy załączonych do ćwiczenia w ISODzie za pomocą promienia spektralnego dla metod Jacobiego, Gaussa-Seidla i SOR (w=0.7, w=1, w = 1.3).

Matrix	Method	Spectral radius
	Is convergence?	Number of itterations


Mac_1	Jacobi	0.5445	+	15
	Liebman	0.2570	+	8
	SOR (w= 0.7)	0.5588	+	17
	SOR (w= 1.0)	0.2570	+	8
	SOR (w= 1.3)	0.7026	+	29


Mac_2	Jacobi	0.6036	+	20
	Liebman	0.3925	+	11
	SOR (w= 0.7)	0.6412	+	22
	SOR (w= 1.0)	0.3925	+	11
	SOR (w= 1.3)	0.4048	+	12


Mac_3
Basic = 20	Jacobi	3.3258	-	70
	Liebman	0.7613	+	30
	SOR (w= 0.7)	0.8227	+	45
	SOR (w= 1.0)	0.7613	+	30
	SOR (w= 1.3)	0.8299	+	49


Mac_3
Basic = 1	Jacobi	0.7705	+	39
	Liebman	0.1671	+	5
	SOR (w= 0.7)	0.4230	+	10
	SOR (w= 1.0)	0.1671	+	5
	SOR (w= 1.3)	0.4058	+	11


Mac_3
Basic = 0.1	Jacobi	0.1021	+	5
	Liebman	0.0132	+	3
	SOR (w= 0.7)	0.3160	+	8
	SOR (w= 1.0)	0.0132	+	3
	SOR (w= 1.3)	0.3115	+	9

Mac_4
b = [1 1 2 3]’	Jacobi	0.6667	+	23
	Liebman	0.4444	+	13
	SOR (w= 0.7)	0.6867	+	25
	SOR (w= 1.0)	0.4444	+	13
	SOR (w= 1.3)	0.3000	+	9

Wnioski. 
Analizująć dane zamieszczone w tablice, można dostać wniosku, że liczba iteracji zależy od promienia spektralnego. Liczba iteracji charakteryzuje szybkość zbieżnośći badanych metod. Im bliższy zera jest promień spektralny – tym mniejsza jest liczba iteracji, co w swoją kolej powoduje szybką zbieżnosć.
Spośród badanych zestawów macierzy, najgorze wyniki ma zestaw mac_3: liczba iteracji dla prawie wszystkich metod jest dość duża, co zapewnia zbieżność, ale nie taką szybką jak dla pozostawych zestawów. Jednak metoda Jacobiego dla zestawu mac_3 jest niezbieżna.
Problem 3.
Sprawdzić jak zmienia się zbieżność po zastosowaniu ściskania macierzy (prekondycjonowania). Należy zaimplementować prekondycjonerów Gaussa-Seidela oraz prekondycjonera Jacobiego.
Należy zaimplementować prekondycjoner lewostronny zgodnie z formułą:
 

gdzie:
Dla prekondycjonera Jacobiego:  
Dla prekondycjonera Gaussa-Seidela:  
Matrix	Method	Spectral radius
	Number of itterations
		without preconditioner	Jacobi	Liebman	without preconditioner	Jacobi	Leibman


Mac_1	Jacobi	0.5445	0.4702	0.2247	15	12	8
	Liebman	0.2570 	0.1499 	0.1010 	8 	7 	6 
	SOR (w= 0.7)	0.5588 	0.4848 	0.4154 	17 	14 	13 
	SOR (w= 1.0)	0.2570 	0.1499 	0.1010 	8 	7 	6 
	SOR (w= 1.3)	0.7026 	0.6070 	0.5053 	29 	20 	11 


Mac_2	Jacobi	0.6036 	0.4024 	0.2630 	20 	11 	8 
	Liebman	0.3925 	0.2012 	0.1383 	11 	7 	6 
	SOR (w= 0.7)	0.6412 	0.5040 	0.4382 	22 	14 	13 
	SOR (w= 1.0)	0.3925 	0.2012 	0.1383 	11 	7 	6 
	SOR (w= 1.3)	0.4048 	0.3704 	0.4629 	12 	11 	14 


Mac_3
Basic = 20	Jacobi	3.3259 	2.1543 	0.3136 	- 	- 	1 
	Liebman	0.7613 	0.5278 	0.0557 	30 	15 	1 
	SOR (w= 0.7)	0.8227 	0.5842 	0.3466 	45 	13 	9 
	SOR (w= 1.0)	0.7613 	0.5278 	0.0557 	30 	15 	1 
	SOR (w= 1.3)	0.8299 	0.6590 	0.3361 	49 	24 	9 


Mac_3
Basic = 1	Jacobi	0.7705 	2.1543 	0.3136 	39 	- 	1 
	Liebman	0.1671 	0.5278 	0.0557 	5 	15 	1 
	SOR (w= 0.7)	0.4230 	0.5842 	0.3466 	10 	13 	9 
	SOR (w= 1.0)	0.1671 	0.5278 	0.0557 	5 	15 	1 
	SOR (w= 1.3)	0.4058 	0.6590 	0.3361 	11 	24 	9 


Mac_3
Basic = 0.1	Jacobi	0.1021 	2.1543 	0.3136 	5 	- 	1 
	Liebman	0.0132 	0.5278 	0.0557 	3 	15 	1 
	SOR (w= 0.7)	0.3160 	0.5842 	0.3466 	8 	13 	9 
	SOR (w= 1.0)	0.0132 	0.5278 	0.0557 	3 	15 	1 
	SOR (w= 1.3)	0.3115 	0.6590 	0.3361 	9 	24 	9 

Mac_4
b = [1 1 2 3]’	Jacobi	0.6667 	0.5000 	0.3333 	23 	14 	9 
	Liebman	0.4444 	0.2500 	0.1111 	13 	8 	6 
	SOR (w= 0.7)	0.6867 	0.5625 	0.4579 	25 	17 	14 
	SOR (w= 1.0)	0.4444 	0.2500 	0.1111 	13 	8 	6 
	SOR (w= 1.3)	0.3000 	0.3000 	0.3000 	9 	9 	8 
 
Listing do problem 3.
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
legend('without preconditioner’, ‘Jacobi’s preconditioner’, ‘Liebman’s preconditioner');
title('SOR (linear scale)');
grid on;
hold off
 
end




Implementacji metod Jacobiego, Liebmanna i SOR są identyczne jak w zadaniu 1 i 2. Kombinację tych metod z różnymi rodzajami prekondycjonerów wyglądają następująco:

%% Jacobi method implementation with Jacobi’s preconditioner

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

%% Jacobi method implementation with Liebman’s preconditioner
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

%% Liebman method implementation with Jacobi’s preconditioner
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

%% Liebmann method implementation with Liebman’s preconditioner
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

%% The method of successive over-relaxation with Jacobi’s preconditioner
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



Wnioski.
Analizując poszczególne zestawy macierzy, można dostać wniosku, że dla mac_1 prekondencjonery polepszały znaczenia promienia spektralnego oraz liczby iteracji. Ponadto, prekondycjoner Liebmana w porównaniu do Jacobiego miał mniejsze znaczenia badanych wielkości.
Jeżeli chodzi o mac_2, to otrzymane wyniki są podobne do mac_1, za wyjątkiem tego, że dla metody SOR (w= 1.3) prekondycjoner Liebmana  prowadzi do najgorszego wyniku.
W ogólności wykorzystanie prekondycjonerów spowodowało szybczą zbieżność niż dla układów równań bez nich. Zastosowanie prekondycjonera Liebmana dawało najlepsze rezultaty jak dla metody Jacobiego, tak i dla Liebmana i SOR, i było najbardziej efektywnym rozwiązaniem problemu. Natomiast prekondycjoner Jacobiego, choć przeważnie nie tak żle poprawiał wyniki, czasami nie miał możliwości zapewnić zbieżnośc badanej metody.




