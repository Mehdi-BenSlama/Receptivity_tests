
delS = 0.01;
S = -1+delS:delS:1;

a = 1;

% Real space
X = tan(pi*S/2)

A0_0 = abs(X)

A0 = fsolve( @(A) A0_solver(A,X,a), 0*A0_0);