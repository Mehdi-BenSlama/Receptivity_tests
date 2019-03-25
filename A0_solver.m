function [FA_0] = A0_solver(A0,a)

    N = numel(A0);
    delS = 2/(N+1);
    S = -1+delS:delS:1-delS;
    X = tan(pi*s/2);
    
    % Left hand side for A0 equation
    LHS = A0.^2 - X.^2 + 2*a;
    
    % Initialize Right Hand Side (RHS)
    RHS = zeros(N,1);
    RHSconstant = gamma(3/4)/sqrt(2)/gamma(5/4);
    
    % Calculate A0 second derivative
    A0dd = SecondDerivative(A0,X);
    
    % Loop to compute integral
    for i = 1:N-1
        
        fi = A0dd * pi/2 ./ (cos(pi * X/2).^2) .* ((S(i) - S)./(X(i) - X)).^0.25
        
        
    end
    RHS(end-1) = 2*(X(end)-X(end-1))^(1/2)*A0dd(end-1); % Head term
    RHS(end) = 0;
    RHS = RHSconstant*RHS;
    
    
    FA_0 = LHS - RHS;
end

