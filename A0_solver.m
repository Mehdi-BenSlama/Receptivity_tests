function [FA_0] = A0_solver(A0,X,a)

    % Left hand side for A0 equation
    LHS = A0.^2 - X.^2 + 2*a;
    
    % Initialize Right Hand Side (RHS)
    RHS = zeros(size(A0));
    RHSconstant = gamma(3/4)/sqrt(2)/gamma(5/4);
    
    % Calculate A0 second derivative
    A0dd = SecondDerivative(A0,X);
    
    % Loop to compute integral
    for i = 1:length(RHS)-2
        RHS(i) = sum(...
            (A0dd(i+1:end-1).*sqrt(X(i+1:end-1)-X(i)) + ...
            A0dd(i+2:end).*sqrt(X(i+2:end)-X(i)))       ...
            .*(X(i+2:end)-X(i+1:end-1))/2 ...
            ) ...
            + 2*(X(i+1)-X(i))^(1/2)*A0dd(i); % Head term
    end
    RHS(end-1) = 2*(X(end)-X(end-1))^(1/2)*A0dd(end-1); % Head term
    RHS(end) = 0;
    RHS = RHSconstant*RHS;
    
    
    FA_0 = LHS - RHS;
end

