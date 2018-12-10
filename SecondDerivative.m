function [ddA] = nonUniformD(A,X)
    %NONUNIFORMD Summary of this function goes here
    %   Detailed explanation goes here

    
    dXp = X(3:end) - X(2:end-1);
    dXm = -X(2:end-1) + X(1:end-2);
    
    
    ddA = 2*(dXp .* (A(2:end-1) - A(1:end-2)) - dXm .* (A(2:end-1) - A(3:end)))./(dXp.^2.*dXm - dXp.*dXm.^2);
    ddA = [0 ddA 0]
    
end

