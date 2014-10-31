function [rho] = get_rho_from_TP(P, T, PDT, rho_guess, fsolve_options)

% [rho, Q] = refpropm('DQ','T',T,'P',P/1e3,'N2O');

% if Q < 1 && Q > 0

[rho, ~, exitflag] = fzero(@(rho) rho_eqns(rho, T, P, PDT),[rho_guess],fsolve_options);

if exitflag ~= 1
    disp('something wrong')
end

function F = rho_eqns(rho, T, P, PDT)

% [P_var] = refpropmJEZ('P','T',T,'D&',rho,'N2O');
if rho > 500
    P_var = interp2(PDT.T,PDT.rho,PDT.P,T,rho,'linear');
else
    P_var = interp2(PDT.T,PDT.rho,PDT.P,T,rho,'linear');
end
P_var = P_var*1e3;

F = (P - P_var)/P;