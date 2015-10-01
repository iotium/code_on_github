% iterative Ax=b solver

function [x, error_flag] = linear_equation_solver(A,b)

converged = 0;
k = 1;
switch_val = 0;

[L,U,P] = lu(A);

y = L\(P*b);

x = U\y;

% disp('new call')

delta_old = 1;

while converged == 0;
    
%     dy = L\(P * ( (A*x) - b ) );
%     
%     dx = U\dy;

dx = A\( (A*x - b));
    
    x = x - dx;
    
    k = k + 1;
    
    delta_new = max((abs(x) > 1e-6 ).*abs(dx./(x + eps)));
%     fprintf('delta = %4.4g\n',delta_new);
    
    if delta_new < 1e-6
        converged = 1;
        error_flag = 0;
    end
    
    if k > 50 && switch_val == 0
        % try new x:
        x = A\b;
        switch_val = 1;
    end
    
%     if k > 10 && abs(delta_new - delta_old)/delta_old < 1e-3
%         disp('stuck')
%     end
    
    if k > 100
%         error(['failed to converge, delta = ' num2str(delta_new)])
%     disp('linear equation solver not converging')
    converged = 1;
    error_flag = 1;
    
    end
    
    delta_old = delta_new;
    
end