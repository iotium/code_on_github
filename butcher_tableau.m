% initialized ODE solver

function [adaptive, a, b, c, bs, s] = butcher_tableau(ode_solver)


% for differential equation solver:
% y(n+1) = y(n) + sum(i = 1 -> s) of b(i)*k(i)
% k(i) = h*f(stuff)
% h = step size
% f = the ode's  (dy/dt)
% f for k(1) = f( t(n), y(n) )
% f for k(2) = f( t(n) + c(2)*h , y(n) + a(2,1)*k(1) )
% f for k(3) = f( t(n) + c(3)*h , y(n) + a(3,1)*k(1) + a(3,2)*k(2) )
% ...
% the bs denotes b_star, and is used for error estimation


switch ode_solver
    case 'RK4'
        % 4th order runge-kutta
        adaptive = 0;
        
        a = [0  0   0   0;
            .5  0   0   0;
            0   .5  0   0;
            0   0   0   1];
        
        c = [0; .5; .5; 1];
        
        b = [1/6; 1/3; 1/3; 1/6];
        
        bs = nan;
        
    case 'euler'
        % 1st order euler
        adaptive = 0;
        
        a = [0];
        
        b = 1;
        
        c = [0];
        
        bs = nan;
        
    case 'RKF'
        adaptive = 1;
        
        % 4th/5th order runge kutta fehlberg
        % a, b, bs (ie b star), and c comprise the butcher tableau
        % c is what's usually shown on the vertical axis
        % b is for the 5th order solution, bs for the 4th
        a = [0          0           0           0           0       0;
            .25         0           0           0           0       0;
            3/32        9/32        0           0           0       0;
            1932/2197   -7200/2197  7296/2197   0           0       0;
            439/216     -8          3680/513    -845/4104   0       0;
            -8/27       2           -3544/2565  1859/4104   -11/40  0];
        
        c = [0;         .25;        3/8;        12/13;      1;      .5];
        
        b = [16/135;    0;          6656/12825; 28561/56430;-9/50;  2/55];
        
        bs = [25/216;   0;          1408/2565;  2197/4104;  -1/5;   0];
        
    case 'CK'
        adaptive = 1;
        
        % cash-karp
        c = [0; 1/5; 3/10; 3/5; 1; 7/8];
        
        b = [37/378      0           250/621     125/594     0               512/1771]';
        
        bs = [2825/27648	0           18575/48384	13525/55296	277/14336       1/4]';
        
        a = [0          0       0           0               0           0;
            1/5        0       0           0               0           0;
            3/40       9/40    0           0               0           0;
            3/10       -9/10   6/5         0               0           0;
            -11/54      5/2     -70/27      35/27           0           0;
            1631/55296	175/512 575/13824	44275/110592	253/4096	0];
        
    case 'DP'
        adaptive = 1;
        
        % dormand-prince
        
        c = [0              1/5     3/10        4/5         8/9         1           1]';
        
        b = [35/384         0       500/1113	125/192	-2187/6784      11/84       0]';
        
        bs = [5179/57600	0       7571/16695	393/640	-92097/339200	187/2100	1/40]';
        
        
        a =[0           0           0           0           0       0       0;
            1/5         0           0           0           0       0       0;
            3/40        9/40        0           0           0       0       0;
            44/45       -56/15      32/9        0           0       0       0;
            19372/6561	-25360/2187	64448/6561	-212/729    0       0       0;
            9017/3168	-355/33     46732/5247	49/176	-5103/18656 0       0;
            35/384      0           500/1113	125/192	-2187/6784	11/84   0];
        
        
        
        
end

s = length(c); % number of stages in the scheme