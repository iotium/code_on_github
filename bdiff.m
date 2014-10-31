
function out = bdiff(in,j,n,t,adaptive)
% n = current step
% j = first step
% in = vector of whatever to diff
% h = time step

k = n - j + 1; % number of points to work with

if adaptive == 0
    
    switch k
        case 1
            out = 0;
        case 2
            % 1st order
            h = t(n) - t(n-1);
            out = (in(n) - in(n-1))/h;
        case 3
            % 2nd order
            h = t(n) - t(n-1);
            out = (3/2*in(n) - 2*in(n-1) + 1/2*in(n-2))/h;
        case 4
            % 3rd order
            h = t(n) - t(n-1);
            out = (11/6*in(n) - 3*in(n-1) + 3/2*in(n-2) - 1/3*in(n-3))/h;
        case 5
            % 4th order
            h = t(n) - t(n-1);
            out = (25/12*in(n) - 4*in(n-1) + 3*in(n-2) - 4/3*in(n-3) +1/4*in(n-4))/h;
        case 6
            % fifth order
            h = t(n) - t(n-1);
            out = (137/60*in(n) - 5*in(n-1) + 5*in(n-2) - 10/3*in(n-3) + 5/4*in(n-4) -1/5*in(n-5))/h;
        otherwise
            % 6th order
            h = t(n) - t(n-1);
            out = (49/20*in(n) - 6*in(n-1) + 15/2*in(n-2) - 20/3*in(n-3) + 15/4*in(n-4) - 6/5*in(n-5) + 1/6*in(n-6))/h;
            
    end
    
else
    switch k
        case 1
            out = 0;
            
        case 2
            % 1st order
            out = (in(n) - in(n-1))/(t(n) - t(n-1) );
            
        case 3
            % 2nd order
            dt_dtau = (3/2*t(n) - 2*t(n-1) + 1/2*t(n-2));
            dy_dtau = (3/2*in(n) - 2*in(n-1) + 1/2*in(n-2));
            out = (1/dt_dtau)*dy_dtau;
            
        case 4
            % 3rd order
            dt_dtau = (11/6*t(n) - 3*t(n-1) + 3/2*t(n-2) - 1/3*t(n-3));
            dy_dtau = (11/6*in(n) - 3*in(n-1) + 3/2*in(n-2) - 1/3*in(n-3));
            out = (1/dt_dtau)*dy_dtau;
            
        case 5            % 4th order
            dt_dtau = (25/12*t(n) - 4*t(n-1) + 3*t(n-2) - 4/3*t(n-3) +1/4*t(n-4));
            dy_dtau = (25/12*in(n) - 4*in(n-1) + 3*in(n-2) - 4/3*in(n-3) +1/4*in(n-4));
            out = (1/dt_dtau)*dy_dtau;
        case 6
            % fifth order
            dt_dtau = (137/60*t(n) - 5*t(n-1) + 5*t(n-2) - 10/3*t(n-3) + 5/4*t(n-4) -1/5*t(n-5));
            dy_dtau = (137/60*in(n) - 5*in(n-1) + 5*in(n-2) - 10/3*in(n-3) + 5/4*in(n-4) -1/5*in(n-5));
            out = (1/dt_dtau)*dy_dtau;
        otherwise
            % 6th order
            dt_dtau = (49/20*t(n) - 6*t(n-1) + 15/2*t(n-2) - 20/3*t(n-3) + 15/4*t(n-4) - 6/5*t(n-5) + 1/6*t(n-6));
            dy_dtau = (49/20*in(n) - 6*in(n-1) + 15/2*in(n-2) - 20/3*in(n-3) + 15/4*in(n-4) - 6/5*in(n-5) + 1/6*in(n-6));
            out = (1/dt_dtau)*dy_dtau;
            
            
    end
    
end