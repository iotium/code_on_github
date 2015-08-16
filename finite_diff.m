function out = finite_diff(in, first_step, last_step, current_step, h, direction, max_order)
% equispaced grids only!
% n = current step
% j = first step
% in = vector of whatever to diff
% h = time step

% max_order = 3;

n = current_step;

switch direction
    case 'central'
        
        k1 = current_step - first_step + 1; % number of points to work with
        k2 = last_step - current_step; % number of points to work with
        k = min([k1 k2]);
        k = min([(max_order + 1) k]);
        
        switch k
            case 1
                out = 0;
            case 2
                % 2nd order
                out = (0.5*in(n+1) - 0.5*in(n-1))/h;
            case 3
                % 4th order
                out = (-1/12*in(n+2) + 2/3*in(n+1) - 2/3*in(n-1) + 1/12*in(n-2))/h;
            case 4
                % 6th order
                out = (1/60*in(n+3) -3/20*in(n+2) + 3/4*in(n+1) - 3/4*in(n-1) + 3/20*in(n-2) - 1/60*in(n-3))/h;
            otherwise
                % 8th order
                out = (-1/280*in(n+4) + 4/105*in(n+3) -1/5*in(n+2) + 4/5*in(n+1) - 4/5*in(n-1) + 1/5*in(n-2) - 4/105*in(n-3) +1/280*in(n-4))/h;
        end
        
    case 'backwards'
        
        k = current_step - first_step + 1; % number of points to work with

        k = min([(max_order + 1) k]);
    
        switch k
            case 1
                out = 0;
            case 2
                % 1st order
                out = (in(n) - in(n-1))/h;
            case 3
                % 2nd order
                out = (3/2*in(n) - 2*in(n-1) + 1/2*in(n-2))/h;
            case 4
                % 3rd order
                out = (11/6*in(n) - 3*in(n-1) + 3/2*in(n-2) - 1/3*in(n-3))/h;
            case 5
                % 4th order
                out = (25/12*in(n) - 4*in(n-1) + 3*in(n-2) - 4/3*in(n-3) +1/4*in(n-4))/h;
            case 6
                % fifth order
                out = (137/60*in(n) - 5*in(n-1) + 5*in(n-2) - 10/3*in(n-3) + 5/4*in(n-4) -1/5*in(n-5))/h;
            otherwise
                % 6th order
                out = (49/20*in(n) - 6*in(n-1) + 15/2*in(n-2) - 20/3*in(n-3) + 15/4*in(n-4) - 6/5*in(n-5) + 1/6*in(n-6))/h;

        end
    case 'forwards'
        
        k = last_step - current_step; % number of points to work with

        k = min([(max_order + 1) k]);
        
        switch k
            case 1
                out = 0;
            case 2
                % 1st order
                out = (in(n+1) - in(n))/h;
            case 3
                % 2nd order
                out = -(3/2*in(n) - 2*in(n+1) + 1/2*in(n+2))/h;
            case 4
                % 3rd order
                out = -(11/6*in(n) - 3*in(n+1) + 3/2*in(n+2) - 1/3*in(n+3))/h;
            case 5
                % 4th order
                out = -(25/12*in(n) - 4*in(n+1) + 3*in(n+2) - 4/3*in(n+3) +1/4*in(n+4))/h;
            case 6
                % fifth order
                out = -(137/60*in(n) - 5*in(n+1) + 5*in(n+2) - 10/3*in(n+3) + 5/4*in(n+4) -1/5*in(n+5))/h;
            otherwise
                % 6th order
                out = -(49/20*in(n) - 6*in(n+1) + 15/2*in(n+2) - 20/3*in(n+3) + 15/4*in(n+4) - 6/5*in(n+5) + 1/6*in(n+6))/h;

        end
        
end
    