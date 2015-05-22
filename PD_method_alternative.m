% alternative to PD method
function [r, w] = PD_method_alternative(varargin)
% finds abscissas and weights to match given moments

if nargin == 2
    moments = varargin{1};
    p = varargin{2};
else
    moments = varargin{1};
    p = 1;
end

N = length(moments)/2;

IC(1:N) = log(logspace(-5,-2,N));
IC(N+1:2*N) = log(0.1*ones(1,N));
% 
% LB = zeros(2*N,1);
% UB = ones(2*N,1);

LB = -100*ones(2*N,1);
UB = 10*ones(2*N,1);

m0 = moments(1);

[x] = fsolve(@(x) error_fn(x), IC,optimset('maxfunevals',1e5,'maxiter',1e4));
% [x, fval] = fmincon(@(x) error_fn(x), IC, [], [], [], [], LB, UB ,[], optimset('maxfunevals',1e4,'maxiter',1e4));


r = exp(x(1:N));
w = m0*exp(x(N+1:end));

[r, ind] = sort(r,'descend');
w = w(ind);



    function E = error_fn(x)
       r = exp(x(1:N));
       w = exp(x(N+1:end));
        
        for i = 1:2*N
            E(i) = ((moments(i)/m0) - sum(r.^((i-1)/p).*w))/(moments(i)/m0);
%             E(i) = (log(moments(i)) - sum(((i-1)/p)*log(r) + log(w)))/log(moments(i));

        end
        
%         E = sum(abs(log(E)));

    end

end