% function moments = gamma_dist_moments(alpha, beta, N)
function moments_out = gamma_dist_moments(varargin)
if nargin == 3
    alpha = varargin{1};
    beta = varargin{2};
    N = varargin{3};
    
    moments_out = beta_moments([1:N]);
    
    
elseif nargin == 4
    alpha = varargin{1};
    beta = varargin{2};
    N = varargin{3};
    p = varargin{4};
    
    moments_out = beta_moments([1:N]/p);
    
%     max_N = ceil(N/p);
%     
%     
%     moms = beta_moments(max_N);
%     
%     moments_out = interp1([1:max_N], log(moms), [1:N]/p,'spline','extrap');
%     
%     moments_out = exp(moments_out(:));
    
    
else
    error('wrong number of inputs')
end


    function mom = beta_moments(indices)
        
        mom = beta.^indices .* gamma( alpha + indices)/gamma(alpha);
        
%         mom = zeros(N_mom,1);
%         
%         for i = 1:N_mom
%             %             if i == 1
%             %                 mom(i) = alpha*beta;
%             %             else
%             %                 mom(i) = mom(i-1) * beta * (alpha + (i-1) );
%             %             end
%             
%             mom(i) = beta^i * gamma(alpha + i)/gamma(alpha);
%         end
    end

end


