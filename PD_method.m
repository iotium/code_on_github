function [r,w] = PD_method(moments)

recover_moments = 0;

moments = moments(:)';

% PD method
% following along Marchisio 2003 and McGraw 1997. Both seem to have
% errors, but notation here matches Marchisio
% clear all

% for 2*Nq moments
Nq = length(moments)/2;

if moments(1) == 0
    
    r = rand(Nq,1);
    w = zeros(Nq,1);
    
else
    
    m_0 = moments(1);
    m = [moments, 0]/m_0;
    
    % 6 moments (0-5 )
    % put a 0 at the end
    % use m_0 = 1 here, then multiply everything at the end
    % m_0 = 1;
    % m = [1, 5, 33.3333, 277.778, 2777.78, 32407.4, 0];
    
    P(:,1) = [1; zeros(2*Nq, 1)];
    
    ivec = [1:(2*Nq + 1)];
    
    % note that it should be m_i-1, but matlab doesn't start at 0
    P(:,2) = (-1).^(ivec-1) .* m(ivec);
    
    for j = 3:(2*Nq + 1)
        for i = 1:(2*Nq + 2 - j)
            P(i,j) = P(1,j-1) * P(i+1,j-2) - P(1,j-2)*P(i+1,j-1);
        end
    end
    
    alpha = zeros(2*Nq,1);
    
    for i = 2:2*Nq
        alpha(i) = P(1,i+1)/( P(1,i) * P(1,i-1) );
    end
    
    ivec = [1:(2*Nq - 1)];
    
    for i = 1:Nq
        a(i) = alpha(2*i) + alpha(2*i - 1);
    end
    
    for i = 1:(Nq-1)
        b(i) = sqrt( alpha(2*i + 1) * alpha(2*i) );
    end
    
    J = diag(a, 0) + diag(b, -1) + diag(b, +1);
    
    [eigvec, eigval] = eig(J);
    
    for i = 1:Nq
        w(i) = m(1) * eigvec(1,i)^2;
        r(i) = eigval(i,i);
    end
    
    % w are the weights,
    % r are the abscissas
    w = m_0 * fliplr(w);
    r = fliplr(r);
    
    %     if sum(abs(imag([r(:); w(:)]))) > 0
    %         disp('imaginary abscissas or weights')
    %     end
    if recover_moments
        % % recover moments (for debugging)
        
        for i = 1:(2*Nq)
            rec_moment(i) = sum( r.^(i-1) .* w );
        end
        fprintf(' recovered moments: ')
        fprintf('%6.6g\t',rec_moment)
        fprintf('\n')
    end
end

