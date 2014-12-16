function moments = gamma_dist_moments(alpha, beta, N)

moments = zeros(N,1);

for i = 1:N
    if i == 1
        moments(i) = alpha*beta;
    else
        moments(i) = moments(i-1) * beta * (alpha + (i-1) );
    end
end
   