% generate distribution and moments


clear all

N_mom = 8;

dist = 'gamma';

if ~strcmp(dist,'random') && ~strcmp(dist, 'mcgraw')
    
    alpha_val = 5;
    beta_val = 1;
    
    syms alpha beta t
    
    switch dist
        case 'gamma'
            % gamma dist
            m = (1 - beta*t)^-alpha;
        case 'weibull'
            % weibull dist
            m = beta^(t/alpha) * gamma( 1 + t/alpha);
        case 'normal'
            % normal dist
            % here mu = alpha
            % sigma = beta
            m = exp(alpha*t + (t^2*beta^2/2) );     
    end
    
    for i = 1:N_mom
        deriv = diff(m, t, i);
        moment(i) = subs(deriv, [alpha beta t], [alpha_val beta_val 0]);
    end
    
    x = logspace(-10,3,1000);
    
    switch dist
        case 'gamma'
            
            % gamma dist
            pdf = beta_val^alpha_val * x.^(alpha_val - 1) .* exp( -beta_val*x )...
                /gamma(alpha_val);
        case 'weibull'
            % weibull dist
            pdf = alpha_val / beta_val * x.^(alpha_val - 1) .* exp( -x.^alpha_val / beta_val );
            
        case 'normal'
            % normal dist
            pdf = exp(- (x - alpha_val).^2/(2 * beta_val^2) )/( beta_val*sqrt(2*pi) );
    end
    
    moments_of_pdf = double(moment);
%     
%     hold on
%     plot(x,pdf)
%     set(gca,'xscale','log')
    
else
    
    switch dist
        case 'random'
    
            moments_of_pdf = rand(N_mom,1);
    
        case 'mcgraw'
            moments_of_pdf = [1 5 33.3333 277.778 2777.78 32407.4];
    end
end

[r,w] = PD_method(moments_of_pdf);

if sum(isnan([r(:); w(:)])) > 0
    disp('nans result from PD method')
end

if sum([r(:); w(:)] < 0 ) > 0
    disp('negative abscissas or weights from PD method')
end

if sum(imag([r(:); w(:)])) > 0
    disp('complex abscissas or weights from PD method')
end
save('moments_for_qmom.mat','moments_of_pdf')


