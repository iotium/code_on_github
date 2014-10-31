% looking for new E's to try
clear all
Erat = linspace(-1,1,5);

Perr = [0.1337 0.04584 0.01902 0.04594 0.06303;
        0.1539 0.07284 0.05040 0.1173 0.1603;
        0.03895 0.04824 0.0780 0.09181 0.09824;
        0.05531 0.04606 0.06150 0.07191 0.07838 
        0.04466 0.02222 0.04139 0.05037 0.1428];

for i = 1:size(Perr,1)
    
    [~, minind] = sort(Perr(i,:));
    
    minind = minind(1:3);
    
    cubic_fit = fit(Erat(:), Perr(i,:)', 'poly3');   
    quad_fit = fit(Erat(:), Perr(i,:)', 'poly2');
    quad_fit3 = fit(Erat(minind)', Perr(i,minind)', 'poly2');
    
    coef = coeffvalues(quad_fit3);
    a = coef(1);
    b = coef(2);
    c = coef(3);
    
    x = -b/2/a;
    y = feval(quad_fit3,x);
    
    mins = [Perr(i,minind(1)), y, Perr(i,minind(3))];
    
    [min_err, ind] = min(mins);
    
    switch ind
        case 1
            fprintf('min is at log(Erat) = %8.7g\n', Erat(i,minind(1)))
        case 2
            fprintf('min is at log(Erat) = %8.7g\n', x)
        case 3
            fprintf('min is at log(Erat) = %8.7g\n', Erat(i,minind(3)))
    end
        
    figure(i)
    plot(Erat(:), Perr(i,:)','ks')
    hold on
    plot(cubic_fit,'k')
    plot(quad_fit,'k:')
    plot(quad_fit3,'k--')
    legend('data','cubic','quadratic','quadratic 3pt')
    
%     fprintf('cubic = %8.7g\t quintic = %8.7g\n', cubic_min, quint_min)
    
end

    