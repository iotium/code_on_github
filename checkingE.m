% checking E results

LD = [8.67; 16.6; 6.7; 6.701; 14; 14.01];
L =[65.05; 62.3; 32; 32.01; 14; 14.01];
r = [3.25; 1.875; 2.375; 2.3751; .5; .501];
FL = [64; 90; 95; 85; 87; 87.01];
Pi = [653; 725; 692.8; 790.8; 593.3; 545.8];
t_LRO = [4.91; 8.74; 5.56; 4.85; 19.88; 4.858];
mi = [20; 8.2; 7.2; 6.25; 0.137; 0.14];
E = [1.3e3; 2.4e2; 5.3e2; 5.9e2; 3.4e2; 7.5e2];
Vtank = pi*r.^2.*L;

[curve, gof] = fit(LD,log(E),'poly1');
c = coeffvalues(curve);
m = c(1);
b = c(2);

fit_LD = linspace(5,20,100);
fit_E = c(1)*fit_LD + c(2);

figure(1)
hold on
plot(LD,E,'ks')
plot(fit_LD,exp(fit_E),'k-')
set(gca,'yscale','log')

xlabel('L/D')
ylabel('E')

legend('Experiment','Fit')