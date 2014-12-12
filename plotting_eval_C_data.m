% plotting eval C data
load eval_C_data

figure(1)
subplot(2,2,1)
contourf(C_nuc_rate, C_rdot, t_min)
set(gca,'xscale','log')
set(gca,'yscale','log')
title('t min')
colorbar

% figure(2)
subplot(2,2,2)

contourf(C_nuc_rate, C_rdot, t_peak)
set(gca,'xscale','log')
set(gca,'yscale','log')
title('t peak')
colorbar

% figure(3)
subplot(2,2,3)

contourf(C_nuc_rate, C_rdot, P_min)
set(gca,'xscale','log')
set(gca,'yscale','log')
title('P min')
colorbar

% figure(4)
subplot(2,2,4)

contourf(C_nuc_rate, C_rdot, P_peak)
set(gca,'xscale','log')
set(gca,'yscale','log')
title('P peak')
colorbar