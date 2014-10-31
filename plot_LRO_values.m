% plotting P_LRO and T_LRO vs t_LRO for my test setup
clear all

NA = 20;
A_i = 0.4e-7;
A_f = 1e-6;
E = 9e2;
A = logspace(log10(A_i),log10(A_f),NA);

m = 7.5688e8;
b = 2.219e2;

for i = 1:NA
    disp(['A_inj = ' num2str(A(i)) ', starting EQ'])
    [P_EQ(i), T_EQ(i), t_EQ(i)] = EQ(A(i),5);
    disp('starting CP')
    [P_CP(i), T_CP(i), t_CP(i)] = CP(A(i),5);
    disp('starting ZK')
    E = m*A(i) + b;
    [P_ZK(i), T_ZK(i), t_ZK(i)] = ZK(A(i),5,E);
end

figure(1)
hold on
plot(1./t_EQ,P_EQ,'k')
plot(1./t_CP,P_CP,'k:')
plot(1./t_ZK,P_ZK,'k--')

figure(2)
hold on
plot(1./t_EQ,T_EQ,'k')
plot(1./t_CP,T_CP,'k:')
plot(1./t_ZK,T_ZK,'k--')