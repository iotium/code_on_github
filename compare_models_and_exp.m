% compare P and T traces to experimental data
function compare_models_and_exp

load their_data

specified_case = 4;

plot_string = cell(3,1);
plot_string{1} = 'k-';
plot_string{2} = 'b--';
plot_string{3} = 'r-.';

for j = 1:3
    
    switch specified_case
        case 1
            % greg's test
            t_LRO_exp = 4.91;
            exp_data.t = greg.t;
            exp_data.p = greg.p*6894.757;
            exp_data.tLRO = t_LRO_exp;
            switch j
                case 1
                    A = 8.66e-5;
                case 2
                    A = 1.07e-4;
                case 3
                    A = 9.347e-5;
            end
            E = 1.3e3;
            
        case 2
            % 2003 284 test
            t_LRO_exp = 8.74;
            exp_data.t = G_284_03.t;
            exp_data.p = G_284_03.p*1e3;
            exp_data.tLRO = t_LRO_exp;
            
            [exp_data.t, it] = unique(exp_data.t);
            exp_data.p = exp_data.p(it);
            switch j
                case 1
                    A = 1.97e-5;
                case 2
                    A = 2.36e-5;
                case 3
                    A = 2.277e-5;
            end
            E = 2.4e2;
        case 3
            % 2013 284 ground
            t_LRO_exp = 5.56;
            exp_data.t = G_284_13.t;
            exp_data.p = G_284_13.p_tank*6894.757;
            exp_data.tLRO = t_LRO_exp;
            switch j
                case 1
                    A = 2.80e-5;
                case 2
                    A = 3.48e-5;
                case 3
                    A = 2.936e-5;
            end
            E = 5.3e2;
        case 4
            % 2013 284 flight
            t_LRO_exp = 4.85;
            exp_data.t = F_284_13.t;
            exp_data.p = F_284_13.p_tank*6894.757;
            exp_data.tLRO = t_LRO_exp;
            switch j
                case 1
                    A = 2.55e-5;
                case 2
                    A = 2.838e-5;
                case 3
                    A = 2.637e-5;
            end
            E = 1e3;
        case 5
            % N2O-2
            t_LRO_exp = 19.88;
            exp_data.t = N2O_2.t;
            exp_data.p = N2O_2.p;
            exp_data.T = N2O_2.T;
            exp_data.tLRO = t_LRO_exp;
            switch j
                case 1
                    A = 1.55e-7;
                case 2
                    A = 2.36e-7;
                case 3
                    A = 1.560e-7;
            end
            E = 3.4e2;
        case 6
            % N2O-11
            t_LRO_exp = 4.858;
            exp_data.t = N2O_11.t;
            exp_data.p = N2O_11.p;
            exp_data.T = N2O_11.T;
            exp_data.tLRO = t_LRO_exp;
            switch j
                case 1
                    A = 6.87e-7;
                case 2
                    A = 1.07e-6;
                case 3
                    A = 6.977e-7;
            end
            E = 7.5e2;
            
    end
    
    [P,T,t] = tank_model(A, specified_case, j, E);
    
    Perr(j) = check_Perror(P, t, exp_data);
    
    figure(1)
    hold on
    plot(t,P/1e6,plot_string{j},'linewidth',2)
    
    if isfield(exp_data,'T')
        figure(2)
        hold on
        plot(t,T,plot_string{j},'linewidth',2)
    end
        
    
end

fprintf('case = %3.2g\t Perr for EQ = %8.7g\t for CP = %8.7g\t for ZK = %8.7g\n',...
    specified_case, Perr(1), Perr(2), Perr(3))

figure(1)
hold on
plot(exp_data.t, exp_data.p/1e6, 'k:', 'linewidth', 2)
legend('EQ','CP','ZK','Exp.')
xlabel('Time [s]')
ylabel('Pressure [MPa]')

if isfield(exp_data,'T')
    figure(2)
    hold on
    plot(exp_data.t, exp_data.T, 'k:', 'linewidth', 2)
    legend('EQ','CP','ZK','Exp.')
    xlabel('Time [s]')
    ylabel('Temperature [K]')
end

function Perr = check_Perror(P, t, exp_data)

for i = 1:length(t)
    Pexp(i) = interp1(exp_data.t, exp_data.p, t(i));
end

Perr = trapz(t, abs(P - Pexp))/(exp_data.tLRO*mean(Pexp));


function [P,T,t] = tank_model(A_inj, specified_case, model, E)

switch model
    case 2
        [~, ~, ~, P, T, t] = CP(A_inj,specified_case);
    case 1
        [~, ~, ~, P, T, t] = EQ(A_inj,specified_case);
    case 3
        [~, ~, ~, P, T, t] = ZK(A_inj, specified_case, E);
end