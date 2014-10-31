% finding CdA to match experiments
% also find E for ZK model

function find_CdAEk

k_norm = 0.5;

% loop through experimental data
for i =  1:6
    
    % loop through the 3 models
    % 1 = EQ
    % 2 = CP
    % 3 = ZK
    
    for j = 1:3
        
        specified_case = i; % the dataset
        
        model = j; % the model
        
        switch specified_case
            case 1
                load('their_data.mat','greg')
                % greg's test
                t_LRO_exp = 4.91;
                A_norm = 9.61e-5;
                E_norm = 1*731;
                exp_data.t = greg.t;
                exp_data.p = greg.p*6894.757;
                exp_data.tLRO = t_LRO_exp;
            case 2
                load('their_data.mat','G_284_03')
                
                % 2003 284 test
                t_LRO_exp = 8.74;
                A_norm = 2.119e-5;
                E_norm = 1*3.2e2;
                exp_data.t = G_284_03.t;
                exp_data.p = G_284_03.p*1e3;
                exp_data.tLRO = t_LRO_exp;
                
                [exp_data.t, it] = unique(exp_data.t);
                exp_data.p = exp_data.p(it);
                
            case 3
                load('their_data.mat','G_284_13')
                
                % 2013 284 ground
                t_LRO_exp = 5.56;
                A_norm = 2.888e-5;
                E_norm = 1*711;
                exp_data.t = G_284_13.t;
                exp_data.p = G_284_13.p_tank*6894.757;
                exp_data.tLRO = t_LRO_exp;
            case 4
                load('their_data.mat','F_284_13')
                
                % 2013 284 flight
                t_LRO_exp = 4.85;
                A_norm = 2.936e-5;
                E_norm = 1*1049;
                exp_data.t = F_284_13.t;
                exp_data.p = F_284_13.p_tank*6894.757;
                exp_data.tLRO = t_LRO_exp;
            case 5
                load('their_data.mat','N2O_2')
                
                % N2O-2
                t_LRO_exp = 19.88;
                A_norm = 1.538e-7;
                E_norm = 1*253;
                exp_data.t = N2O_2.t;
                exp_data.p = N2O_2.p;
                exp_data.tLRO = t_LRO_exp;
            case 6
                load('their_data.mat','N2O_11')
                
                % N2O-11
                t_LRO_exp = 4.858;
                A_norm = 6.977e-7;
                E_norm = 1*1e3;
                exp_data.t = N2O_11.t;
                exp_data.p = N2O_11.p;
                exp_data.tLRO = t_LRO_exp;
        end
        
        
        
        switch j
            

            case 1
                % EQ model
                                
                x = fzero(@(x) EQ_model_error(x, A_norm, specified_case, exp_data), ...
                [0.75 1.25], optimset('display','iter','tolx',1e-5));
            
                A_inj = x*A_norm;

            case 2
                % CP model
                x = fminsearch(@(x) model_error(x, A_norm, k_norm, ...
                    specified_case, exp_data, model), ...
                    [1; 1], optimset('display','iter'));

                A_inj = x(1)*A_norm;
                k = x(2)*k_norm;
                               
           case 3
                 % ZK model

                x = fminsearch(@(x) model_error(x, A_norm, E_norm, ...
                    specified_case, exp_data, model), ...
                    [1; 1], optimset('display','iter'));

                A_inj = x(1)*A_norm;
                E = x(2)*E_norm;
            
                
        end

    end
end

    Fs = [];
    y = [];
    load gong.mat
    sound(y)

end

function t_err = EQ_model_error(x, A_norm, specified_case, exp_data)

A_inj = x*A_norm;

[~, ~, t_LRO] = EQ(A_inj, specified_case);

t_err = (t_LRO - exp_data.tLRO)/exp_data.tLRO;

end

% function that calculates the combined pressure and t_LRO error
function [combined_error] = model_error(x, A_norm, E_norm, specified_case, exp_data, model)

% un-normalize variables
A_inj = A_norm*x(1);
E = E_norm*x(2);

switch model
    case 2
        [~, ~, t_LRO, P, ~, t] = CP(A_inj, specified_case, E);
    case 3
        [~, ~, t_LRO, P, ~, t] = ZK(A_inj, specified_case, E);
end

for i = 1:length(t)
    Pexp(i) = interp1(exp_data.t, exp_data.p, t(i));
end

P_err = trapz(t, abs(P - Pexp))/(exp_data.tLRO*mean(Pexp));

t_err = abs(t_LRO - exp_data.tLRO)/exp_data.tLRO;

disp(['P err = ' num2str(P_err) ', t_err = ' num2str(t_err) ', x(1) = ' num2str(x(1)) ', x(2) = ' num2str(x(2))])

combined_error = P_err + 100*t_err;

end