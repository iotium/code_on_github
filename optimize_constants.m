% optimize constants to fit data

function optimize_constants

A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
lb = log([1     1   .01   0.001    .001 ]);     
ub = log([10    1e4 1e4     10     10 ]);
x0 = mean([lb; ub]);
% x0 = lb;
options = optimset('display','iter','UseParallel','Always','OutputFcn',@output_fn);
matlabpool
[x, fval, exitflag] = fmincon( @(x) cost_fn(x), x0, A, b, Aeq, beq, ...
    lb, ub, nonlcon, options);

save optimization_results

function E = cost_fn(x)
    
    x = exp(x);
    
    load data_for_finding_constants

try

    [t_sim, P_sim, fill_level_f] = bubble_growth(x(1), x(2), x(3), x(4), x(5), 'intermediate_save_for_optimization_run');


    P_exp = interp1(t/t(end), P, t_sim/t_sim(end));

    E = trapz(t_sim/t_sim(end), abs(P_exp - P_sim))/mean(P);
    
    if fill_level_f > 0.03
        E = E*mean(P);
    end

catch err
    
    E = mean(P);
    
end
    
function stop = output_fn(x, optimValues, state)

stop = false;

save intermediate_optimization_results


