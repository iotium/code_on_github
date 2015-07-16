% get liquid density given temperature and pressure
function D = get_D_from_TP(T, P, guesses, fluid)

D = fsolve(@(D) density_eqn(D), guesses.rho_l, optimset('display','off'));

    function out = density_eqn(D)
        P_var = 1e3*refpropm('P', 'T', T, 'D', D, fluid);

        out = (P_var - P)/P;

    end

end