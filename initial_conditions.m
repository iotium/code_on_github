% function returns initial conditions for a given "specified case"
% a number is passed that corresponds to a specific set of experimental data

function [Ti, fill_level, V_tank, L_tank, ...
    A_inj, Cd, Po, T_air, rho_w, cv_w, t_w, D, k_w, fluid] = initial_conditions(specified_case)

switch specified_case
    case 1
        % greg's data
                fluid = 'N2O';

        Ti = 288;           % [K] initial temperature
        fill_level = 0.64;  % [] initial fill_level ratio (by volume)
        %         E = 1.3e3;          % [] heat transfer multiplier
        V_tank = 0.0354;    % [m^3] tank volume
        L_tank = 65.05*0.0254;       % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 9e-5;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 298;        % [K] air temperature
        rho_w = 2700;       % [kg/m^3] density of wall material (aluminum)
        cv_w = 904;         % [J/kg.K] specific heat of wall (aluminum)
        t_w = 0.0254*1/16;  % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 167;          % [W/m.K] thermal conductivity of wall
        
    case 2
        % 2003 284 data
                fluid = 'N2O';

        Ti = 292.5;           % [K] initial temperature
        fill_level = 0.9;        % [] initial fill_level ratio (by volume)
        %         E = 2e2;          % [] heat transfer multiplier
        V_tank = 0.011276;   % [m^3] tank volume
        L_tank = 62.3*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 2700;       % [kg/m^3] density of wall material (aluminum)
        cv_w = 904;         % [J/kg.K] specific heat of wall (aluminum)
        t_w = 0.0254*1/8;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 167;          % [W/m.K] thermal conductivity of wall
    case 3
        % 2013 284 ground
                fluid = 'N2O';

        Ti = 290.7;           % [K] initial temperature
        fill_level = 0.95;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 0.00929;   % [m^3] tank volume
        L_tank = 32*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 2700;       % [kg/m^3] density of wall material (aluminum)
        cv_w = 904;         % [J/kg.K] specific heat of wall (aluminum)
        t_w = 0.0254*1/8;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 167;          % [W/m.K] thermal conductivity of wall
        
    case 4
        % 2013 284 flight
                fluid = 'N2O';

        Ti = 296.5;           % [K] initial temperature
        fill_level = 0.85;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 0.00929;   % [m^3] tank volume
        L_tank = 32*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 2700;       % [kg/m^3] density of wall material (aluminum)
        cv_w = 904;         % [J/kg.K] specific heat of wall (aluminum)
        t_w = 0.0254*1/8;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 167;          % [W/m.K] thermal conductivity of wall
    case 5
        % N2O test 2 from my data
                fluid = 'N2O';

        Ti = 283.7;           % [K] initial temperature
        fill_level = 0.87;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.80e-4;   % [m^3] tank volume
        L_tank = 0.356;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 1360;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 1250;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*1/4;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 0.195;          % [W/m.K] thermal conductivity of wall
        
    case 6
        % N2O test 11 from my data
                fluid = 'N2O';

        Ti = 280.4;           % [K] initial temperature
        fill_level = 0.87;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.80e-4;   % [m^3] tank volume
        L_tank = 0.356;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 1360;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 1250;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*1/4;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 0.195;          % [W/m.K] thermal conductivity of wall
    case 7
        % CO2 test #24 with quartz tube
                fluid = 'CO2';

        Ti = 290.5;           % [K] initial temperature
        fill_level = 0.80;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.8083e-4;   % [m^3] tank volume
        L_tank = 14.05*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 1360;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 1250;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*1/4;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 0.195;          % [W/m.K] thermal conductivity of wall

        d_inj = 0.022 * 0.0254;
        A_inj = pi/4 * (d_inj^2);
        Cd = 0.6545;

    case 8
        % CO2 test #48 with quartz tube
                fluid = 'CO2';

        Ti = 16.2+273.15;           % [K] initial temperature
        fill_level = 0.83;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.8083e-4;   % [m^3] tank volume
        L_tank = 14.05*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 1360;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 1250;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*1/4;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 0.195;

        d_inj = 0.022 * 0.0254;
        A_inj = pi/4 * (d_inj^2);
        Cd = 0.6545;

    case 9
        % CO2 test #236 with glass gauge
                fluid = 'CO2';

        Ti = 18.5+273.15;           % [K] initial temperature
        fill_level = 0.80;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.233e-3;   % [m^3] tank volume
        L_tank = 25.25*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 800;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 500;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*0.551;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 16.3;

        d_inj = 0.060 * 0.0254;
        A_inj = pi/4 * (d_inj^2);
        Cd = 0.6545;
    case 10
        % CO2 test #291 with glass gauge
                fluid = 'CO2';

        Ti = 20.5+273.15;           % [K] initial temperature
        fill_level = 0.82;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.233e-3;   % [m^3] tank volume
        L_tank = 25.25*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 800;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 500;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*0.551;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 16.3;

        d_inj = 0.060 * 0.0254;
        A_inj = pi/4 * (d_inj^2);
        Cd = 0.6545;

    case 11
        % N2O test #289 with glass gauge
                fluid = 'N2O';
        Ti = 19+273.15;           % [K] initial temperature
        fill_level = 0.82;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.233e-3;   % [m^3] tank volume
        L_tank = 25.25*0.0254;     % [m] tank length
        Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 800;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 500;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*0.551;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 16.3;

        d_inj = 0.089 * 0.0254;
        A_inj = pi/4 * (d_inj^2);
        Cd = 0.6545;

    case 12
        % CO2 test #257 with glass gauge (t_LRO = 5.96)
        fluid = 'CO2';
        Ti = 17.6+273.15;           % [K] initial temperature
        fill_level = 0.825;        % [] initial fill_level ratio (by volume)
        %         E = 2.1e4;          % [] heat transfer multiplier
        V_tank = 1.233e-3;   % [m^3] tank volume
        L_tank = 25.25*0.0254;     % [m] tank length
%         Cd = 1;         % [] injector Cd
        %         A_inj = 3.15e-7;       % [m^2] injector area
        Po = 1e5;           % [Pa] external pressure
        T_air = 293;        % [K] air temperature
        rho_w = 800;       % [kg/m^3] density of wall material (polycarb)
        cv_w = 500;        % [J/kg.K] specific heat of wall (polycarb)
        t_w = 0.0254*0.551;   % [m] wall thickness
        D = sqrt(4/pi*V_tank/L_tank);
        % [m] tank diameter
        k_w = 16.3;
        
        d_inj = 0.089 * 0.0254;
        A_inj = pi/4 * (d_inj^2);
%         Cd = 0.6545;
Cd = 0.80;
                
end