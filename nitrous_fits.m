% curve fits to saturation properties

function varargout = nitrous_fits(varargin)

T = varargin{1};

% values from IHS:
Tc = 309.57;
Pc = 7251e3;
rho_c = 452;

for i = 2:nargin
    property = varargin{i};
    
    switch property
        case 'P'
            b1 = -6.71893;
            b2 = 1.35966;
            b3 = -1.3779;
            b4 = -4.051;
            v = 1 - Tr;
            P = Pc*exp( 1/Tr*(b1*v + b2*v^1.5 + b3*v^2.5 + b4*v^5 ) );
            varargout{i} = P;
        case 'rho_liq'
            b1 = 1.72328;
            b2 = -0.83950;
            b3 = 0.51060;
            b4 = -0.10412;
            v = 1 - Tr;
            rho_liq = rho_c*exp( b1*v^(1/3) + b2*v^(2/3) + b4*v^(4/3) );
            varargout{i} = rho_liq;
        case 'rho_vap'
            b1 = -1.009;
            b2 = -6.28792;
            b3 = 7.50332;
            b4 = -7.90463;
            b5 = 0.629427;
            v = 1/Tr - 1;
            rho_vap = rho_c*exp( b1*v^(1/3) + b2*v^(2/3) + b4*v^(4/3) + ...
                b5*v^(5/3) );
            varargout{i} = rho_vap;
        case 'h_liq'
            b1 = -200;
            b2 = 116.043;
            b3 = -917.225;
            b4 = 794.779;
            b5 = -589.587;
            v = 1 - Tr;
            h_liq = 1e3*(b1 + b2*v^(1/3) + b3*v^(2/3) + b4*v + b5*v^(4/3));
            varargout{i} = h_liq;
        case 'Dh_LV'
            b1 = -200;
            b2 = 116.043;
            b3 = -917.225;
            b4 = 794.779;
            b5 = -589.587;
            v = 1 - Tr;
            h_liq = 1e3*(b1 + b2*v^(1/3) + b3*v^(2/3) + b4*v + b5*v^(4/3));
           b1 = -200;
            b2 = 440.055;
            b3 = -459.701;
            b4 = 434.081;
            b5 = -485.338;
            v = 1 - Tr;
            h_vap = 1e3*(b1 + b2*v^(1/3) + b3*v^(2/3) + b4*v + b5*v^(4/3));
            
            Dh_LV = h_vap - h_liq;
            varargout{i} = Dh_LV;
        case 'h_vap'
            b1 = -200;
            b2 = 440.055;
            b3 = -459.701;
            b4 = 434.081;
            b5 = -485.338;
            v = 1 - Tr;
            h_vap = 1e3*(b1 + b2*v^(1/3) + b3*v^(2/3) + b4*v + b5*v^(4/3));
            varargout{i} = h_vap;
        case 'cp_liq'
            b1 = 2.49973;
            b2 = 0.023454;
            b3 = -3.80136;
            b4 = 13.0945;
            b5 = -14.5180;
            v = 1 - Tr;
            cp_liq = 1e3*b1*( 1 + b2/v + b3*v + b4*v^2 + b5*v^3 );
            varargout{i} = cp_liq;
        case 'cp_vap'
            b1 = 132.632;
            b2 = 0.052187;
            b3 = -0.364923;
            b4 = -1.20233;
            b5 = 0.536141;
            v = 1 - Tr;
            cp_vap = 1e3*b1*(1 + b2*v^(-2/3) + b3*v^(-1/3) + b4*v^(1/3) + b5*v^(2/3) );
            varargout{i} = cp_vap;
        case 'mu_liq'
            b1 = 1.6089;
            b2 = 2.0439;
            b3 = 5.24;
            b4 = 0.0293423;
            q = (Tc - b3)/(T - b3);
            mu_liq = 1e-3*b4*exp(b1*(q-1)^(1/3) + b2*(q-1)^(4/3));
            varargout{i} = mu_liq;
        case 'mu_vap'
            b1 = 3.3281;
            b2 = -1.18237;
            b3 = -0.055155;
            v = (1/Tr - 1);
            mu_vap = 1e-6*exp(b1 + b2*v^(1/3) + b3*v^(4/3));
            varargout{i} = mu_vap;
        case 'k_liq'
            b1 = 72.35;
            b2 = 1.5;
            b3 = -3.5;
            b4 = 4.5;
            v = (1 - Tr);
            k_liq = 1e-3*b1*(1 + b2*v^(1/3) + b3*v^(2/3) + b4*v);
            varargout{i} = k_liq;
        case 'k_vap'
            b1 = -7.0887;
            b2 = -0.276962;
            b3 = 2.88672;
            b4 = 16.6116;
            b5 = -11.8221;
            v = 1 - Tr;
            k_vap = 1e-3*exp(b1 + b2*v^(-2/3) + b3*v^(-1/3) + b4*v^(1/3) + ...
                b5*v^(2/3) );
            varargout{i} = k_vap;
        case 'ST'
            b1 = 69.31;
            b2 = 1.19346;
            b3 = 0;
            v = 1 - Tr;
            ST = 1e-3*b1*v^b2 * (1 + b3*v);
            varargout{i} = ST;
        case 'cp_IDvap'
            b1 = -0.169903;
            b2 = 0.099053;
            b3 = 1.20822;
            b4 = -0.248324;
            cp_IDvap = 1e3*(b1 + b2*Tr^(-1/2) + b3*Tr^(1/2) + b4*Tr);
            varargout{i} = cp_IDvap;
        case 'h_IDvap'
            b1 = -209.559;
            b2 = 61.3277;
            b3 = -52.5969;
            b4 = 249.352;
            b5 = -38.4368;
            h_IDvap = 1e3*(b1 + b2*Tr^(1/2) + b3*Tr + b4*Tr^(3/2) + b5*Tr^2);
            varargout{i} = h_IDvap;
    end
end





