function [sys] = choose_an_estuary(n_system_in);
% choose_an_estuary.m  7/27/2006  Parker MacCready
%
% this sets up the estuarine system to be to be solved by eta2d.m

if nargin==0;
    n_system_list = [1:.1:100];
else
    n_system_list = n_system_in;
end

for jj = 1:length(n_system_list)
    n_system = n_system_list(jj);
    
    % These apply to all systems, unless overriden in the selection
    L = 400e3;              % channel length (m) [400e3]
    nx = 400;               % number of points to use [400]
    x = linspace(-L,0,nx);  % x-axis (m) with the last point at the mouth
    % note that positive-x points to the mouth, at x=0
    xkm = x/1000;           % x-axis (km) for plotting
    dx = x(2)-x(1);         % grid spacing (m)
    Kv_calc = 1;            % set to 0 to specify Km and Ks
    Kh_calc = 1;            % set to 0 to specify Kh
    Socn = 30;              % oceanic salinity (psu)

    Ut = 1 * ones(1,nx);    % default tidal current amplitude (m s-1)
    Qr = 100;               % default river volume flux (m3 s-1)
    H = 10 * ones(1,nx);    % default channel depth (m)
    B = 1000 * ones(1,nx);  % default channel width (m)

    disp_on = logical(1);
    switch round(n_system*10)/10
        case 1
            labtext = 'Diffusive';
            Ut = 2 * Ut;
            Kh_calc = 0;
            Kh = 200 * ones(1,nx); 
            Kv_calc = 0;
            Km = 100e-4 * ones(1,nx);
            Ks = Km;
        case 2
            labtext = 'Chatwin';
            Kh_calc = 0;
            Kh = 0 * ones(1,nx);    
            Kv_calc = 0;
            Km = 14e-4 * ones(1,nx);
            Ks = Km/3;
        case 3
            labtext = 'HR65';
            Ut = 0.4 * Ut;
            A = H.*B;               
            ubar = Qr./A;        
            Kh_calc = 0;
            Kh = ubar.*(x+75e3);
            Kh(Kh<=0) = 0;
        case 4
            labtext = 'Constant Km, Ks, Kh';
            Kv_calc = 0;
            Km = 10e-4 * ones(1,nx);
            Ks = Km/3;
            Kh_calc = 0;
            Kh = 50 * ones(1,nx);
        case 5
            labtext = 'James River';
            Ut = .5 * Ut;
            H = 7 * ones(1,nx);
            B = 2500 * ones(1,nx);
            Qr = 113;
        case 6
            labtext = 'Double H near mouth';
            H = H.*(1 + exp(x/20e3));
        case 7
            labtext = 'Double U_{T} near mouth';
            Ut = 0.5*Ut.*(1 + exp(x/20e3));
        case 8
            labtext = 'Default';
        case 8.1
            labtext = 'Default (Ut/2)';
            Ut = Ut/2;
        case 8.2
            labtext = 'Default (Ut/8)';
            Ut = Ut/8;
        case 8.3
            labtext = 'Default (Ut/2, Kh=0)';
            Kh = 0*Ut;
            Kh_calc=0;
            Ut = Ut/2;
        case 8.4
            labtext = 'Default (Ut/4), Sill';
            Ut = Ut/4;
            sill_fun = (1 + .4*exp(-((xkm+80)/20).^2));
            H = H ./ sill_fun;
            Ut = Ut .* sill_fun;
        case 9
            labtext = 'Willapa Bay (low flow)';
            L = 100e3;
            x = linspace(-L,0,nx);
            xkm = x/1000;
            dx = x(2)-x(1);
            Ut = 1 * ones(1,nx);
            H = linspace(5,15,nx);
            B = linspace(1e3,4e3,nx);
            Qr = 10;
        case 10
            labtext = 'Willapa Bay (medium flow)';
            L = 100e3;
            x = linspace(-L,0,nx);
            xkm = x/1000;
            dx = x(2)-x(1);
            Ut = 1 * ones(1,nx);
            H = linspace(5,15,nx);
            B = linspace(1e3,4e3,nx); 
            Qr = 100;
        case 11
            labtext = 'Willapa Bay (high flow)';
            L = 100e3;
            x = linspace(-L,0,nx);
            xkm = x/1000;
            dx = x(2)-x(1);
            Ut = 1 * ones(1,nx);
            H = linspace(5,15,nx);
            B = linspace(1e3,4e3,nx);
            Qr = 1e3;
        case 12
            labtext = 'North SF Bay (Low Flow)';
            Ut = 0.85 * Ut;
            % Monismith et al. (2002) JPO give Urms = 0.6 m/s (p.3012)
            xkm_data = [min(xkm) -150, -100:10:0];
            A_data = 1e4 * [1 1 1 1 1 1 3 2 2 5 4 7 9];
            H_data = [10 10 10 10 10 10 10 12 12 10 20 20 30];
            % NOTE: last H should be 100
            H = interp1(xkm_data,H_data,xkm);
            A = interp1(xkm_data,A_data,xkm);
            B = A./H;
            Qr = 100;
        case 12.1
            labtext = 'North SF Bay (high Flow)';
            Ut = 0.85 * Ut;
            % Monismith et al. (2002) JPO give Urms = 0.6 m/s (p.3012)
            xkm_data = [min(xkm) -150, -100:10:0];
            A_data = 1e4 * [1 1 1 1 1 1 3 2 2 5 4 7 9];
            H_data = [10 10 10 10 10 10 10 12 12 10 20 20 30];
            % NOTE: last H should be 100
            H = interp1(xkm_data,H_data,xkm);
            A = interp1(xkm_data,A_data,xkm);
            B = A./H;
            Qr = 1000;
        case 14
            labtext = 'Columbia River (Q_{R} = 4000) Neap';
            xkm_data = [min(xkm) -100 -50 -25 0];
            A_data = 1e3 * [5 5 15 35 35];
            H_data = 15 * ones(size(xkm_data));
            H = interp1(xkm_data,H_data,xkm);
            A = interp1(xkm_data,A_data,xkm);
            B = A./H;
            Qr = 4000;
        case 15
            labtext = 'Columbia River (Q_{R} = 4000) Spring';
            Ut = 2* Ut;
            xkm_data = [min(xkm) -100 -50 -25 0];
            A_data = 1e3 * [5 5 15 35 35];
            H_data = 15 * ones(size(xkm_data));
            H = interp1(xkm_data,H_data,xkm);
            A = interp1(xkm_data,A_data,xkm);
            B = A./H;
            Qr = 4000;
        case 16
            labtext = 'Delaware Estuary (low flow)';
            Ut = 0.8 * Ut;
            H = 14 * ones(1,nx);
            A = 3.2e5 * exp((xkm+20)/25);
            A(A>3.2e5) = 3.2e5;
            A(A<1e4) = 1e4;
            B = A./H;
            Qr = 100;
        case 16.1
            labtext = 'Delaware Estuary (high flow)';
            Ut = 0.8 * Ut;
            H = 14 * ones(1,nx);
            A = 3.2e5 * exp((xkm+20)/25);
            A(A>3.2e5) = 3.2e5;
            A(A<1e4) = 1e4;
            B = A./H;
            Qr = 1000;
        case 17
            labtext = 'Hudson River (Low Flow)';
            Ut = 0.9 * Ut;
            % Ut = (0.67 0.9 1.14 HP) = +/- 26.7%
            H = 14 * ones(1,nx);
            A = 16e3 * ones(1,nx);
            B = A./H;
            Qr = 100;
        case 17.1
            labtext = 'Hudson River (High Flow)';
            Ut = 0.9 * Ut;
            H = 14 * ones(1,nx);
            A = 16e3 * ones(1,nx);
            B = A./H;
            Qr = 1000;
        case 17.2   % for the paper figure
            labtext = 'Hudson River: Q_{R} = 300 m^{3} s^{-1}';
            Ut = 0.9 * Ut;
            H = 14 * ones(1,nx);
            A = 16e3 * ones(1,nx);
            B = A./H;
            Qr = 300;
        otherwise
            disp_on = logical(0);
    end
    if disp_on & nargin==0
        disp(['n_system = ',num2str(n_system),' ',labtext])
    end
end

% ********************************************************************
% These apply to all systems
g = 9.8;                % gravity (m s-2)
beta = 7.7e-4;          % paremeter in equation of state (psu-1)
rho0 = 1000;            % background density (kg m-3)
c = sqrt(g*beta*Socn*H);
A = H.*B;               % channel cross-sectional area (m2)
ubar = Qr./A;           % section-averaged velocity (m s-1)
lowsig = .001; % lowest value of sigma to define the salt intrusion

% Initial calculation of the diffusivities
phi = zeros(size(Ut)); % default initial stratification
if Kh_calc; [Kh] = KH_calc(Ut,B,x); end
if Kv_calc; [Km,Ks] = KV_calc(phi,Ut,c,H); end
% Kv: vertical eddy viscosity (m2 s-1)
% Ks: vertical eddy diffusivity (m2 s-1)
% Kh: along-channel eddy diffusivity (m2 s-1)

% Load the system varibles into the structure "sys"
sys.L = L; sys.nx = nx; sys.x = x; sys.xkm = xkm; sys.dx = dx;
sys.H = H; sys.B = B; sys.Qr = Qr; sys.ubar = ubar; sys.g = g; sys.beta = beta;
sys.Socn = Socn; sys.rho0 = rho0; sys.Ut = Ut; sys.Km = Km; sys.Ks = Ks; sys.Kh = Kh;
sys.Kv_calc = Kv_calc; sys.Kh_calc = Kh_calc; sys.n_system = n_system; sys.labtext = labtext;
sys.A = A; sys.c = c; sys.lowsig = lowsig;
