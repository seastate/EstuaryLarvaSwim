
%   These lines of code calculate arrays for the current timestep from the
%   unsteady data arrays in Parker MacCready's eta2d model, unless the
%   current data are outdated.


%   In this version, elements from Parker's plot_tidal.m code are included
%   to approximate tidal fluctuations.


%   Since computing tidal fluctuations, some of these calculations have to
%   be redone every timestep

% while t > t_next_time
if t > t_next_time

    t_index = t_index + 1;
    t_next_time = Tsave_vec_exact(t_index);
    
    ii = t_index;

    ubar = ubar_mat(ii,:); ue = ue_mat(ii,:);
    sigma = Sigma(ii,:); sigmax = Sigmax(ii,:);
    Ks = Ks_mat(ii,:);
    sys = choose_an_estuary(n_system);
    %
    L= sys.L; nx = sys.nx; x = sys.x; xkm = sys.xkm; dx = sys.dx;
    H = sys.H; B = sys.B; Qr = sys.Qr; ubar = sys.ubar; g = sys.g; beta = sys.beta;
    Socn = sys.Socn; rho0 = sys.rho0; Ut = sys.Ut; Km = sys.Km; Ks = sys.Ks; Kh = sys.Kh;
    Kv_calc = sys.Kv_calc; Kh_calc = sys.Kh_calc; n_system = sys.n_system; labtext = sys.labtext;
    A = sys.A; c = sys.c; lowsig = sys.lowsig;
    clear sys

    ubar = ubar_mat(ii,:); ue = ue_mat(ii,:);
    sigma = Sigma(ii,:); sigmax = Sigmax(ii,:);
    Ks = Ks_mat(ii,:);

    
    %   New lines from plot_tidal.m
    
    % derive bulk quantities
    ihi = find(sigma>=lowsig);
    ihi = ihi(1);
    if ihi > 1;
        ilo = ihi-1;
        if FreeMat_flag==1
            xlo = interplin1(sigma(ilo:ihi),x(ilo:ihi),lowsig);
        else
            xlo = interp1(sigma(ilo:ihi),x(ilo:ihi),lowsig);
        end
    else
        xlo = x(1);
    end
    x_new = linspace(xlo,x(end),100);
    if FreeMat_flag==1
        sigmax_new = interplin1(x,sigmax,x_new);
    else
        sigmax_new = interp1(x,sigmax,x_new);
    end
%     phi_new = interplin1(x,phi,x_new);
    sigmax_mean = mean(sigmax_new);
%     phi_mean = mean(phi_new);
    sigma_int = dx*trapz(sigma);

    % try smoothing things out in the tail
    sigmin = eps;
    tail_list = find(sigma<sigmin);
    if ~isempty(tail_list);
        ntail = tail_list(end);
        sigma(1:ntail) = 0;
        sigmax(1:ntail) = 0;
    end

    % set up to plot results
    % make the lower limit for plotting
    nlow = max(find(sigma<=lowsig));
    if isempty(nlow); nlow = 1; end
    if nlow > 10
        xkm_low = xkm(nlow-10);
    else
        xkm_low = min(xkm);
    end
    xkm_hi = max(xkm);

    % calculate the diffusive fraction of up-estuary salt flux
    warning off MATLAB:divideByZero
    nu = Kh.*sigmax ./ (ubar.*sigma);
    warning on MATLAB:divideByZero
    nu(1:nlow) = NaN;
    xx = linspace(xkm(nlow+1),0,20); % x-axis with fewer points
    if FreeMat_flag==1
        nu = interplin1(xkm(nlow+1:end),nu(nlow+1:end),xx);
    else
        nu = interp1(xkm(nlow+1:end),nu(nlow+1:end),xx);
    end
    

    % Calculate the streamfunction in x-zeta (psi)
    nz = 100;
    zeta = linspace(-1,0,nz)';   % dimensionless vertical coordinate
    % make the polynomial expressions  on the meshgrid
    [X,Z] = meshgrid(xkm,zeta);
    Z2 = Z.^2; Z3 = Z.^3; Z4 = Z.^4; Z5 = Z.^5;
    P1 = 0.5 - 1.5*Z2;
    P2 = 1 - 9*Z2 - 8*Z3;
    P3 = -(7/120) + 0.25*Z2 - 0.125*Z4;
    P4 = -(1/12) + 0.5*Z2 -0.75*Z4 - 0.4*Z5;
    UBAR = ones(nz,1)*ubar;
    UE = ones(nz,1)*ue;
    UP = UBAR.*P1 + UE.*P2;
%     U = UBAR + UP;


    %   New lines from plot_tidal.m
    
    UT = ones(nz,1) * Ut;

    % calculate the salinity field
    KS = ones(nz,1)*Ks;
    KH = ones(nz,1)*Kh;
%     KV = ones(nz,1)*Kv;
    H2 = ones(nz,1)*(H.*H);
    SBAR = ones(nz,1)*sigma*Socn;
    SBARX = ones(nz,1)*sigmax*Socn;
    SP = H2.*SBARX.*(UBAR.*P3 + UE.*P4)./KS;
%     S = SBAR + SP;

    % make a matrix of the actual depth
    for iii = 1:nx
        ZZ(:,iii) = H(iii)*zeta;
    end

    
    
    if FreeMat_flag==1
        this_qr = interplin1(td,Qr_vec,tsd(ii));
        this_ut = interplin1(td,Ut0_vec,tsd(ii));
        this_so = interplin1(td,Socn_vec,tsd(ii));

    else
        this_qr = interp1(td,Qr_vec,tsd(ii));
        this_ut = interp1(td,Ut0_vec,tsd(ii));
        this_so = interp1(td,Socn_vec,tsd(ii));
    end

    
end
    
UTIDE = 1.5 * UT .* (1-Z2) * cos(omega*t);
U = UBAR + UP + UTIDE;


u_surface = U(end,:);   % surface velocity (m s-1)
%

%   New lines from plot_tidal.m
STIDE = -1.5 * (1/omega) * SBARX .* UT .* (1-Z2) * sin(omega*t);
S = SBAR + SP + STIDE;

S(S<0) = 0;

% end





