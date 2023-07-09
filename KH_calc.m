function [Kh] = KH_calc(Ut,B,x);
% KH_calc.m  3/9/2005  Parker MacCready
% 
% This calculates Kh using various parameterization of
% tidally-averaged mixing

Tt = 12.42*3600;    % tidal period (s)
Lt = Ut*Tt/pi;      % tidal excursion (m)
% use the Banas et al (2004) parameterization
% modified to include the case when B > Lt
Kh1 = 0.035 * Ut .* min(B,Lt);

% add a mouth-stirring increase to Kh
Ltm = Lt(end); Utm = Ut(end);
eps1 = 0.1; % an efficiency factor to account for the fact
% that often the mouth region has a lot of plume water in it [0.1]
if 0
    Khm = eps1 * ( Ltm^2 - (sqrt(2*B(end)/pi))*(Ltm^1.5) ) / Tt; % Khm (m2 s-1)
    Kh2 = Khm*exp(x/Ltm);
else
    Khm = eps1 * (Ltm*Utm/pi) * (1 - sqrt(2*B(end)/(pi*Ltm))); % Khm (m2 s-1)
    if Khm<=0; Khm = 0; end
    Kh2 = Khm*(1 + x/Ltm);
end

% combine, using the local maximum
Kh = max(Kh1,Kh2);
