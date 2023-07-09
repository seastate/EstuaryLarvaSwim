function [Km,Ks] = KV_calc(phi,Ut,c,H,flag);
% KV_calc.m  9/21/2005  Parker MacCready
% 
% This calculates Km and Ks using various parameterization of
% tidally-averaged mixing
%
% "flag" is used only in Peters_Ks.m where we wish to plot the curve of Ks
% versus RiL for the case Alpha3 = 0.
% ...use flag = 1 to do this...

if nargin == 4
    flag = 0;
end

% set parameters
% note: 0 <= Alpha3 <= 1
Alpha0 = 0.5 * 0.065;     % [0.065/2]
Alpha1 = 0.065/3;    % [.065/3]
Alpha2 = 3.33;        % [10/3]
if flag == 1
    Alpha3 = 0;     % [0]
else
    Alpha3 = 0.3;   % [0.3]
end
Cd = 2.6e-3;
% some calculations
c2 = c.^2;
Ut2 = Ut.^2;
Ks0 = Alpha1*Cd*Ut.*H;
RiL = c2.*phi./Ut2;
% calculate vertical diffusivities
Km = Alpha0*Cd*Ut.*H;
% note that if phi=0 then RiL=0 and Ks=Ks0
Ks = Ks0.*(Alpha3 + (1-Alpha3)./(1+Alpha2*RiL));
