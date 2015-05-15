%% must be developed!!!
%% Contribution of macropores to derivative of compartment (1/T)
% SWAP-32 uses dFdhMp and I use dSINKmacr_dh:
%   -remember that the calculation of dFdh multiplies dSINKmacr_dh by dz
%    and you should account for it in the following calculations:
dSINKmacr_dh    = zeros(P.nz,1);
% see:
%   ./MacroRate.for:72:               dFdhMp(ic)      = 0.0
%% 
macr            = zeros(P.nz,1);