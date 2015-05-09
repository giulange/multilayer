function [dKi_dhi,Ci,teta_pm1] = multilayer_dKi_dhi(h, Psh, i, Ksat)
% [dKi_dhi,Ci,teta_pm1] = multilayer_dKi_dhi(h, Psh, i, Ksat)
% 
% DESCRIPTION
%   This function computes the conductivity derivative to the pressure head
%   (dKi_dhi) following the definition given in SWAP-32, Appendix 5, page
%   239.
%   It also returns the capillary Capacity (Ci) and the Moisture fraction
%   if required.
%   The h passed in input is the pressure head at all nodes and at
%   iteration-step "p" before/after the solution of the tri-diagonal
%   system.
%   At first iteration-step you have to provide the correct pressure head h
%   (??which should be the pressure head from previous time-step).
%   It is similar to the "dhconduc" function in functions.for, SWAP-32
%   fortran code.
% 
% INPUTs
%   h:          Pressure head during the Newton-Raphson iteration scheme.
%   Psh:        Soil hydraulic characteristics at every node of the soil
%               grid.
%   i:          One or more nodes at which outputs are computed.
%   Ksat:       Conductivity at saturation (I should delete this term as it
%               is already provided in Psh).
% 
% OUTPUT
%   dKi_dhi:    Derivative of hydraulic conductivity to the pressure head
%               at node(s) i.
%   Ci:         Capillary Capacity at node(s) i.
%   teta_pm1:   Moisture fraction at node(s) i.

%% main

% Moisture at "p" iteration:
teta_pm1= fnteta( h, Psh, i );
% Relative saturation (SWAP-32, Eq. 2.7, page 28):
% check that:
%   -P.h1(=h_pm1) is what you need to define starting condition!
%   -all nodes???? or 2:P.nz-1????
Se      = ( teta_pm1 - Psh.tetar(i)) ./ (Psh.tetas(i) - Psh.tetar(i));

% lambda shape parameter:
l       = Psh.bita; % unsure !!!!!!!!

% empirical shape factor (SWAP-32 manual, Eq.2.5, page 28):
%   see "dhconduc" function in functions.for
m       = 1-1./Psh.en(i);
dK_dSe  = Ksat .* Se.^(l-1) .* (1-(1-Se.^(1./m)).^m) .* ( l + (1-Se.^(1./m)).^(m-1) .* ((2+l).*Se.^(1./m)-l) );

% Capacity at "p" iteration, which substitutes the calculation of the
% derivative of theta to pressure head (SWAP-32 manual, Eq.2.8, page 28):
Ci       = multilayer_capacity(  h,  Psh, i );

% Derivative of hydraulic conductivity to the pressure head (SWAP-32,
% Appendix 5, page 239):
dKi_dhi = dK_dSe .* Ci ./ (Psh.tetas-Psh.tetar);
%% end
return