function [dKi_dhi,Ci,teta_p] = multilayer_dKi_dhi(h, Psh, nodes, dt)
% [dKi_dhi,Ci,teta_p] = multilayer_dKi_dhi(h, Psh, nodes, dt)
% 
% NOTE
%   This function is similar to the dhconduc function of SWAP-32, used in
%   headcalc.for and implemented in functions.for.
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
%   nodes:      One or more nodes at which outputs are computed.
% 
%   dt:         Simulation time step.
% 
% OUTPUT
%   dKi_dhi:    Derivative of hydraulic conductivity to the pressure head
%               at node(s) i.
%   Ci:         Capillary Capacity at node(s) i.
%   teta_p:     Moisture fraction at node(s) i.

%% pre
Ksat        = Psh.k0(nodes);
% lexp shape parameter:
l           = Psh.bita(nodes);
% empirical shape factor (SWAP-32 manual, Eq.2.5, page 28):
thetar      = Psh.tetar(nodes);
thetas      = Psh.tetas(nodes);
alfamg      = Psh.alfvg(nodes);
n           = Psh.en(nodes);
m           = 1-1./Psh.en(nodes);
h_enpr      = Psh.h_enpr(nodes);
%% main
% Moisture at "p" iteration:
% teta_p      = fnteta( h, Psh, nodes );
teta_p      = multilayer_fnteta_vgm( h, Psh, nodes );

% Approximation of the differential water Capacity [cm-1] at "p" iteration,
% substituted by the calculation of the derivative of theta to pressure
% head (SWAP-32 manual, Eq.2.8, page 28). The following calculation
% accounts for modification of VG retention curve:
Ci          = multilayer_capacity(  h,  Psh, nodes, dt );% dimocap

% Relative saturation (SWAP-32, Eq. 2.7, page 28):
% check that:
%   -h_pm1 is what you need to define starting condition!
%   -all nodes???? or 2:P.nz-1????
Se          = ( teta_p - thetar) ./ (thetas - thetar);% relsat
%% conductive derivative to relative water content:

% % % **implemented by Giuliano:
% % %  -see Appendix 5, page 240:
% % dK_dSe  = Ksat .* Se.^(l-1) .* (1-(1-Se.^(1./m)).^m) .* ( l + (1-Se.^(1./m)).^(m-1) .* ((2+l).*Se.^(1./m)-l) );
% % %         |___term5_______|    |_____term6_________|          |____term3__________|    |_____term2_________|
% % % 
% % %   correct the possible Inf values to prevent numerical instability:
% % %       > ask to Antonio if it is correct to assume this kind of
% % %         modification, and if it is right to apply it to all compartments
% % % dK_dSe( isinf(dK_dSe) ) = 1.0d-12;

% **alternative, as implemented in SWAP:
% %  -see "dhconduc" function in functions.for in SWAP-32;
Sc      = (  1 + abs( alfamg .* h_enpr ).^n  ) .^ (-m);% s_enpr
term0   = ( Se.* Sc ).^ (1./m);
term4   = 1-( 1-(Sc) .^ (1./m) ).^m;
term6   = 1-(1-term0).^m;
term2   = (2+l).*term0 -l;
term3   = (1-term0).^(m-1);
dK_dSe  = Ksat .* Se.^(l-1) .* term6 .* ( l + term3 .* term2 ) ./ term4.^2;
%         |___term5_______|
%% hydraulic conductivity derivative to pressure head (SWAP-32, Appendix 5,
% page 239):
dKi_dhi = dK_dSe .* Ci ./ (thetas-thetar);
% adaptation in case of dry conditions
dKi_dhi(Se < 0.001d0) = 0;
% ?
dKi_dhi(Se > Sc)  = 1.0d-12;
%% end
return