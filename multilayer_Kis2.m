%% NOTES
% The script calculates the internodal conductivities at i-1/2 (P.Kim2) and
% at i+1/2 (P.Kip2).
% The letter "s" in Kis2 stays for "sign" which can be minus ("m") or plus
% ("p").
% I SHOULD AVOID 2 Kis2 AND CALCULATE ONLY ONE Kmean!!
%% **NODAL CONDUCTIVITIES**
% hydraulic conductivity, treated implictly (SWAP-32, Eq. 2.6, page
% 28):
% NOTES:
%   -you should pass both h and teta to multilayer_conductivity_node and
%   internally use ifc to compute P.K(i).
x               = P.teta;
notTeta         = P.sh.ifc~=1 & P.sh.ifc~=3;
x(notTeta)      = P.h(notTeta);
P.K             = multilayer_conductivity_node( x, P.sh, 1:P.nz );
if W.SwMacro
    P.K = P.FrArMtrx .* P.K;
end
%% **INTERNODAL CONDUCTIVITIES**
%% Internodal conductivity at intermediate nodes:
for ii = 2:P.nz-1% avoid loop defining multilayer_conductivity_internode.m more generally!!!
    P.Kim2(ii)    = multilayer_conductivity_internode( P.K(ii-1),  P.K(ii),   W.Kmeth, P.nodes.dz(ii-1),   P.nodes.dz(ii) );
    P.Kip2(ii)    = multilayer_conductivity_internode( P.K(ii),    P.K(ii+1), W.Kmeth, P.nodes.dz(ii),     P.nodes.dz(ii+1) );
end
%% Internodal conductivity at top    node (only P.Kip2 can be defined here):
P.Kip2(1)         = multilayer_conductivity_internode( P.K(1),     P.K(2),    W.Kmeth, P.nodes.dz(1),      P.nodes.dz(2) );
%P.Kim2(1)        --> computed in multilayer_boundtop.m
%% Internodal conductivity at bottom node (only P.Kim2 can be defined):
P.Kim2(P.nz)    = multilayer_conductivity_internode( P.K(P.nz-1),P.K(P.nz), W.Kmeth, P.nodes.dz(P.nz-1), P.nodes.dz(P.nz) );
P.Kip2(P.nz)    = P.K(P.nz);% <--- this is used in Fi(NN) for free drainage
clear x notTeta