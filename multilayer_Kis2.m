%% NOTES
% The script calculates the internodal conductivities at i-1/2 (Kim2) and
% at i+1/2 (Kip2).
% The letter "s" in Kis2 stays for "sign" which can be minus ("m") or plus
% ("p").
%% init
Kim2            = zeros(P.nz,1);
Kip2            = zeros(P.nz,1);
%% **INTERNODAL CONDUCTIVITIES**
%% Internodal conductivity at intermediate nodes:
for ii = 2:P.nz-1% avoid loop defining multilayer_conductivity_internode.m more generally!!!
    Kim2(ii)    = multilayer_conductivity_internode( K(ii-1), K(ii), W.Kmeth, dz_im1(ii-1), dz_i(ii-1) );
    Kip2(ii)    = multilayer_conductivity_internode( K(ii), K(ii+1), W.Kmeth, dz_i(ii-1), dz_ip1(ii-1) );
end
%% Internodal conductivity at top    node (only Kip2 can be defined here):
Kip2(1)         = multilayer_conductivity_internode( K(1), K(2), W.Kmeth, dz(1), dz(2) );
%Kim2(1)        --> computed in multilayer_boundtop.m
%% Internodal conductivity at bottom node (only Kim2 can be defined):
Kim2(P.nz)      = multilayer_conductivity_internode( K(P.nz-1), K(P.nz), W.Kmeth, dz(P.nz-1), dz(P.nz) );
%Kip2(P.nz)     --> maybe should be defined in multilayer_boundbot.m !!!
