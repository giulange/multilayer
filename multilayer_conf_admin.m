W.tolle1            = 0.0001;
W.dtmin             = 0.00000001;
W.dtmax             = 0.5;
W.multmin           = 1.2;
W.multmax           = 0.7;

% Nj_shp:           Set how number of iterations grows according to
%                   increasing W.tmax values. For small values of W.tmax or
%                   for fast convergence case studies W.Nj_shp can converge
%                   to 2.6, otherwise it must decrease towards 1.0 (rarely
%                   below 1.0). A value of 1.5 is advised for most
%                   applications.
W.Nj_shp            = 2.0;

% files:            Set all the filenames used to save the following
%                   variables:
%                   { O.C2, O.h22, O.fluxsurf, O.fluxbot, O.runoff }
% 
% h22:              Print of flux at intermediate nodes.
%                   [nodes x times x montecarlo]
O.files.h22         = 'potential.txt';
% 
% C2:               Print of solutes concentrations.
%                   [nodes x times x montecarlo x solutes]
O.files.C2          = 'concentration.txt';
% 
% flux:             Print of fluxes at top (=surface) and bottom boundaries
%                   and of runoff (runon is not yet considered).
%                   [1 x times x montecarlo]
O.files.flux        = 'flux.txt';
