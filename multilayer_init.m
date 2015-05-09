%% initialization of enviroment variables, made before Monte Carlo Simulation

% We need a standard:
%   possible indexes are --> { node, time, montecarlo, solute? }
%   These are the rules:
%   ---to-be-saved---
%   (1) node, (2) time, (3) montecarlo, (4) solute;
%   (1) node, (2) time, (3) montecarlo;
%   (1) 1,    (2) time, (3) montecarlo;
%   ---temporary---
%   (1) node, (2) time;
%   (1) 1,    (2) time;
%   (1) node, (2) 1;
%   other...

warning('You have to decide whether to allocate for P.nz or for P.nz+1 !!')

% ****TO BE SAVED****
%   > O.C2          --> solutes concentrations
%                       [nodes x times x montecarlo x solutes]
%   NOTE: I have to substitute C2 with a variable specific for any solute
%         we would like to implement in the model.
O.C2                            = NaN( P.nz, P.Nj, M.nvp, 2 );
%   > O.h22         --> flussi ai nodi intermedi
%                       [nodes x times x montecarlo]
O.h22                           = NaN( P.nz, P.Nj, M.nvp );
%   > O.fluxsurf    --> flusso al contorno superiore
%                       Water flux through soil surface?
%                       [1 x times x montecarlo]
O.fluxsurf                      = NaN( 1, P.Nj, M.nvp );
%   > O.fluxbot     --> flusso al contorno inferiore
%                       [1 x times x montecarlo]
O.fluxbot                       = NaN( 1, P.Nj, M.nvp );
%   > O.runoff      --> runoff
%                       [1 x times x montecarlo]
O.runoff                        = NaN( 1, P.Nj, M.nvp );
%   > O.pond        --> ponding
%                       [1 x times x montecarlo]
O.pond                          = NaN( 1, P.Nj, M.nvp );
