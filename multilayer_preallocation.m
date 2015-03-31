%% pre-allocate all variables required within any Monte Carlo Simulation

% ****TEMPORARY****
% --nodes
P.nodes.num                     = NaN( W.nlay+0,    1 );
P.nodes.thickness               = NaN( W.nlay+0,    1 );
P.nodes.cumsum                  = NaN( W.nlay+1,    1 );
P.nodes.soillayer               = NaN( P.nz+1,      1 );
P.nodes.z                       = NaN( P.nz+1,      1 );
P.nodes.dz                      = NaN( P.nz+1,      1 );
% --soil-grid
P.sh.dap                        = NaN(P.nz,1);
P.sh.tetas                      = NaN(P.nz,1);
P.sh.tetar                      = NaN(P.nz,1);
P.sh.alfrs                      = NaN(P.nz,1);
P.sh.fi                         = NaN(P.nz,1);
P.sh.alfvg                      = NaN(P.nz,1);
P.sh.en                         = NaN(P.nz,1);
P.sh.alfvg2                     = NaN(P.nz,1);
P.sh.en2                        = NaN(P.nz,1);
P.sh.ifr                        = NaN(P.nz,1);
P.sh.k0                         = NaN(P.nz,1);
P.sh.k0macr                     = NaN(P.nz,1);
P.sh.bita                       = NaN(P.nz,1);
P.sh.bita2                      = NaN(P.nz,1);
P.sh.ifc                        = NaN(P.nz,1);
P.sh.tetafc                     = NaN(P.nz,1);
% --scalars:
% ---counters:
P.j                             = NaN;
P.jstar                         = NaN;
P.k                             = NaN;
P.SS                            = NaN;
P.niter                         = NaN;
% ---others:
P.dpt                           = NaN;
P.op                            = NaN;
P.teta_hsurf                    = NaN;
P.Emax                          = NaN;
P.teta_hbot                     = NaN;
% --vectors:
P.h2                            = NaN(P.nz,1);
% ---time:
P.time                          = NaN(1,P.Nj);
P.km_max                        = NaN(1,P.Nj);
P.fluxsurf_max                  = NaN(1,P.Nj);
P.km                            = NaN(1,P.Nj);
% ---others:
P.Cinput                        = NaN(2,1);
% --arrays: % check with Antonio --> I would DELETE the P.Nj dimension!!
P.C1                            = NaN(P.nz,P.Nj,2);
P.h1                            = NaN(P.nz,2,P.Nj); % **check with Antonio
P.h1star                        = NaN(P.nz,2,P.Nj); % **check with Antonio
P.ECstar                        = NaN(P.nz,P.Nj);
P.teta                          = NaN(P.nz,P.Nj);
P.kond                          = NaN(P.nz,P.Nj);
P.cap                           = NaN(P.nz,P.Nj);
P.sink                          = NaN(P.nz,P.Nj);
P.kp                            = NaN(P.nz,P.Nj);
P.flux                          = NaN(P.nz,P.Nj);
