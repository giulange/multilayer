%% pre-allocate all variables required within each Monte Carlo Simulation

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
% --hydr
P.teta                          = NaN(P.nz,1);
P.kond                          = NaN(P.nz,1);
P.cap                           = NaN(P.nz,1);
P.sink                          = NaN(P.nz,1);
P.kp                            = NaN(P.nz,1);
P.flux                          = NaN(P.nz,1);
P.h1                            = NaN(P.nz,1);
P.h1star                        = NaN(P.nz,1);
% --others:
P.ECstar                        = NaN(P.nz,1);
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