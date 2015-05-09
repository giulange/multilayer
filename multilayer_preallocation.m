%% pre-allocate all variables required within each Monte Carlo Simulation

% -GEOMETRY
% % --nodes
% P.nodes.num                     = NaN( W.nlay+0,    1 );
% P.nodes.thickness               = NaN( W.nlay+0,    1 );
% P.nodes.cumsum                  = NaN( W.nlay+1,    1 );
% P.nodes.soillayer               = NaN( P.nz+1,      1 );
% P.nodes.z                       = NaN( P.nz+1,      1 );
% P.nodes.dz                      = NaN( P.nz+1,      1 );
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
% --hydr (put in .sh.)
P.teta                          = NaN(P.nz,1);
P.kond                          = NaN(P.nz,1);
P.cap                           = NaN(P.nz,1);
P.sink                          = NaN(P.nz,1);
P.kp                            = NaN(P.nz,1);
P.flux                          = NaN(P.nz,1);
P.h1                            = NaN(P.nz,1);
P.h1star                        = NaN(P.nz,1);
P.h2                            = NaN(P.nz,1);
% --others:
P.ECstar                        = NaN(P.nz,1);
P.C1                            = NaN(P.nz,2);

% -TIME:
P.time                          = NaN(1,P.Nj); % **USEFULL on Nj
P.km_max                        = NaN(1,P.Nj); % **USELESS on Nj
P.fluxsurf_max                  = NaN(1,P.Nj); % **USELESS on Nj
P.km                            = NaN(1,P.Nj); % **USELESS on Nj
% [p;nr_breaked;fl_noconv;n_noconv;iL;bt_breaked]
P.iter                          = NaN(6,P.Nj);

% -SCALARS:
% --counters:
P.j                             = 1;
P.dt                            = NaN;
P.k                             = NaN;
% P.niter                         = NaN(1,P.Nj);% might be useful!
% --others:
P.Emax                          = NaN;
P.teta_hsurf                    = NaN;
P.teta_hbot                     = NaN;
%   initialized to zero (SWAP initialise pond in SoilWater(1), line 173)
% P.pond                          = 0;
P.dpt                           = NaN;
P.op                            = NaN;
P.ETp0                          = NaN;
P.Ep                            = NaN;
P.Tp                            = NaN;
P.Droot                         = NaN;
