%% pre-allocate all variables required within each Monte Carlo Simulation

% to be rationalized:
macr                            = zeros(P.nz,1);
dSINKmacr_dh                    = zeros(P.nz,1);

% time:
P.j                             = 1;
P.dt                            = NaN;
P.tidx                          = 0;

% -GEOMETRY
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
P.K                             = NaN(P.nz,1);
P.cap                           = NaN(P.nz,1);
P.sink                          = NaN(P.nz,1);
P.kp                            = NaN(P.nz,1);
P.flux                          = NaN(P.nz,1);
P.h1                            = NaN(P.nz,1);
P.h1star                        = NaN(P.nz,1);
P.h2                            = NaN(P.nz,1);
P.Kim2                          = NaN(P.nz,1);
P.Kip2                          = NaN(P.nz,1);
% Fraction of horizontal area of soil matrix per compartment (-):
P.FrArMtrx                      = NaN(P.nz,1);
Fi                              = NaN(P.nz,1);
dKim2_dhi                       = NaN(P.nz,1);
dKip2_dhi                       = NaN(P.nz,1);
dKim2_dhim1                     = NaN(P.nz,1);
dKip2_dhip1                     = NaN(P.nz,1);
dF_dh                           = zeros(P.nz,P.nz);

% --others:
P.ECstar                        = NaN(P.nz,1);
P.C1                            = NaN(P.nz,2);

% -TIME:
P.time                          = NaN(1,P.Nj); 
P.tidx_jm1                      = 0;
P.L                             = 0;
P.flEndOfDay                    = false;
P.km_max                        = NaN(1,P.Nj); 
P.fluxsurf_max                  = NaN(1,P.Nj); 
P.km                            = NaN(1,P.Nj); 
% [p;nr_breaked;fl_noconv;n_noconv;iL;bt_breaked]
P.iter                          = NaN(6,P.Nj);

% -SCALARS:
% --counters:
P.k                             = NaN;
% P.niter                         = NaN(1,P.Nj);% might be useful!
% --others:
P.rnf                           = 0;
P.ktec                          = 1;
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
