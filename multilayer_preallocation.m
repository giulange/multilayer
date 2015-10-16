%% pre-allocate all variables required within each Monte Carlo Simulation

% initialize vegetation
init_crop_vars

% to be rationalized:
macr                            = zeros(P.nz,1);
dSINKmacr_dh                    = zeros(P.nz,1);
% ArMpSs:                       Area fraction of macropores at soil
%                               surface. I set it to zero so that it is not
%                               accounted (untill I will implement
%                               macropore module.
ArMpSs                          = 0.d0;%        [-]

% time:
P.sdate                         = datenum(W.sdate);
P.edate                         = datenum(W.edate);
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
P.sh.h_enpr                     = NaN(P.nz,1);
P.sh.ifc                        = NaN(P.nz,1);
P.sh.tetafc                     = NaN(P.nz,1);
% --hydr (put in .sh.)
P.teta                          = NaN(P.nz,1);
P.K                             = NaN(P.nz,1);
P.cap                           = NaN(P.nz,1);
P.sink                          = NaN(P.nz,1);
P.kp                            = NaN(P.nz,1);
P.q                             = NaN(P.nz+1,1);
P.h_jm1                         = NaN(P.nz,1);
P.h                             = NaN(P.nz,1);
P.Kim2                          = NaN(P.nz,1);
P.Kip2                          = NaN(P.nz,1);
P.q                             = NaN(P.nz,1);
% Fraction of horizontal area of soil matrix per compartment (-):
P.FrArMtrx                      = NaN(P.nz,1);
Fi                              = NaN(P.nz,1);
dKim2_dhi                       = NaN(P.nz,1);
dKip2_dhi                       = NaN(P.nz,1);
dKim2_dhim1                     = NaN(P.nz,1);
dKip2_dhip1                     = NaN(P.nz,1);
dF_dh                           = zeros(P.nz,P.nz);
% crop parameters:
P.Kc                            = zeros(1,W.tmax);
P.Ke                            = zeros(1,W.tmax);
P.Droot                         = zeros(1,W.tmax);
P.LAI                           = zeros(1,W.tmax);
% irrigation parameters:
P.irri                          = zeros(1,W.tmax);
% cirr stores NH4+ (col=1) and NO3- (col=2) concentrations of irrigation water: 
P.cirr                          = zeros(2,W.tmax);% zero because it is multiplied by net irrigation depth! 
P.irTYPE                        = NaN(1,W.tmax);
P.irDEPTH                       = NaN(1,W.tmax);
P.h_from                        = NaN(1,W.tmax);
P.h_to                          = NaN(1,W.tmax);
P.t_DRY                         = 0;
%   -critical root zone pressure head, useful for irrigation requirement:
P.hrz_cm                        = NaN(1,W.tmax);
%   -water requirement of the critical root zone, [cm]:
P.dstor                         = NaN(1,W.tmax);
%   -other:
P.ECstar                        = NaN(P.nz,1);
% solute parameters:
P.cml                           = NaN(P.nz,2);
P.samini                        = NaN(1,2);
P.cmsy                          = NaN(P.nz,2);
P.solbal                        = NaN(2,P.Nj);
P.Ndtsolute                     = zeros(1,P.Nj);% need 0 as first value!

% -TIME:
P.time                          = NaN(1,P.Nj); 
P.tidx_jm1                      = 0;
P.L                             = 0;
P.flEndOfDay                    = false;
P.flStartOfDay                  = true;
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
P.ptra                          = NaN;
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
