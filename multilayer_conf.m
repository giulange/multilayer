%%   Project Info
% ----------------------------------

proj.description    = 'MARWA CORRECTED NOVEMBRE 2011';

% ipath:                The path where all the input files required by
%                       multilayer simulation are in.
proj.ipath          = '/home/giuliano/git/multilayer';

% opath:                The path in which all the output printed files are
%                       stored during multilayer simulation. A check is
%                       performed to ensure that new simulations do not
%                       overwrite old ones, if not needed.
proj.opath          = '/home/giuliano/git/multilayer/output/';
%%   WATER INPUT
% ----------------------------------

% *SIMULATION & TIME*
% sdate:            Start date of simulation run. YYYY-MM-DD
W.sdate             =   '2013-01-01';       % 0
% edate:            End   date of simulation run. YYYY-MM-DD
W.edate             =   '2013-05-01';       % N
% timestep          Time step that can be used for input (meteo? only
%                   meteo?) data. 
%                   It can be day or hour. The total number of timesteps is
%                   given by:
%                       N = length( W.sdate:W.timestep:W.edate );
%                   The N number can be used as the maximum threshold with
%                   which define W.tp (and maybe other variables?).
%                   At this moment only day can be used (a future version
%                   will be available on a more detailed time step).
W.timestep          = 'day';

% tp:               Time steps at which print results. A value represents a
%                   fraction of day (e.g. 0.5 is half a day, and so on).
%                   The time step resolution to build the final time prints
%                   is the hour.
%                   -----------------------------------
%                       Time [day]          W.tp [-]
%                       [sdate,edate]       [0,...,N]
W.tp                = {
                        '2013-01-01,00'     1/24
                        '2013-01-02,00'     1.00
                      };
%                   -----------------------------------

% hin:              Initial condition for soil potential at any depth.
%                   -----------------------------------
%                       Depth [cm]      W.hin [?]
%                       [0,botlim]      [-100,+100]
W.hin               = [  
                        0               -100
                      ];
%                   -----------------------------------

% MTCL:                 Montecarlo simulation
%                           0:  no
%                           1:  yes
W.MTCL              = 0;

% isol:                 Indice per la simulazione del trasporto dei soluti:
%                           0:  non simulato
%                           1:  trasporto di soluti con Jury
%                           2:  trasporto di soluti con CDE
%                               (advection-dispersion)
W.isol              = 2;

% iosm:                 ?
W.iosm              = 0;
W.ads               = 1;
W.iveg              = 1;
W.dtin              = 0.00001;

% *SOIL GRID GEOMETRY*
% nlay:                 Number of soil layers.
W.nlay              = 3;
% zint:                 Soil layer bottom boundaries. The bottom depth of
%                       the lowest layer is the "botlim" value used to
%                       define the bottom Z-limit of all depth-dependent
%                       variables.
W.zint              = [25, 60, 300]; % e.g. botlim = 300 [cm]
% -------------------------------------------------------------------------
% VERTICAL DISCRETIZATION: [W.sg] --> sg='soil geometry'
% -------------------------------------------------------------------------
% type:                 Type of vertical discretization used to build the
%                       soil grid geometry.
%                       The following can be selected:
%                           > 1     --> regular grid spacing, according to
%                                       the 1/2 distance at soil layers
%                                       saddle points. No more parameters
%                                       must be defined apart from 'nz',
%                                       'nlay' and 'zint'.
%                           > 2     --> "sublayers" spacing, in which each
%                                       layer can be divided into more
%                                       layers within which the spacing is
%                                       regular. More parameters must be
%                                       defined
%                           > 3     --> any other kind of geometry we would
%                                       like to implement!
W.sg.type           = 1;
 
% regular:              It creates a soil grid with regular node spacing,
%                       at least within the same soil layer.
%                       Node spacing is quite similar between soil layers,
%                       but a small adjustment is made in order to fit
%                       layer thickness with the putative number of nodes
%                       and to make the bottom layer boundary as the
%                       bisector of the two adjacent nodes.
%                       Two following parameters are required:
%                        -numnodes: The number of nodes for the whole soil
%                                   profile defining the geometry of the
%                                   soil column during simulation of water
%                                   flow. It includes the top and the botom
%                                   boundaries.
%                       %# numnodes #%
W.sg.regular        = 100;

% sublayers:            Define a specific grid in which each sublayer can
%                       be different from the others but within which nodes
%                       spacing is regular.
%                       #CHECK WITH ANTONIO#
%                       You have to define equal nodes spacing across soil
%                       layers boundaries (i.e. the soil layer bottom
%                       boundary must bisect the two adjacent nodes!).
%                       Remeber that the cumulative sum of all hSubLay-s
%                       must equal the bottom depth of the last soil layer,
%                       and that:
%                           hSubLay = hNode x nNodes
% 
%                       User has to define the following parameters:
%                        -SoilLay:  The number of soil layer at which
%                                   one or more sub-layers can be defined.
%                        -SubLay:   The global number of sublayer in soil
%                                   profile.
%                        -hSubLay:  The thickness of current sublayer,
%                                   which can be discretized into one or
%                                   more nodes.
%                        -hNode:    The thickness of each node constituting
%                                   the sublayer.
%                        -nNodes:   The number of regularly spaced nodes
%                                   that discretise the sublayer.
%                     ---------------------------------------------------
%                       %#  SoilLay SubLay hSubLay   hNode  nNodes #%
%                              [-]   [-]    [cm]     [cm]    [-]
%                     ---------------------------------------------------
W.sg.sublayers      = [         1     1      5.0      1.0     5  
                                1     2     15.0      3.0     5
                                1     3      5.0      1.0     5 % 25 cm
                                2     4      3.0      1.0     3
                                2     5     30.0     10.0     3
                                2     6      2.0      2.0     1 % 35 cm
                                3     7     10.0      2.0     5
                                3     8     10.0      5.0     2
                                3     9    200.0     40.0     5
                                3    10     20.0      5.0     4 % 240 cm
                       ];                  % TOT=               % 300 cm
%                     ---------------------------------------------------

% anotherkind:          Whatever we want to implement (we should check what
%                       HYDRUS makes as a suggestion).
W.sg.anotherkind    = [];
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% SOIL GRID NODES WITH HYDRAULIC CHARACTERISTICS:
% -------------------------------------------------------------------------
% W.crc.? -->   Curva ritenzione/conducibilità: es. W.crc.dap, W.crc.tetas,
%               ecc.
%                   -------------------------------------------------------
% #SoilLay#             1        2       3      4       5       6       7
%                   -------------------------------------------------------
W.dap               = [ 1.1,     1.1,    1.1     ];
W.tetas             = [ 0.340,   0.310,  0.300   ];
W.tetar             = [ 0.000,   0.000,  0.000   ];
W.alfrs             = [ 0.000,   0.000,  0.000   ];
W.fi                = [ 0.000,   0.000,  0.000   ];
W.alfvg             = [ 0.120,   0.140,  0.150   ];
W.en                = [ 1.120,   1.140,  1.250   ];
W.alfvg2            = [ 0.0000,  0.0000, 0.0000  ];
W.en2               = [ 0.000,   0.000,  0.000   ];
W.ifr               = [ 1,       1,      1       ];
W.k0                = [ 50.00,   20.00,  20.00   ];
W.k0macr            = [ 0.000,   0.000,  0.000   ];
W.bita              = [ 0.5,     0.5,    0.5     ];
W.bita2             = [ 9999,    9999,   9999    ];
W.ifc               = [ 1,       1,      1       ];
%                   -------------------------------------------------------

% ??
W.hfc               = -333; % --> check with Antonio!!
W.vpr               = 0.5;
W.tetal             = 0.0002;
W.bital             = 15.0;
W.hsurfmax          = 0.0;
% -------------------------------------------------------------------------

W.itopvar           = 1;
W.ibotvar           = 0;

% iCtopvar          Indice per la lettura dei dati di concentrazione al
%                   contorno superiore:
%                       0:  valore di concentrazione Cinput
%                           (solute_CDE_inp.txt)
%                       1:  condizioni al contorno superiore variabili
%                           (Ctopbound_inp.txt)
W.iCtopvar          = 1;

% itbc:             Kind of top boundary condition:
%                       0   --> flux
%                       1   --> potential
W.itbc              = 0;

% ibbc:             Kind of top boundary condition:
%                       0   --> flux
%                       1   --> potential
%                       2   --> gradient
W.ibbc              = 2;

% Can I use W.hqsurf instead of the two following???
W.hsurf             = 9999;
W.qsurf             = 0.1;

W.hbot              = -0.0;
W.qbot              = 9999;
W.grad              = 1.0;

W.inhin             = 0;
% W.inhin     = 1;
% W.hin       = load( fullfile(proj.ipath,'initial_inp.txt') )

% W.tmax              = 120;
% W.ntp               = 131;
%%   TOP BOUNDARY INPUT
% ----------------------------------
B.top.description   = 'MARWA CORRECTED NOVEMBRE 2011';%--> 'a discretization...'
B.top.ntbc          = 126;
%%   BOTTOM BOUNDARY INPUT
% ----------------------------------
B.bot.description   = 'MARWA CORRECTED botboundary';%--> 'a discretization...'
B.bot.nbbc          = 10;
%%   CONCENTRATION TOP BOUNDARY INPUT
% ----------------------------------
% 
% INPUT di SOLUTI
%   { FR=fertirr, SD=solido, UR=urea, 
%     ORG_rp=organico a mineralizzazione rapida,
%     ORG_sw=organico a mineralizzazione lenta }
% Si assume che FR sia liquido e rappresenti l'apporto in superficie
% (C_input nell'equazione ADE) mentre le altre forme si considerano
% distribuite su uno spessore dL ed entrano nell'ADE come sink-source.
% Tstar è la temperatura in °C per il calcolo di Kmineralizzazione, sia
% rapida che lenta.

B.Ctop.description  = 'MARWA CORRECTED Ctopbound';%--> 'a discretization...'
B.Ctop.nCtop        = 126; % = B.top.ntbc !!!
% measurement units??
B.Ctop.KhUR         = 1.000; % 
B.Ctop.KvUR         = 1.000;
B.Ctop.KmORG_rp     = 0.020; % andrebbe messo il punto tra "KmORG" ed "rp"
B.Ctop.KmORG_sw     = 0.002; % andrebbe messo il punto tra "KmORG" ed "sw"
B.Ctop.dL           = 30;   % [cm]
%%   VEGETATION INPUT
% ----------------------------------
V.description       = 'prova Lodi Arm_Art vegetation';%--> 'the plant used was...'
V.nET               = 126;
% extf:             Esponente legge Beers.
V.extf              = 0.6;
% ifs:              Indicatore funz. sink {Feddes, vanGen.}
V.ifs               = 1;
% ifg:              Indicatore funz. distribuz. appar.rad.
V.ifg               = 1;
% hI:               Pot. stress idrico Feddes
V.hI                = -1;
% hII:              Pot. stress idrico Feddes
V.hII               = 10;
% hIIIH:            Pot. stress idrico Feddes
V.hIIIH             = -400;
% hIIIL:            Pot. stress idrico Feddes
V.hIIIL             = -600;
% hIV:              Pot. stress idrico Feddes
V.hIV               = -8000;
% hw50:             Pot. idrico dimezzamento traspirazione van Genuchten
V.hw50              = -1000;
% pw1:              Esponente stress idrico van Genuchten
V.pw1               = 3;
% hs50:             Pot.osmotico dimezzamento traspirazione van Genuchten
V.hs50              = -1500;
% ps1:              Esp.stress osmotico van Genuchten
V.ps1               = 3;
% aMH:              Par. stress osmotico Mass & Hofmann
V.aMH               = -760;
% bMH:              Par. stress osmotico Mass & Hofmann
V.bMH               = 0.000794;
% rda:              Par. distribuzione radici logistica
V.rda               = 1.027;
% rdb:              Par. distribuzione radici logistica
V.rdb               = 15.016;
% rdc:              Par. distribuzione radici logistica
V.rdc               = 0.074;
% zc:               Par. distribuzione radici doppia-lineare
V.zc                = 25;
% g0:               Par. distribuzione radici doppia-lineare
V.g0                = 0.032;
% gzc:              Par. distribuzione radici doppia-lineare
V.gzc               = 0.008;
% Drf:              Par. distribuzione radici doppia-lineare
V.Drf               = 85;
%%   SOLUTE Jury INPUT -- translate definitions in variables that can be used by the program
% ----------------------------------
if W.isol==1
S.J.description     = 'prova puglianello solute transport';
S.J.decad           = 0;
S.J.retard          = 0;

% tetasst:          ??
%                   --------------------------------------
%                       Depth [cm]      S.J.tetasst [??]
%                       [0,botlim]      [0.000,+10.000]
S.J.tetasst         = [
                        0               0.500
                      ];
%                   --------------------------------------

% sigma:            ??
%                   --------------------------------------
%                       Depth [cm]      S.J.tetasst [??]
%                       [0,botlim]      [0.000,+10.000]
S.J.sigma           = [
                        0               0.200
                      ];
%                   --------------------------------------

% Rcoeff:           ??
%                   --------------------------------------
%                       Depth [cm]      S.J.Rcoeff [-]
%                       [0,botlim]      [a,b]
S.J.Rcoeff          = [
                        0               1.000
                      ];
%                   --------------------------------------

% Dcoeff:           ??
%                   --------------------------------------
%                       Depth [cm]      S.J.Dcoeff [-]
%                       [0,botlim]      [a,b]
S.J.Dcoeff          = [
                        0               0.000
                      ];
%                   --------------------------------------

end
%%   SOLUTE CDE INPUT
% ----------------------------------
if W.isol==2
S.CDE.description   = 'MARWA CORRECTED';

% tCinput:          Tempo di immissione del soluto.
%                   [??]
S.CDE.tCinput       = 0.000;

% tCinput_end:      ??
S.CDE.tCinput_end   = 0.500;

% Cinput:           Input concentration.
%                   [g cm-3]
S.CDE.Cinput        = 0.040;
% Topt:             The optimum temperature for the XXX process
%                   [°C]
S.CDE.Topt          = 25;

% ***Freundlich Isotherm***
% The Freundlich isotherm is here used on NX, where X can be 'H' for
% AMMONIA (NH4) or 'O' for NITRATE (NO3).
% Two coefficients must be defined below to use Freundlich method: kf1 and
% kf2.
% It is assumed that the isotherm is made by two segments with two
% different slopes (i.e. kf1 and kf2).
% 
% NX.kf1:           Pendenza del primo tratto dell'isoterma di Freundlich.
%                   [cm3 g-1]
%                      'NH'   'NO'
S.CDE.NX.Kf1        = [0.100, 0.000];
% 
% NX.kf2:           Pendenza del secondo tratto dell'isoterma di
%                   Freundlich.
%                   [cm3 g-1]
%                      'NH'   'NO'
S.CDE.NX.Kf2        = [1.000, 1.000];

% NX.Kr:            Fattore di attingimento radicale.
%                   [-]
S.CDE.NX.Kr         = [1.000, 1.000];


% Cin.NH:           Initial concentration of AMMONIA, defined at different
%                   depths.
%                   --------------------------------------------
%                       Depth [cm]      S.CDE.Cin.NH [g cm-3]
%                       [0,botlim]      [0.0000,+10.0000]
S.CDE.Cin.NH        = [
                        0               0.0000
                      ];
%                   --------------------------------------------


% Cin.NO:           Initial concentration of NITRATE, defined at different
%                   depths.
%                   --------------------------------------------
%                       Depth [cm]      S.CDE.Cin.NO [g cm-3]
%                       [0,botlim]      [0.000,+10.000]
S.CDE.Cin.NO        = [
                        0               0.0002
                        30              0.0000
                      ];
%                   --------------------------------------------


% lambda:           Dispersivity of??, defined at different depths.
%                   --------------------------------------------
%                       Depth [cm]      S.CDE.lambda [cm]
%                       [0,botlim]      [a,b]
S.CDE.lambda        = [
                        0               3.0
                      ];
%                   --------------------------------------------


% Knitr:            Nitrification coefficient.
%                   --------------------------------------------
%                       Depth [cm]      S.CDE.Knitr [-]
%                       [0,botlim]      [a,b]
S.CDE.Knitr         = [
                        0               1.0000
                        30              0.0100
                      ];
%                   --------------------------------------------


% Kimmob:           Immobilization coefficient.
%                   --------------------------------------------
%                       Depth [cm]      S.CDE.Kimmob [-]
%                       [0,botlim]      [a,b]
S.CDE.Kimmob        = [
                        0               0.0400
                        30              0.0300
                      ];
%                   --------------------------------------------


% Kdenitr:          Denitrification coefficient.
%                   --------------------------------------------
%                       Depth [cm]      S.CDE.Kdenitr [-]
%                       [0,botlim]      [a,b]
S.CDE.Kdenitr       = [
                        0               0.0400
                        30              0.0200
                      ];
%                   --------------------------------------------

end
%%   EC DATA -- modify according to the new definition of depth-dependent variables
% ----------------------------------
% 
% Lettura valori di potenziale osmotico
%   (100 nodi x numero tempi di misura)
% Nel file fnsink, gli V.ifs>3 si riferiscono ai casi di attingimento con
% stress salino, accoppiato o non allo stress idrico. In presenza di stress
% salino, occorre far leggere un file di input con i dati di EC misurati
% per i 100 nodi di calcolo.
if W.iosm==1 && V.ifs>3
% load DATA
EC.matrix           = load(fullfile(proj.ipath,'EC_data.txt'));
EC.z                = EC.matrix(2:end,1);
EC.t                = EC.matrix(1,2:end);
EC.data             = EC.matrix(2:end,2:end);
end
%%   MONTECARLO
% ----------------------------------
% to A.Basile:  Why not to consider the in-between combination of
%               parameters?
% i.e.          The cartesian product is not by "rows" of montecarlo.txt,
%               but by every single parameter to be stochastically varied.

if W.MTCL==1
% nlay:             Define how the list of stochastically defined
%                   parameters in M.data are distributed in the W.nlay
%                   strata. If a soil layer is missing, set its value to
%                   zero (e.g. [22,22,0] means that the third layer will
%                   not be considered as stochastic, but the configuration
%                   given in W is taken for that soil layer).
%                   -------------------------------------------------------
% #SoilLay#             1        2       3      4       5       6       7
%                   -------------------------------------------------------
M.nlay              = [ 22,     33,     11];

% combinatorial:    A Montecarlo combinatorial calculus is performed to set
%                   all the possible combinations of each soil layer
%                   stochastic repetitions with others.
%                       0:  sequntial (non-combinatorial).
%                       1:  combinatorial
M.combinatorial     = 1;

% nvp:              Number of stochastic simulations.
%                   Its value must be set only for non-combinatorial
%                   Montecarlo simulations (it will be ignored anyway).
%                   The program runs taking the first nvp stochastic
%                   repetitions from M.data for each soil layer (defined in
%                   W.nlay).
%                   If M.nvp=22 and M.nlay=[22,22,11] it means that only
%                   the first 11 repetitions are taken from M.data
%                   (considering that M.nlay(3)=11 is the limiting factor
%                   for using all M.nvp=22 repetitions).
M.nvp               = 22;

% list:             List of variables that must be simulated with
%                   stochastic Montecarlo and that are listed/defined in
%                   the montecarlo.txt file. For instance:
%                       = { 'W.tetas', 'W.zint' };
M.list              = { 'W.tetas','W.tetar','W.alfvg','W.en','W.k0' };
% data:             Configuration for each Montecarlo repetition.
M.data              = load(fullfile(proj.ipath,'montecarlo.txt'));

% nnc:              Number of nonconvergences (niter>10,dt<=W.dtmin)
M.nnc               = 0;
% -----------------------------
end