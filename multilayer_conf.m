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

% SOIL GRID GEOMETRY:
% nlay:                 Number of soil layers.
W.nlay              = 3;
% zint:                 Soil layer bottom boundaries.
W.zint              = [25, 60, 300];
% VERTICAL DISCRETIZATION: [W.sg] --> sg='soil geometry'
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
W.sg.regular        = [     50          ];
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
% 
%                       %#  SoilLay SubLay hSubLay   hNode  nNodes #%
%                              [-]   [-]    [cm]     [cm]    [-]
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
% anotherkind:          Whatever we want to implement (we should check what
%                       HYDRUS makes as a suggestion).
W.sg.anotherkind    = [];

% SOIL GRID NODES WITH HYDRAULIC CHARACTERISTICS:
% W.crc.? --> curva ritenzione/conducibilità: es. W.crc.dap, W.crc.tetas, ecc. 
W.dap               = [1.1, 1.1, 1.1];
W.tetas             = [0.340, 0.310, 0.300];
W.tetar             = [0.000, 0.000, 0.000];
W.alfrs             = [0.000, 0.000, 0.000];
W.fi                = [0.000, 0.000, 0.000];
W.alfvg             = [0.120, 0.140, 0.150];
W.en                = [1.120, 1.140, 1.250];
W.alfvg2            = [0.0000, 0.0000, 0.0000];
W.en2               = [0.000, 0.000, 0.000];
W.ifr               = [1, 1, 1];
W.hfc               = -333; % --> check with Antonio!!
W.k0                = [50.00, 20.00, 20.00];
W.k0macr            = [0.000, 0.000, 0.000];
W.bita              = [0.5, 0.5, 0.5];
W.bita2             = [9999, 9999, 9999];
W.ifc               = [1, 1, 1]; 

% ??
W.vpr               = 0.5;
W.tetal             = 0.0002;
W.bital             = 15.0;
W.hsurfmax          = 0.0;

W.itopvar           = 1;
W.ibotvar           = 0;

% iCtopvar              Indice per la lettura dei dati di concentrazione al
%                       contorno superiore:
%                           0:  valore di concentrazione Cinput
%                               (solute_CDE_inp.txt)
%                           1:  condizioni al contorno superiore variabili
%                               (Ctopbound_inp.txt)
W.iCtopvar          = 1;

% itbc                  ???
W.itbc              = 0; % (0=flux; 1=potential)
W.ibbc              = 2; % (0=flux; 1=potential; 2=gradient)

% Can I use W.hqsurf instead of the two following???
W.hsurf             = 9999;
W.qsurf             = 0.1;

W.hbot              = -0.0;
W.qbot              = 9999;
W.grad              = 1.0;

W.inhin             = 0;
% W.inhin     = 1;
% W.hin       = load( fullfile(proj.ipath,'initial_inp.txt') )

W.tmax              = 120;
W.ntp               = 131;
%%   TOP BOUNDARY INPUT
% ----------------------------------
B.top.description   = 'MARWA CORRECTED NOVEMBRE 2011';%--> 'a discretization...'
B.top.ntbc          = 126;
% **DATA**
% times of flux/potential at top boundary.
B.top.thqstar       = 0:1:B.top.ntbc-1;
% flux/potential at top boundary.
B.top.hqstar        = [ -0.961	-0.361	0	-1.325	-0.314	-0.489	-0.489	-0.472	-0.489	0	-0.911	-0.489	-0.489	-0.492	-0.389	-0.417	0	-0.407	0	-0.647	0	-0.833	0	-0.6	0	-0.69	0	-0.685	0	-0.518	0	0	-0.68	0	-0.661	0	-0.638	0	0	-0.578	0	-0.623	0	-0.516	-0.6	0	-0.52	0	-0.451	0	-0.643	0	-0.487	0	0	-0.526	0	-0.463	-0.5	0	-0.726	0	-0.611	0	-0.549	-0.7	0	-0.402	0	-0.598	0	-0.793	0	0	-0.537	0	0	-0.527	0	-0.5	0	-0.575	0	0	-0.499	0	-0.512	0	0	-0.499	0	-0.48	0	-0.508	0	0	-0.562	0	0	0	-0.518	0	0	-0.526	0	0	0	-0.463	0	0	-0.537	0	-0.5	0	-0.401	0	0	-0.544	0	0	0	-0.568	0	0	-0.405	0 ];

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
% Tstar � la temperatura in �C per il calcolo di Kmineralizzazione, sia
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
V.extf              = 0.6;  % esponente legge Beers
V.ifs               = 1;    % indicatore funz. sink {Feddes, vanGen.} 
V.ifg               = 1;    % indicatore funz. distribuz. appar.rad.
V.hI                = -1;   % pot. stress idrico Feddes
V.hII               = 10;   % pot. stress idrico Feddes
V.hIIIH             = -400; % pot. stress idrico Feddes
V.hIIIL             = -600; % pot. stress idrico Feddes
V.hIV               = -8000;% pot. stress idrico Feddes
V.hw50              = -1000;% pot.idrico dimezzamento traspirazione van Genuchten
V.pw1               = 3;    % esponente stress idrico van Genuchten
V.hs50              = -1500;% pot.osmotico dimezzamento traspirazione van Genuchten
V.ps1               = 3;    % esp.stress osmotico van Genuchten
V.aMH               = -760; % par. stress osmotico Mass & Hofmann
V.bMH               = 0.000794;% par. stress osmotico Mass & Hofmann
V.rda               = 1.027;% par. distribuzione radici logistica
V.rdb               = 15.016;% par. distribuzione radici logistica
V.rdc               = 0.074;% par. distribuzione radici logistica
V.zc                = 25;   % par.distribuzione radici doppia-lineare
V.g0                = 0.032;% par.distribuzione radici doppia-lineare
V.gzc               = 0.008;% par.distribuzione radici doppia-lineare
V.Drf               = 85;   % par.distribuzione radici doppia-lineare
%%   SOLUTE Jury INPUT
% ----------------------------------
if W.isol==1
S.J.description     = 'prova puglianello solute transport';
S.J.decad           = 0;
S.J.retard          = 0;
end
%%   SOLUTE CDE INPUT
% ----------------------------------
if W.isol==2
S.CDE.description   = 'MARWA CORRECTED';
S.CDE.tCinput       = 0.000;
S.CDE.tCinput_end   = 0.500;
S.CDE.Cinput        = 0.040;
% Topt:                 The optimum temperature (°C) for the XXX process
S.CDE.Topt          = 25;
% NX:                   Where X = {'H':NH4, 'O':NO3}
S.CDE.NX.Kf1        = [0.100, 1.000]; % ['NH','NO']
S.CDE.NX.Kf2        = [1.000, 1.000];
S.CDE.NX.Kr         = [0.000, 1.000];
end
%%   EC DATA
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
M.nlay              = [22, 33, 11];

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

M.data              = load(fullfile(proj.ipath,'montecarlo.txt'));

% ----Servono questi sotto, e quali?----
M.tetasum           = 0;
M.tetasumSQ         = 0;
M.concsum           = 0;
M.concsumSQ         = 0;
M.fluxsum           = 0;
M.fluxsumSQ         = 0;
% Number of nonconvergences (niter>10,dt<=W.dtmin)
M.nnc               = 0;
% -----------------------------
end