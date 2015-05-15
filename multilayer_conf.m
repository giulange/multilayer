%% Documentation about filling this config file
% -----------------------------------------------
% 
% There exist two sets of input data:
%   (1) varying on depth of soil profile (e.g. I_depth.txt);
%   (2) varying on time of simulation (e.g. I_time.txt).
% 
% The program accepts inputs in two ways:
%   (A) they are provided in the config file;
%   (B) they are defined in an external file: one file for (1) and another
%       for (2).
% 
% I pointed out the following schema for data inputs of type (A):
% 
%   (1) This kind of parameters are defined giving the value of the
%       parameter together with the depth at which the value occurs. For
%       instance (given a soil profile defined in Z=[0,300] cm):
%           Z [cm]      W.hin [?]
%        [  0           -100    ;
%           50          -80     ;
%           250         -50     ];
%       I need to define at least one record, the one for the top of the
%       pedon limit (i.e. at Z=0).
%       All other values in between are taken into consideration according
%       to a step function in which the value at current Zi is given by
%       value of the parameter defined at the previous depth.
%       For instance, in the example above W.hin at depth 30 cm is equal to
%       -100, while at depth 75 is equal to -80.
% 
%   (2) The same approach is followed by the kind of variables, except here
%       we refer to the time variable.
%       The Time column must be defined using the "YYYY-MM-DD" format for
%       daily inputs or the "YYYY-MM-DD,HH" for the hourly inputs.
%       Give a look at the example below, assuming that the simulation
%       interval is between 2013-01-01 and 2013-12-31:
%           Time [day]       B.Ctop.Tstar [C]
%       {  '2013-01-01'     18.6
%          '2013-03-01'     23.4
%          '2013-06-01'     25.2
%          '2013-09-01'     20.2
%          '2013-11-01'     19.2    }
%%  Project Info
% ----------------------------------

proj.description    = 'MARWA CORRECTED NOVEMBRE 2011';

% ipath:            The path where all the input files required by
%                   multilayer simulation are in.
proj.ipath          = '/home/giuliano/git/multilayer';

% opath:            The path in which all the output printed files are
%                   stored during multilayer simulation.
%                   A check is performed to ensure that new simulations do
%                   not overwrite old ones, if not needed.
proj.opath          = '/home/giuliano/git/multilayer/output_sol/';

% video:            Set whether to graph (TRUE) the solute transport
%                   simulation in time or not (FALSE).
%                   It is made using the soil grid configuration as
%                   background.
%                   Now it uses the pcolor built-in function, but in future
%                   version I have to enhance it (the nodes are not well
%                   defined, above all the bottom ones!).
proj.video          = false;
%% PROGRAM MODULES:
% -------------------------------------------------------------------------
% MTCL:             Montecarlo simulation
%                   	*0: no
%                       *1:	yes
W.MTCL              = 0;
% ------------------
% wt_mod:           Flag to select the module used to solve the water
%                   transport.
%                       *0:     Basic old-fashioned solver (but working
%                               fine) based on Thomas algorithm.
%                       *1:     New solver based on Newton-Raphson
%                               algorithm (to be comleted and tested).
W.wt_mod            = 1;
% ------------------
% isol:             Flag to set the model of solute transport.
%                   	*0: none solute transport
%                       *1: Jury solute transport
%                       *2: ADE  solute transport (advection-dispersion)
W.isol              = 2;
% ------------------
% iveg:             Flag to set the presence/absence of crop(s) during the
%                   simulation.
%                   	*0: bare soil (no crop)
%                       *1: presence of crop
W.iveg              = 1;
%%   GENERAL | <== from W.*
%%   WATER INPUT
% ----------------------------------

% -----------------------------------------------------
% *RESUME
% (1) SIMULATION, TIME & PRINT
% (2) SOIL GRID GEOMETRY and VERTICAL DISCRETISATION
% (3) HYDRAULIC CHARACTERISTICS of SOIL GRID NODES
% (4) INITIAL CONDITIONS
% (5) TOP BOUNDARY CONDITIONS
% (6) BOTTOM BOUNDARY CONDITIONS
% -----------------------------------------------------

% iosm:             ??
W.iosm              = 0;
% ads:              Adsorbimento soluti
W.ads               = 1;

% -------------------------------------------------------------------------
% (1) SIMULATION, TIME & PRINT:
% -------------------------------------------------------------------------
% On of the most important parameters is W.timestep because its value
% influences the LENGTH of ALL time-dependent parameters and
% variables (with the exception of W.tp, whose length is not affected but
% only its valorization, contrariwise to all other parameters and
% variables).
% ------------------
% sdate:            Start date of simulation run. YYYY-MM-DD
W.sdate             =   '2013-01-01';       % 0
% ------------------
% edate:            End   date of simulation run. YYYY-MM-DD
W.edate             =   '2013-04-30';       % N
% ------------------
% timestep          Time step that is used as reference frame to build all
%                   input data (any other time-dependent parameter).
%                   This means that user has to provide meteo data (and any
%                   other time-dependent external data) using the timestep
%                   here defined.
%                   It can be:
%                       *2.00   --> 2-days;
%                       *1.00   --> one day;
%                       *1/2	--> half a day (12 hours);
%                       *1/24   --> one hour;
%                       *x      --> x timestep.
%                   The total number of timesteps is given by:
%                       N = length( W.sdate:W.timestep:W.edate );
%                   Default is 1.00 and range can be [1/24, 10.00].
%                   At this moment only day can be used (a future version
%                   will be available on a more detailed time step).
W.timestep          = 1;
% ------------------
% dtin:             Initial time step of simulation.
%                   During simulation this value is adapted to accomodate
%                   the valorization of W.timestep (i.e. W.dtin =
%                   W.dtin/W.timestep).
%                   Default value is 1e-5 and range is [].
W.dtin              = 1e-5;
% ------------------
% tolle1:           Maximum difference in pressure head per compartment,
%                   during tridiagonal system solution by Thomas algorithm
%                   (if W.wt_mod==0).
%                   Range is [1e-9 1.0].
W.tolle1            = 1e-4;
% ------------------
% dtmin:            Minimum timestep during simulation.
%                   During simulation this value is adapted to accomodate
%                   the valorization of W.timestep (i.e. W.dtin =
%                   W.dtin/W.timestep).
%                   Default value is 1e-6 and range is [1e-7, 1e-2]
W.dtmin             = 1e-6;
% ------------------
% dtmax:            Maximum timestep during simulation.
%                   During simulation this value is adapted to accomodate
%                   the valorization of W.timestep (i.e. W.dtin =
%                   W.dtin/W.timestep).
%                   The range is [0.01, 0.5]
W.dtmax             = 0.5;
% ------------------
% tp:               Time steps at which print results. A value represents a
%                   fraction of day (e.g. 0.5 is half a day, and so on).
%                   The time step resolution used to build the table of
%                   time-prints is the hour.
%                   During simulation this value is adapted to accomodate
%                   the valorization of W.timestep (i.e. 
%                   W.tp = W.tp*W.timestep).
%                   The time format is 'YYYY-MM-DD,HH'.
% 
%                   -----------------------------------
%                       Time [day]          W.tp [-]
%                       [sdate,edate]       [0,...,N]
W.tp                = {
                        '2013-01-01,00'     1/24
                        '2013-01-02,00'     1.00
                      };
%                   -----------------------------------
% tptole:           This tollerance is used to assign a simulation timestep
%                   to a user defined print-timestep when the first one
%                   does not use the exact W.tp value.
W.tptole            = 1e-02;
% ------------------
% multmin:          ??
W.multmin           = 1.2;
% ------------------
% multmax:          ??
W.multmax           = 0.7;
% 
% THE FOLLOWING PARAMETERS ARE SET IN CASE THE NEWTON-RAPHSON ALGORITHM IS
% SELECTED (if W.wt_mod==1).
% ------------------
% dtfactor:         Multiplicative (divisor) factor by which dt is raised
%                   (reduced) to increment (decrement) convergence speed.
%                   It is applied either when convergence is reached in
%                   less than five Newton-Raphson steps (dtfactor is
%                   multiplicative), or when no convergence can be reached
%                   with current dt (dtfactor is divisor).
%                   Default is 2.0 and range is [1.1,5.0]
W.dtfactor          = 2.0;
% ------------------
% maxit:            Maximum number of steps during the Newton-Raphson
%                   iteration scheme at every simulation timestep.
%                   It is doubled internally when the current timestep is
%                   equal to the minimum user-defined timestep.
%                   Default value is 30 and range is [25, 45].
W.maxit             = 30;
% ------------------
% maxbacktrack:     Maximum number of backtracks during a single
%                   Newton-Raphson iteration. It is used to make a step
%                   size able to reduce 0.5*sum(Fi^2).
%                   [See SWAP-32 manual,page 37]
%                   Default value is 3 and range is [1, 3].
W.maxbacktrack      = 3;
% ------------------
% CritDevBalCp:     Critical water balance deviation for i-th compartment.
%                   It is used to exit the backtrack iteration (which is
%                   implemented to reduce the Newton-Raphson step) when
%                   the F-function at each compartment (=Fi) is below
%                   CritDevBalCp.
%                   The F-function is the difference between:
%                       (i) the derivative of water storage to time, and
%                       (ii)the I/O water terms (head pressure gradient,
%                           internodal conductivity, root, drain and
%                           macropore).
%                   Deafult is 1.0d-6.
W.CritDevBalCp      = 1.0d-06;
% ------------------
% CritDevBalTot:    Critical water balance deviation for the whole soil
%                   profile [cm d-1].
%                   It is used to exit the Newton-Raphson iteration scheme,
%                   which happens when the sum(Fi) is below CritDevBalTot.
%                   Deafult is 1.0d-5.
W.CritDevBalTot     = 1.0d-05;
% ------------------
% CritDevh1Cp:      Convergence criterium for Richards equation: relative
%                   difference in pressure heads [-].
W.CritDevh1Cp       = 1.0d-02;
% ------------------
% CritDevh2Cp:      Convergence criterium for Richards equation: absolute
%                   difference in pressure heads [cm].
W.CritDevh2Cp       = 1.0d-01;
% ------------------
% CritDevPondDt:    Maximum water balance error of ponding layer [cm].
%                   Range is [1.0d-6, 0.1].
W.CritDevPondDt     = 1.0d-04;
% ------------------
% SwkImpl:          Switch for explicit/implicit solution Richards equation
%                   with hydraulic conductivity:
%                       *0 = explicit solution
%                       *1 = implicit solution
W.SwkImpl           = 1;
% ------------------
% SwMacro:          Switch for simulation of macropore flow:
%                       *0: no
%                       *1: yes (not implemented yet)
W.SwMacro           = 0;
% ------------------
% SwSnow:           Switch for simulation of snow accumulation and melt.
%                       *0 = no
%                       *1 = yes (not implemented yet)
W.SwSnow            = 0;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% (2) SOIL GRID GEOMETRY and VERTICAL DISCRETISATION:
% -------------------------------------------------------------------------
% NOTES:
%   sg='soil geometry'
% 
% nlay:             Number of soil layers.
W.nlay              = 3;
% ------------------
% zint:             Soil layer bottom boundaries.
%                   The bottom depth of the lowest layer is the "botlim"
%                   value used to define the bottom Z-limit of all
%                   depth-dependent variables [i.e. botlim = W.zint(end)].
%                   -------------------------------------------------------
% #SoilLay#             1       2       3      4       5       6       7
%                   -------------------------------------------------------
W.zint              = [ 25,     60,     300 ];
% ------------------
% type:             Type of vertical discretization used to build the soil
%                   grid geometry.
%                   The following can be selected:
%                       *1  --> regular grid spacing, according to the 1/2
%                               distance at soil layers saddle points.
%                               Apart from 'nlay' and 'zint' one more
%                               parameter must be defined: the number of
%                               nodes (in W.sg.regular).
%                       *2	--> "sublayers" spacing, in which each layer
%                               can be divided into more layers within
%                               which the spacing is regular. More
%                               parameters must be defined (in
%                               W.sg.sublayers).
%                       *3	--> any other kind of geometry we would like to
%                               implement! (empty at the moment).
W.sg.type           = 1;
% ------------------
% regular:          It creates a soil grid with regular node spacing, at
%                   least within the same soil layer.
%                   Node spacing is quite similar between soil layers, but
%                   a small adjustment is made in order to fit layer
%                   thickness with the putative number of nodes and to make
%                   the bottom layer boundary as the bisector of the two
%                   adjacent nodes.
%                   The following parameter is required:
%                       -numnodes:  The number of nodes for the whole soil
%                                   profile defining the geometry of the
%                                   soil column during simulation of water
%                                   flow. It includes the top and the botom
%                                   boundaries.
%                  %# numnodes #%
W.sg.regular        = 100;
% ------------------
% sublayers:        Define a specific grid in which each sublayer can be
%                   different from the others but within which nodes
%                   spacing is regular.
%                   #CHECK WITH ANTONIO#
%                   You have to define equal nodes spacing across soil
%                   layers boundaries (i.e. the soil layer bottom boundary
%                   must bisect the two adjacent nodes!).
%                   Remeber that the cumulative sum of all hSubLay-s must
%                   equal the bottom depth of the last soil layer, and
%                   that:
%                       hSubLay = hNode x nNodes
% 
%                   User has to define the following parameters:
%                       -SoilLay:   The number of soil layer at which
%                                   one or more sub-layers can be defined.
%                       -SubLay:    The global number of sublayer in soil
%                                   profile.
%                       -hSubLay:   The thickness of current sublayer,
%                                   which can be discretized into one or
%                                   more nodes.
%                       -hNode:     The thickness of each node constituting
%                                   the sublayer.
%                       -nNodes:    The number of regularly spaced nodes
%                                   that discretise the sublayer.
% 
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
% ------------------
% anotherkind:      Whatever we want to implement (we should check what
%                   HYDRUS makes as a suggestion).
W.sg.anotherkind    = [];
% ------------------
% plotme:           Set whether to plot (TRUE) the soil grid configuration
%                   or not (FALSE).
W.sg.plotme         = false;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% (3) HYDRAULIC CHARACTERISTICS of SOIL GRID NODES:
% -------------------------------------------------------------------------
% W.crc.? -->   Curva ritenzione/conducibilità: es. W.crc.dap, W.crc.tetas,
%               ecc.

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% #SoilLay#             1        2       3      4       5       6       7
%--------------------------------------------------------------------------
% dap:              Bulk density [g cm-3]
W.dap               = [ 1.1,     1.1,    1.1     ];
% ------------------
% tetas:            Saturated water content [cm3 cm-3]
W.tetas             = [ 0.340,   0.310,  0.300   ];
% ------------------
% tetar:            Residual water content [cm3 cm-3]
W.tetar             = [ 0.000,   0.000,  0.000   ];
% ------------------
% ifr:              Flag to select the type of retention function:
%                       *1  --> van Genuchten,      unimodal
%                       *2  --> Ross and Smettem,   bimodal
%                       *3  --> Durner,             bimodale
%                       *4  --> Ross and Smettem,   unimodale
W.ifr               = [ 1,       1,      1       ];
% ------------------
% ifc:              Flag to select the type of conductivity function:
%                       *1  --> Mualem unimodal [van Genuchten].
%                       *2  --> Mualem bimodal [Durner].
%                       *3  --> Bimodal, twofold exponential using "fi"
%                               weight + "beta" and "beta2" parameters for
%                               the first and second part of conductivity
%                               function, respectively.
%                       *4  --> Mualem bimodal, [Ross and Smettem for
%                               interconnected pore systems].
%                       *5  --> Mualem bimodal, [Ross and Smettem for
%                               independent system of pores].
%                       *6  --> Mualem unimodal, [Ross and Smettem? when
%                               ifr=4].
W.ifc               = [ 1,       1,      1       ];
% ------------------
% alfrs:            Ross and Smettem alpha empirical shape factor (unused
%                   if retention function is unimodal) [cm-1]
W.alfrs             = [ 0.000,   0.000,  0.000   ];
% ------------------
% fi:               Weight in bimodal retention functions. It must be
%                   always provided, but it is unused in case of unimodal
%                   function.
W.fi                = [ 0.000,   0.000,  0.000   ];
% ------------------
% alfvg:            van Genuchten "alpha" empirical shape factor.
W.alfvg             = [ 0.120,   0.140,  0.150   ];
% ------------------
% en:               van Genuchten "n" empirical shape factor
W.en                = [ 1.120,   1.140,  1.250   ];
% ------------------
% alfvg2:           van Genuchten "alpha" empirical shape factor in case:
%                       *ifr={2,3}, AND
%                       *ifc={2,4,5}.
W.alfvg2            = [ 0.0000,  0.0000, 0.0000  ];
% ------------------
% en2:              van Genuchten "n" empirical shape factor in case:
%                       *ifr={2,3}, AND
%                       *ifc={2,4,5}.
W.en2               = [ 0.000,   0.000,  0.000   ];
% ------------------
% k0:               Saturated conductivity [cm day-1]
W.k0                = [ 50.00,   20.00,  20.00   ];
% ------------------
% k0macr:           Saturated conductivity in macropores [cm day-1].
W.k0macr            = [ 0.000,   0.000,  0.000   ];
% ------------------
% bita:             It is:
%                       *Mualem tau exponent, if ifc={1,2}
%                       *first beta parameter of exponential, if ifc={3}
W.bita              = [ 0.5,     0.5,    0.5     ];
% ------------------
% bita2:            Second beta parameter of exponential, if ifc=3. It must
%                   be always provided, but it is unused in case of
%                   unimodal function.
W.bita2             = [ 9999,    9999,   9999    ];
%--------------------------------------------------------------------------
% #SoilLay#             1        2       3      4       5       6       7
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% ------------------
% Kmeth:            Method for the calculation of the internodal
%                   conductivity. Options are:
%                       *1  --> arithmic mean
%                       *2  --> weighted arithmic mean
%                       *3  --> geometric mean
%                       *4  --> weighted geometric mean
W.Kmeth             = 1;
% ------------------
% hfc:              ??
W.hfc               = -333; % --> check with Antonio!!
% ------------------
% vpr:              Relative humidity (fraction) ?? <-- read from meteo
W.vpr               = 0.85;
% ------------------
% tetal:            Threshold teta for the exponential-type conductivity
%                   function.
%                   Default value is 0.0002
W.tetal             = 0.0002;
% ------------------
% bital:            ??
%                   Default value is 15
W.bital             = 15.0;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% (4) INITIAL CONDITIONS:
% -------------------------------------------------------------------------
% inhin:            Flag to set initial condition:
%                       *0  --> fixed
%                               same potential along whole soil profile
%                       *1  --> variable
%                               potential varies with depth
% W.inhin             = 0; ==> now it is unused (redundant)
% ------------------
% hin:              Initial condition for soil potential at any depth.
%                   -----------------------------------
%                       Depth [cm]      W.hin [?]
%                       [0,botlim]      [-100,+100]
W.hin               = [
                        0               -100
                      ];
%                   -----------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% (5) TOP BOUNDARY CONDITIONS:
% -------------------------------------------------------------------------
% itopvar:          Flag to set top boudary conditions:
%                       *0  --> fixed
%                               The value used is that of hsurf/qsurf
%                               (B.top.hqstar is ignored).
%                       *1  --> variable
%                               B.top.hqstar is used (hsurf/qsurf ignored).
%                   This parameter is now set in multilayer_check_and_load.
% W.itopvar           = 1;
% ------------------
% itbc:             Set INITIAL top boundary condition:
%                       0   --> flux      controlled [qsurf]
%                       1   --> potential controlled [hsurf]
%                   This parameter should be deleted, because it is set
%                   at runtime by the program according to the other
%                   conditions (such as pond, rain, etc.) found there.
W.itbc              = 0;
% ------------------
% Can I use W.hqsurf instead of the two following???
% qsurf:            Top boundary FLUX. Negative in case of infiltation.
%                   If unused its values is NaN.
W.qsurf             = 0.1;
% ------------------
% hsurf:            Top boundary POTENTIAL. Negative in unsaturated
%                   condition.
%                   If unused its values is NaN.
W.hsurf             = NaN;
% ------------------
% iCtopvar          Indice per la lettura dei dati di concentrazione al
%                   contorno superiore:
%                       0:  valore di concentrazione Cinput
%                           (solute_CDE_inp.txt)
%                       1:  condizioni al contorno superiore variabili
%                           (Ctopbound_inp.txt)
W.iCtopvar          = 1;
% ------------------
% hsurfmax:         Max potential allowed at soil surface (ponding at 0
%                   node)
W.hsurfmax          = -0.0;
% ------------------
% pondmax:          Maximum amount of ponding on soil surface before runoff
%                   starts. [cm]
%                   Default value is 0.0 and range is [0.0, 2.0] for well
%                   maintained agricultural fields.
%                   [It should be the same of hsurfmax, which is used in
%                   THOMAS algorithm, and to PndmxMp used in SWAP-32]
W.pondmax           = 1;% [cm]
% ------------------
% rsro:             Drainage resistance for surface runoff. [0.001, 1.0]
%                   See Eq. 4.2, page 71, SWAP-32 manual.
W.rsro              =  0.5;
% ------------------
% rsroExp:          Exponent in drainage equation of surface runoff.
%                   [0.1, 10.0]
%                   See Eq. 4.2, page 71, SWAP-32 manual.
W.rsroExp           =  1.0;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% (6) BOTTOM BOUNDARY CONDITIONS:
% -------------------------------------------------------------------------
% swbotb:           Switch to set the bottom bounday condition in the
%                   Newton-Raphson implementation (this might substitute
%                   ibotvar and ibbc).
%                   This parameter must be set in case W.wt_mod==1.
%                   The following options (should be/)are available:
%                     1  Prescribe groundwater level
%                     2  Prescribe bottom flux (Newmann condition??)
%                     3  Calculate bottom flux from hydraulic head of deep
%                        aquifer (Cauchy condition)
%                     4  Calculate bottom flux as function of groundwater
%                        level (Dirichlet condition??)
%                     5  Prescribe soil water pressure head of bottom
%                        compartment
%                     6  Bottom flux equals zero
%                     7  Free drainage of soil profile
%                     8  Free outflow at soil-air interface
%                   This parameter is used in SWAP-32 too.
W.SwBotB            = 7;
% ------------------
% ibotvar:          Flag to set bottom boudary conditions:
%                       *0  --> fixed
%                               The value used is that of hbot/qbot
%                               (B.bot.hqstar is ignored).
%                       *1  --> variable
%                               B.bot.hqstar is used (hbot/qbot ignored).
W.ibotvar           = 0;
% ------------------
% ibbc:             Kind of bottom boundary condition:
%                       0   --> flux
%                       1   --> potential
%                       2   --> fixed gradient
W.ibbc              = 2;
% ------------------
% qbot:             Bottom boundary POTENTIAL. Negative in unsaturated
%                   condition.
%                   If unused its values is NaN.
W.qbot              = 9999;
% ------------------
% hbot:             Bottom boundary POTENTIAL. Negative in unsaturated
%                   condition.
%                   If unused its values is NaN.
W.hbot              = -0.0;
% ------------------
% grad:             Gradient value assigned to bottom boundary.
W.grad              = 1.0;
% -------------------------------------------------------------------------
%%   CLIMATIC INPUT [B.top]
% ----------------------------------

% **IMPORTANT
% Now relative humidity is fixed in W.vpr! I should read it here in meteo
% file as a dynamic variable.



% NOTE
% We should provide meteo info in input!!
% possible list is
%   -Tmin
%   -Tmax
%   -Rain
%   -ET0

% Now I retrieve climatic data for the sake of developing the multilayer
% simulation model.
% After I have to build a proper I/- interface for meteo input.

% *USING SOILCONSWEB CLIMATIC DATA:
% Select random elements:
iRow        = 1;
iCol        = 1;
% -------------------------------------------------------------------------
% database in MySQL:
db_name     = 'clim_vars_multilayer';
table_res   = 'h24_320m';   % { h24 , h24_320m }

% run the following in case connection does not work:
javaclasspath( '/home/giuliano/.m2/repository/mysql/mysql-connector-java/5.1.17/mysql-connector-java-5.1.17.jar' )

% extract meteo data
query_str   = sprintf( [ ...
                'SELECT date,eto_hs,rad_int,rain_cum,temp_max,temp_min ',...
                'FROM clim_vars_multilayer.%s ',...
                'WHERE row=%d AND col=%d AND ',...
                'date>="%s" AND date<="%s"'
              ], table_res, iRow, iCol, W.sdate, W.edate );
% database name in MySQL:
driver      = 'com.mysql.jdbc.Driver';
db_url      = ['jdbc:mysql://localhost:3306/' db_name];
% set connection:
conn        = database( db_name,'root','pedology-life-2014', ...
                        driver, db_url, ...
                        'Vendor','MySQL', 'Server','localhost' );
% set returning data type:
% setdbprefs('DataReturnFormat','numeric')
%   -to get also the data in DATE format you can use:
setdbprefs('DataReturnFormat','structure')

% get cursor:
curs        = exec( conn, query_str );
curs        = fetch( curs );
meteo_in    = curs.Data;
% -------------------------------------------------------------------------

% assign to multilayer climatic parameters:
B.top.eto       = 0.1*meteo_in.eto_hs(1:end)';%     [cm d-1]
B.top.radi      = meteo_in.rad_int(1:end)';%        []

% ---CORRECT!!!
% B.top.rain      = 0.1*meteo_in.rain_cum(1:end)';%   [cm d-1]
B.top.rain      = -1*[ -0.961	-0.361	0	-1.325	-0.314	-0.489	-0.489	-0.472	-0.489	0	-0.911	-0.489	-0.489	-0.492	-0.389	-0.417	0	-0.407	0	-0.647	0	-0.833	0	-0.6	0	-0.69	0	-0.685	0	-0.518	0	0	-0.68	0	-0.661	0	-0.638	0	0	-0.578	0	-0.623	0	-0.516	-0.6	0	-0.52	0	-0.451	0	-0.643	0	-0.487	0	0	-0.526	0	-0.463	-0.5	0	-0.726	0	-0.611	0	-0.549	-0.7	0	-0.402	0	-0.598	0	-0.793	0	0	-0.537	0	0	-0.527	0	-0.5	0	-0.575	0	0	-0.499	0	-0.512	0	0	-0.499	0	-0.48	0	-0.508	0	0	-0.562	0	0	0	-0.518	0	0	-0.526	0	0	0	-0.463	0	0	-0.537	0	-0.5	0	-0.401	0	0	-0.544	0	0];%	0	-0.568	0	0	-0.405	0 ]';
B.top.rain      = B.top.rain(1:120);
% B.top.rain(4:8:50)      = 3.2;
% ---

B.top.tmax      = meteo_in.temp_max(1:end)';%       [°C]
B.top.tmin      = meteo_in.temp_min(1:end)';%       [°C]
% clean
clear curs conn driver db_url db_name meteo_in driver iCol iRow table_res query_str
%%   TOP BOUNDARY INPUT
% ----------------------------------
B.top.description   = 'MARWA CORRECTED NOVEMBRE 2011';%--> 'a discretization...'
B.top.isirri        = false;
B.top.irri          = 0;                    %   [cm d-1]

% revameth:         Flag to set the method("meth") of reduction("r") of
%                   soil evaporation("eva") on daily basis.
%                       *0  no function - peva is reduced to maximum Darcy
%                           flux;
%                       *1  Black function - peva is reduced to maximum
%                           Darcy flux and to maximum Black (1969).
%                       *2  Boesten/Stroosnijder function - peva is reduced
%                           to maximum Darcy flux and to maximum
%                           Boesten/Stroosnijder (1986).
%                           (not implemented yet)
B.top.revameth      = 0;
% Rsigni:           Minimum rainfall to reset the model of evaporation
%                   reduction according to Black (see SWAP-32 manual page
%                   63), expressed in [cm d-1].
%                   Deafult value is 0.5 and range is [0.0, 10.0].
B.top.Rsigni        = 0.5;% [cm d-1]

% revacoef:         Soil evaporation coefficient of Black or
%                   Boesten/Stroosnijder reduction function [cm d1/2].
%                   Default value is 0.35 and range is [0.00, 1.00].
B.top.revacoef      = 0.54;% grapevine
%%   BOTTOM BOUNDARY INPUT
% ----------------------------------
B.bot.description   = 'MARWA CORRECTED botboundary';%--> 'a discretization...'

% bot.hqstar:       Flux/Potential at bottom boundary.
%                   The time format is 'YYYY-MM-DD,HH'.
% 
%                   -----------------------------------
%                       Time [day]          hqstar [-]
%                       [sdate,edate]       [0,...,N]
B.bot.hqstar        = {
                        '2013-01-01,00'     0.01
%                         '2013-01-02,00'     0.05
%                         '2013-01-02,12'     0.00
%                         '2013-01-04,00'     0.01
%                         '2013-01-05,00'     0.02
%                         '2013-01-05,12'     0.00
%                         '2013-01-06,00'     0.01
%                         '2013-01-07,00'     0.02
%                         '2013-01-08,00'     0.01
%                         '2013-01-08,12'     0.00
                      };
%                   -----------------------------------
%%   SOLUTE TRANSPORT
% -----------------------------------

% NOTE: decide whether they are discrete OR continuous !!

% INPUT di SOLUTI
%   { FR=fertirr, SD=solido, UR=urea, 
%     ORG_rp=organico a mineralizzazione rapida,
%     ORG_sw=organico a mineralizzazione lenta }
% Si assume che FR sia liquido e rappresenti l'apporto in superficie
% (C_input nell'equazione ADE) mentre le altre forme si considerano
% distribuite su uno spessore dL ed entrano nell'ADE come sink-source.
% Tstar è la temperatura in °C per il calcolo di Kmineralizzazione, sia
% rapida che lenta.

% *DEVELOP -----------------------
%  It is set in multilayer_chack_and_load.m as mean air Temperature, but we
%  need an algorithm to derive soil Temperature (i.e. Tstar) from air
%  Temperature!!
% B.Ctop.Tstar        = [ NaN ];
% --------------------------------

% ------------------
% measurement units??
B.Ctop.KhUR         = 1.000; % 
% ------------------
B.Ctop.KvUR         = 1.000;
% ------------------
B.Ctop.KmORG_rp     = 0.020; % andrebbe messo il punto tra "KmORG" ed "rp"
% ------------------
B.Ctop.KmORG_sw     = 0.002; % andrebbe messo il punto tra "KmORG" ed "sw"
% ------------------
% should be defined for each single fertilizer/compound!!!!
B.Ctop.dL           = 30;   % [cm]
% ------------------

% ------------------
% isFert:           This flag sets the use of fertilizers(TRUE) or
%                   compounds(FALSE) below.
S.isfert            = false;
% ------------------
% fertilizers:      The list of fertilizers can be provided as they are
%                   (sold by vendors).
%                   isfert = true
%                   The following formulations are directly available in
%                   the multilayer program:
%                   --organic--       http://h2g2.com/edited_entry/A2339624
%                       *o1    manure
%                          *a   horse (!!write the content!!)
%                          *b   cow
%                          *c   pig
%                          *d   sheep
%                          *e   chicken
%                          *f   rabbit
%                          *g   guano (by seabird/bat)
%                       *o2    sewage
%                       *o3    crop residue    http://en.wikipedia.org/wiki/Crop_residue 
%                          *a   wheat(frumento)
%                          *b   barley(orzo)
%                          *c   oats(avena)
%                          *d   peas(piselli)
%                          *e   
%                       *o4    peat(torba)
%                       *o5    green manure(sovescio)
%                          *a   clover(trifoglio)
%                          *b   vetch(veccia)
%                          *c   soybean(soia)
%                          *d   alfalfa(erba medica)
%                          *e   millet(miglio)
%                          *f   (loietto)
%                          *g   
%                   --mineral--
%                       *m1    straight N
%                          *a   anhydrous ammonium nitrate  [NH4NO3]
%                          *b   urea                        [CO(NH2)2]
%                          *c   calcium ammonium nitrate    [Ca(NO3)2*NH4NO3*10H2O] 
%                       *m2    straight P
%                          *a   fluorapatite                [Ca5(PO4)3F]
%                          *b   hydroxyapatite              [Ca5(PO4)3OH]
%                       *m3    binary NP
%                          *a   monoammonium phosphate      [NH4H2PO4]
%                          *b   diammonium phosphate        [(NH4)2HPO4]
%                       *m4    complex NPK
%                          *    provide directly the X-Y-Z ratings for
%                               N-P-K respectively.
%                       *...
%                   -----------------------------------------------------
%                       Time [day]          type        quantity [kg ha-1]
%                       [sdate,edate]       [0,...10]   [0,...,10]
S.fertilizers       = {
                        '2013-03-01,00'     'm1a'       100.0
                      };
%                   -----------------------------------------------------
% ------------------


% iCtopvar deve poter prevedere impulsi!!

% compounds:        You have to provide the amount of each basic compound
%                   in mg cm-x (the program will convert to g cm-x).
%                   This parameter must be set if S.isfert = false.
S.compounds         = {
%      ----------------------------------------------------------------------------
%      __YYYY-MM-DD,HH__    _____________[mg cm-2]_______________  _[mg cm-3 H2O]_ 
%                           _UREA_  ___organic___   ____solid____   _fertigation_
%       Time                UR      O_rp    O_sw    SD_nh   SD_no   FR_nh   FR_no   
        '2013-01-02,00'     0.000   4.320   25.65   0.000   0.000   0.000   0.000
        '2013-01-14,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
        '2013-01-23,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
        '2013-01-28,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
        '2013-02-02,00'     0.010   0.010   0.010   0.000   0.000   0.000   0.000
        '2013-02-16,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
        '2013-02-23,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
        '2013-03-16,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
        '2013-03-23,00'     0.010   0.010   0.010   0.000   0.000   0.000   0.000
        '2013-04-09,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
        '2013-04-17,00'     0.010   0.010   0.000   0.000   0.000   0.000   0.000
                      };
%      ----------------------------------------------------------------------------

% % % % Cstar.UR:         Fertilization on Top.
% % % %                   The time format is 'YYYY-MM-DD,HH'.
% % % % 
% % % %                   --------------------------------------------
% % % %                       Time [day]          Cstar.UR [g cm-2]
% % % %                       [sdate,edate]       [0,...,N]
% % % B.Ctop.Cstar.UR     = {
% % %                         '2013-01-14,00'     1e-5
% % %                         '2013-01-23,00'     1e-5
% % %                         '2013-01-28,00'     1e-5
% % %                         '2013-02-02,00'     1e-5
% % %                         '2013-02-16,00'     1e-5
% % %                         '2013-02-23,00'     1e-5
% % %                         '2013-03-16,00'     1e-5
% % %                         '2013-03-23,00'     1e-5
% % %                         '2013-04-09,00'     1e-5
% % %                         '2013-04-17,00'     1e-5
% % %                       };
% % % %                   --------------------------------------------
% % % 
% % % % Cstar.ORG.rp:     Fertilization on Top.
% % % %                   The time format is 'YYYY-MM-DD,HH'.
% % % % 
% % % %                   ----------------------------------------------
% % % %                       Time [day]          Cstar.ORG.rp [g cm-2]
% % % %                       [sdate,edate]       [0,...,N]
% % % B.Ctop.Cstar.ORG.rp = {
% % %                         '2013-01-02,00'     4.32e-3
% % %                         '2013-01-14,00'     1e-5
% % %                         '2013-01-23,00'     1e-5
% % %                         '2013-01-28,00'     1e-5
% % %                         '2013-02-02,00'     1e-5
% % %                         '2013-02-16,00'     1e-5
% % %                         '2013-02-23,00'     1e-5
% % %                         '2013-03-16,00'     1e-5
% % %                         '2013-03-23,00'     1e-5
% % %                         '2013-04-09,00'     1e-5
% % %                         '2013-04-17,00'     1e-5
% % %                       };
% % % %                   ----------------------------------------------
% % % 
% % % 
% % % % Cstar.ORG.sw:     Fertilization on Top.
% % % %                   The time format is 'YYYY-MM-DD,HH'.
% % % % 
% % % %                   ----------------------------------------------
% % % %                       Time [day]          Cstar.ORG.sw [g cm-2]
% % % %                       [sdate,edate]       [0,...,N]
% % % B.Ctop.Cstar.ORG.sw = {
% % %                         '2013-01-02,00'     2.565e-2
% % %                         '2013-02-02,00'     1e-5
% % %                         '2013-03-23,00'     1e-5
% % %                       };
% % % %                   ----------------------------------------------
% % % 
% % % 
% % % % Cstar.NH.FR:      Fertigation on top.
% % % %                   The time format is 'YYYY-MM-DD,HH'.
% % % % 
% % % %                   -------------------------------------------------
% % % %                       Time [day]          Cstar.NH.FR [g cm-3 H2O]
% % % %                       [sdate,edate]       [0,...,N]
% % % B.Ctop.Cstar.NH.FR  = {
% % %                         '2013-01-01,00'     0.000
% % %                       };
% % % %                   -------------------------------------------------
% % % 
% % % % Cstar.NO.FR:      Fertigation on top.
% % % %                   The time format is 'YYYY-MM-DD,HH'.
% % % % 
% % % %                   -------------------------------------------------
% % % %                       Time [day]          Cstar.NO.FR [g cm-3 H2O]
% % % %                       [sdate,edate]       [0,...,N]
% % % B.Ctop.Cstar.NO.FR  = {
% % %                         '2013-01-01,00'     0.000
% % %                       };
% % % %                   -------------------------------------------------
% % % 
% % % % Cstar.NH.SD:      Solid fertilization with NH4 on top distributed on a
% % % %                   user defined soil thickness (dL).
% % % %                   The time format is 'YYYY-MM-DD,HH'.
% % % % 
% % % %                   ---------------------------------------------
% % % %                       Time [day]          Cstar.NH.SD [g cm-2]
% % % %                       [sdate,edate]       [0,...,N]
% % % B.Ctop.Cstar.NH.SD  = {
% % %                         '2013-01-01,00'     0.000
% % %                       };
% % % %                   ---------------------------------------------
% % % 
% % % % Cstar.NO.SD:      Solid fertilization with NO3 on top distributed on a
% % % %                   user defined soil thickness (dL).
% % % %                   The time format is 'YYYY-MM-DD,HH'.
% % % % 
% % % %                   ---------------------------------------------
% % % %                       Time [day]          Cstar.NO.SD [g cm-2]
% % % %                       [sdate,edate]       [0,...,N]
% % % B.Ctop.Cstar.NO.SD  = {
% % %                         '2013-01-01,00'     0.000
% % %                       };
% % % %                   ---------------------------------------------
%%   VEGETATION INPUT
% ----------------------------------

% NOTES
% We have to build a crop module, in the meanwhile we use a "simplified"
% version such that present here.
% In a full-developed crop module we should consider all the parameters
% regarding the crop and necessary to run the most complete simulation as
% possible.

V.description       = 'prova Lodi Arm_Art vegetation';%--> 'the plant used was...'
V.nET               = 126;

% Ke:               Potential soil evaporation reducing factor, according
%                   to SWAP-32 crop file and above all FAO paper 56,
%                   Chapter 7, page 142.
%                   Set Ke=1.00 for all simulation period to switch for the
%                   non-use of the Ke parameter.
%                   Default is Ke=1.00, to unuse it.
%                   This parameter should be better implemented: for
%                   instance adding an additional flag parameter we can
%                   switch between the use of FAO approach (to be
%                   implemeted yet) and the use of values provided below by
%                   the user.
% 
%                   -----------------------------------
%                       Time [day]          V.Ke [-]
%                       [sdate,edate]       [0,...,N]
V.Ke                = {
                        '2013-01-01,00'     1.00
%                         '2013-02-01,00'     0.50
                      };
%                   -----------------------------------

% Kc:               It might be computed using FAO paper 56: We should
%                   partitionate Kc in Kbc and Ke and calculate each
%                   coefficient following FAO paper 56, Chapter 7.
% 
%                   -----------------------------------
%                       Time [day]          V.Ke [-]
%                       [sdate,edate]       [0,...,N]
V.Kc                = {
                        '2013-01-01,00'     0.15
                        '2013-03-22,00'     0.60
                      };
%                   -----------------------------------

% *LIGHT EXTINTION
% extf:             Esponente legge Beers.
V.extf              = 0.6;
% kdif:             Extinction coefficient for diffuse visible light [-].
%                   Default is 0.80 and range is [0.40, 1.10].
V.kdif              = 0.80;
% kdir:             Extinction coefficient for direct visible light [-].
%                   Default is 0.90 and range is [0.00, 2.00].
V.kdir              = 0.90;
% avhhb:            Interception coefficient "a" of Von Hoyningen-Hune
%                   (1983) and Braden (1985).
%                   See SWAP-32, Eq.2.52, page 54.
%                   Default value for agricultural crops is 0.025 [cm d-1]
%                   and range is [0.000, 1,000].
%                   Note that you have to use all top boundary water terms
%                   in millimiters for consistency!
V.avhhb             =      0.025;% [cm d-1]

% ifs:              Flag to set the type of sink function {Feddes, vanGen.}
%                       *1      ->
%                       *2      ->
%                       *3      ->
%                       *4      ->
%                       *5      ->
%                       *6      ->
%                       *7      ->
V.ifs               = 1;
% ifg:              FLag to set the type of root distribution function
%                       *1      ->
%                       *2      ->
V.ifg               = 1;

% ***Io farei in questo modo:***
% %                       hI  hII hIIIH   hIIIL   hIV
% hLim                = [ -1  10  -400    -600    -8000   ];
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
% ***** CHANGE IT!! ***** 

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

% ***** CHANGE IT!! ***** 
% %                       rda     rdb     rdc
% V.rd                = [ 1.027   15.016  0.074   ];
% rda:              Par. distribuzione radici logistica
V.rda               = 1.027;
% rdb:              Par. distribuzione radici logistica
V.rdb               = 15.016;
% rdc:              Par. distribuzione radici logistica
V.rdc               = 0.074;
% ***** CHANGE IT!! ***** 

% ***** CHANGE IT!! ***** 
% %                       zc  g0      gzc     Drf
% V.dl                = [ 25  0.032   0.008   85  ];
% zc:               Par. distribuzione radici doppia-lineare
V.zc                = 25;
% g0:               Par. distribuzione radici doppia-lineare
V.g0                = 0.032;
% gzc:              Par. distribuzione radici doppia-lineare
V.gzc               = 0.008;
% Drf:              Par. distribuzione radici doppia-lineare
V.Drf               = 85;
% ***** CHANGE IT!! *****
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