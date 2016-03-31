%% (-1) DEFINITIONS
% DVS:              It stays for development stage.
%                   range is [0, 2]
%                   It is used in place of MMDD dates when the parameter
%                   can be defined independently from hard dates, but same
%                   patterns are possible on sliding dates, and even in
%                   case of shrink/stretch of time intervals.
%                   At least two DVI-s must be defined, with:
%                       -the first DVS=0, start of crop;
%                       -the last  DVS=2, end   of crop.
%% (00) DESCRIPTION
V.description       = { 
    'Chestnut';
};
%% (01) PLANT GROWTH

% -------------------------------------------------------------------------
% * Part 01: Crop development 
% -------------------------------------------------------------------------
% IDEV        = 1;      % length of crop cycle: 1 = fixed, 2 = variable
% * If fixed growth length (IDEV = 1), specify:                                                
% LCC         = 154;      % Length of the crop cycle [1..366 days, I]

% * If variable growth length (IDEV = 2), specify:                                                
% TSUMEA = 1050.0       % Temperature sum from emergence to anthesis [0..10000 C, R]
% TSUMAM = 1000.0       % Temperature sum from anthesis to maturity  [0..10000 C, R]

% Tbase:            Start value of temperature sum [°C].
V.Tbase             = 0.0;% --> might be used to include the "heat unit concept/model"
% ------------------

% -------------------------------------------------------------------------
% * Part 02: Light Extinction
% -------------------------------------------------------------------------
% ------------------
% extf:             Esponente legge Beers.
%                   [Only used in THOMAS algorithm].
V.extf              = 0.6;
% ------------------
% kdif:             Extinction coefficient for diffuse visible light [-].
%                   Default is 0.80 and range is [0.40, 1.10].
V.kdif              = 0.60;
% ------------------
% kdir:             Extinction coefficient for direct visible light [-].
%                   Default is 0.90 and range is [0.00, 2.00].
V.kdir              = 0.75;
% ------------------

% -------------------------------------------------------------------------
% * Part 03: Leaf area index or soil cover fraction
% -------------------------------------------------------------------------
% ------------------
% islai:            Choice between:
%                     *1: Leaf Area Index       [implemented]
%                     *0: Soil Cover Fraction   [not implemented]
V.islai             = 1;
% ------------------
% swLAI:            Switch for type of LAI calculation
%                     *0:   LAI is defined with couples of { DVS, LAI }.
%                           DVS stays for development stage (see
%                           explanation given at the top of file).
%                           At least two couples are necessary to linearly
%                           interpolate middle values.
%                     *1:   #to-be-implemented#
%                           ?N? parameters must be defined in V.LAI:
%                             p1?, p2?, ..., pN?
V.swLAI             = 0;
% ------------------
% LAI:              
% #0 :: swLAI       Leaf area index (if islai=1) or soil cover fraction (if
%                   islai=0) as function of dev. stage [0..2 -].
%                     *If islai=1, list leaf area index [0..12 ha/ha]
%                     *If islai=0, list soil cover fraction [0..1 m2/m2]
%                   At least two {DVI,LAI} couples must be defined, with:
%                       -the first couple with DVS=0,
%                       -the last  couple with DVS=2.
%                   The following list is filled with data managed from
%                   Giovanna's table.
% *                     DVS     LAI/SCF
V.LAI               = [ 
                            NaN
];
                   
% -------------------------------------------------------------------------
% * Part 04: Crop coefficients
% -------------------------------------------------------------------------
% ------------------
% Ke:               I DON'T KNOW IF IT MUST BE DEFINED HERE!!
%                   Potential soil evaporation reducing factor, according
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
%                   See http://www.fao.org/docrep/x0490e/x0490e0c.htm#evaporation%20component%20(ke%20eto)
% 
%                   -----------------------------------
%                       Time [day]          V.Ke [-]
%                       [sdate,edate]       [0,...,N]
V.Ke                = {
%                         '01-01'             1.00
%                         '02-01'             0.50
                      };
%                   -----------------------------------
% ------------------
% Kc:               It might be computed using FAO paper 56: We should
%                   partitionate Kc in Kbc and Ke and calculate each
%                   coefficient following FAO paper 56, Chapter 7.
% 
%                   -----------------------------------
%                       Time [day]          V.Kc [-]
%                       [sdate,edate]       [0,...,N]
V.Kc                = {
                        '01-01'             1.00
                        '05-07'             1.00
                      };
%                   -----------------------------------
% ------------------

% -------------------------------------------------------------------------
% * Part 05: rooting depth, Root density distribution and root growth
% -------------------------------------------------------------------------
%   We should develop a full crop module from which roots are simulated!
% ------------------
% ifg:              FLag to set the type of root distribution function
%                       *1  linear distribution (constant with depth).
%                       *2  logistic distribution.
%                       *3  Prasad-type distribution (triagle).
%                       *4  two-linear distribution.
V.ifg               = 1;
% ------------------
% ifs:              Flag to set the type of sink function {Feddes, vanGen.}
%                   !!splittare tra acqua e soluti con funct distinti!!
%                       *1  FEDDES water reduction factor with uniform root
%                           distribution.
%                       *2  FEDDES water reduction factor with logistic
%                           root distribution.
%                       *3  van Genuchten water reduction factor with
%                           uniform root distribution.
%                       *4  Maas&Hoffman salinity reduction factor with
%                           uniform root distribution.
%                       *5  van Genuchten salinity reduction factor with
%                           uniform root distribution.
%                       *6  van Genuchten multiplicative water and salinity
%                           reduction factor with uniform root
%                           distribution.
%                       *7  Multiplicative FEDDES water and Maas&Hoffman
%                           salinity reduction factors with uniform root
%                           distribution.
V.ifs               = 1;
% ------------------
% swDroot:          Switch for type of root calculation:
%                     *0:   V.Droot is defined with couples of { date, root
%                           depth }. At least two couples are necessary to
%                           linearly interpolate middle values.
%                     *1:   Three parameters must be defined in V.Droot:
%                             L0, Lm, {r OR Tmax}
% 
V.swDroot           = 0;
% ------------------
% Droot:            Parameters to calculate the rooting depth [cm] at every
%                   simulation timestep.
% 
% #0 :: swDroot     The program will linearly interpolate between values
%                   provided in the following table.
%                   If only ONE value is provided, it is assumed that the
%                   same rooting depth is considered all over the
%                   simulation time.
%                   The rooting depth for the first day of simulation must
%                   be provided!
%                   If the latest date is before the end of simulation, it
%                   is assumed that the rooting depth of the latest date
%                   will be constant(fixed) till the end of simulation.
%                   ----------------------------------------
%                       Time [day]          Root Depth [cm]
%                       [sdate,edate]       [0,...,zmax]
V.Droot             = {
                        '01-01'             100.00
                        '05-07'             100.00
                      };
%                   ----------------------------------------
% %                       DVS                 Root Depth [cm]
%                         0.00                5.00
%                         0.30                20.00
%                         0.50                50.00
%                         0.70                80.00
%                         1.00                90.00
%                         2.00                100.00
%  
% #1 :: swDroot     See this paper to include also the "heat unit
%                   concept/model":
%                       "Modeling of Carbon Dioxide Transport and
%                       Production in Soil. 1.Model Development", Simunek &
%                       Suarez, Water Resources Research, vol.29, No.2,
%                       pages 487-497, 1993.
%                   The following three parameters are used to calculate
%                   rooting depth by means of the classical Verhulst-Pearl
%                   logistic growth function (see Eq. 2.22 and 2.21 in
%                   HYDRUS-1D manual, page 20).
% %                   ----------------------------------------
% %                       L0      Lm      r      
% %                       [cm]    [cm]	[T-1]  
% V.Droot             = [ 3.5,    80,     0.04005 ];
% %                   ----------------------------------------
% ------------------
% DrootE:           Root fraction interval considered for irrigation [-].
%                   This parameter has two thresholds between which the
%                   influence of water stress is considered for irrigation
%                   purpose.
%                   For instance if DrootE=[0.25,0.80] and V.flirri=3, it
%                   means that the h_from threshold is considered only for
%                   roots between the fraction 0.25 and the 0.80 (which for
%                   roots of 100 cm means between 25 and 80 cm).
V.DrootE            = [NaN, NaN];

% -------------------------------------------------------------------------
% * Part 06: yield response
% -------------------------------------------------------------------------
% ...nothing at the moment!


% -------------------------------------------------------------------------
% * Part 07: soil water extraction by plant roots
% -------------------------------------------------------------------------
% ***prendere da tabella di SWAP

% ------------------
% %                       hI    hII     hIIIH   hIIIL   hIV
% h_feddes            = [ -15	-30     -325	-600	-8000   ];
% 
% hI:               Feddes Critical pressure head
V.hI                = -10;
% hII:              Feddes Critical pressure head
V.hII               = -25;
% hIIIH:            Feddes Critical pressure head
V.hIIIH             = -1000;
% hIIIL:            Feddes Critical pressure head
V.hIIIL             = -1000;
% hIV:              Feddes Critical pressure head
V.hIV               = -15000;
% ***** CHANGE IT!! ***** 
% ------------------
% %                       hw50    pw1     hs50    sp1
% h_vanGen            = [ -1000   3       -1500   3   ];
% 
% hw50:             Pot. idrico dimezzamento traspirazione van Genuchten
V.hw50              = -1000;
% pw1:              Esponente stress idrico van Genuchten
V.pw1               = 3;
% hs50:             Pot.osmotico dimezzamento traspirazione van Genuchten
V.hs50              = -1500;
% ps1:              Esp.stress osmotico van Genuchten
V.ps1               = 3;
% ------------------
% %                       aMH     bMH
% V.h_mashof          = [ -760    0.000794 ];
% 
% aMH:              parameter stress osmotico Mass & Hofmann
V.aMH               = -760;
% bMH:              parameter stress osmotico Mass & Hofmann
V.bMH               = 0.000794;
% ***** CHANGE IT!! ***** 
% ------------------
% rd:               Parameters for logistic distribution of roots
%                         rda     rdb     rdc
% V.rd                = [ 1.027   15.016  0.074   ];
% rda:              "a" parameter, logistic distribution of roots
V.rda               = 1.027;
% rdb:              "b" parameter, logistic distribution of roots
V.rdb               = 15.016;
% rdc:              "c" parameter, logistic distribution of roots
V.rdc               = 0.074;
% ***** CHANGE IT!! ***** 
% ------------------
% dl:               Parameters for double-linear roots distribution.
%                         zc  g0      gzc     Drf
% V.dl                = [ 25  0.032   0.008   85  ];
% zc:               parameter distribuzione radici doppia-lineare
V.zc                = 25;
% g0:               parameter distribuzione radici doppia-lineare
V.g0                = 0.032;
% gzc:              parameter distribuzione radici doppia-lineare
V.gzc               = 0.008;
% Drf:              parameter distribuzione radici doppia-lineare
V.Drf               = 85;
% ***** CHANGE IT!! *****
% ------------------

% -------------------------------------------------------------------------
% * Part 08: salt stress
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% * Part 09: interception
% -------------------------------------------------------------------------
% ------------------
% avhhb:            Interception coefficient "a" of Von Hoyningen-Hune
%                   (1983) and Braden (1985).
%                   See SWAP-32, Eq.2.52, page 54.
%                   Default value for agricultural crops is 0.025 [cm d-1]
%                   and range is [0.000, 1,000].
%                   Note that you have to use all top boundary water terms
%                   in millimiters for consistency!
V.avhhb             =      0.25;% [mm d-1]
% ------------------

% -------------------------------------------------------------------------
% * Part 10: Root density distribution and root growth
% -------------------------------------------------------------------------
% ...developed in Part 05!
%% (02) IRRIGATION
% -------------------------------------------------------------------------
% (2) IRRIGATION:
% -------------------------------------------------------------------------

% * Part 1: General

% How/whether to implement "irrigazione di soccorso"!

% ------------------
% flirri:           Flag for irrigation application:
%                    *0 --> none:
%                           No irrigation is considered during the
%                           simulation.
%                    *1 --> fixed in both Time & iDEPTH:
%                           It means that the schedule and volumes of
%                           irrigation supplying are known, fixed and
%                           unmodified by the simulation program.
%                    *2 --> fixed in Time (iDEPTH is undefined):
%                           The program will calculate water requirements
%                           by irrigation, hence the schedule is known but
%                           the volume is calculated during the simulation
%                           time.
%                    *3 --> on-demand (both Time & iDEPTH are undefined):
%                           The program builds both the schedule and the
%                           volumes during the simulation.
%                           The date in V.irri must be provided to set the
%                           validity of pressure heads limits ("from" &
%                           "to"): for instance the definition of one row
%                           starting from the first day of simulation means
%                           that a fixed value for each pressure head limit
%                           is set.
V.flirri            = 3;

% * If V.flirri > 0, continue ....

% irrintervals:     ALLOWED irrigation scheduling interval.
%                       irSTART     irEND
%                       [mm-dd]     [mm-dd]
V.irrintervals      = {
                        '04-01',    '08-31'
%                         '05-20',    '08-31'
%                         '04-01',    '04-20'
%                         '05-20',    '08-31'
%                         '04-01',    '04-20'
%                         '05-20',    '08-31'
                      };

% ------------------
% irri:             Table for the definition of GROSS irrigation provided
%                   by farmer [cm d-1]. 
%                   This parameter can assume a value even if the crop is
%                   absent while in SWAP you can have an irrigation
%                   schedule only if the crop is present.
%                   Five/six parameters of seven must be defined according
%                   to the value assumed by V.flirri (see the '*/#' below
%                   to understand which one to configure).
%                   It is equal to SWAP-32 "gird", which is gross
%                   irrigation depth.
%                   The following is the list of parameters given in
%                   V.irri:
%                     *Time:    format is 'YYYY-MM-DD,HH'. Valorization
%                               changes according to the value of flirri:
%                               1: discrete irrigation dates with iDEPTH
%                                  volumes 
%                               2: discrete irrigation dates with h_to
%                                  thresholds
%                               3: continuous dates defining h intervals
%                                  [from,to]
%                     *iDEPTH:  amount of gross irrigation [cm].
%                     *cirr:    NH4+ / NO3- concentration of irrigation
%                               water. 
%                               Set to ZERO to let the program work
%                               properly (indeed it is multiplied by net
%                               irrigation depth!).
%                     *Type:    switch for type of irrigation
%                                 *0:   surface
%                                 *1:   sprinkler
%                                 *2:   drip        |\
%                                 *3:   fertigation |--> they require depth
%                     *rDEPTH:  depth at which irrigation water is given ==> not implemented yet! 
%                     *h_from:  pressure head threshold below which
%                               irrigation is performed by the program.
%                     *h_to:    the value of pressure head to which each
%                               soil compartment is restored.
% 
%                   -----------------------------------
V.irri              = {
% defined on flirri switch
% 1 *               *       *               *       *
% 2 *                       *               *       *                   * 
% 3 #                       *               *       *       *           * 
%   Time            iDEPTH  cirr            Type    rDEPTH  h_from      h_to 
%   [day]           [cm]    [mg cm-3]       [-]     [cm]    [cm]        [cm] 
%   [sdate,edate]   [0,100] [0,1d03]        [0,3]   [0,200] [-1000,0]   [-500,0]  
%                           ---------------
%                           NO3-    NH4+
    '04-01'         NaN     0       0       1       NaN     -400        -300
    '05-01'         NaN     0       0       1       NaN     -500        -300
    '06-01'         NaN     0       0       1       NaN     -700        -300
                      };
%                   -----------------------------------

% * Specify pressure head at field capacity
% * required for timing options  TCS = 2, 3, or 4 and depth option DCS = 1, else dummy 
% V.hFC = -100.0;   % soil hydraulic pressure head [-1000.0 .. 0.0,cm, R]
