%% PLANT GROWTH

% * Part 1: Crop development 
% IDEV        = 1;      % length of crop cycle: 1 = fixed, 2 = variable
% * If fixed growth length (IDEV = 1), specify:                                                
% LCC         = 168;      % Length of the crop cycle [1..366 days, I]

% * If variable growth length (IDEV = 2), specify:                                                
% TSUMEA = 1050.0       % Temperature sum from emergence to anthesis [0..10000 C, R]
% TSUMAM = 1000.0       % Temperature sum from anthesis to maturity  [0..10000 C, R]

% Tbase:            Start value of temperature sum [°C].
V.Tbase             = 10.0;

% * Part 2: Light Extinction
% ------------------
% kdif:             Extinction coefficient for diffuse visible light [-].
%                   Default is 0.80 and range is [0.40, 1.10].
V.kdif              = 0.80;
% ------------------
% kdir:             Extinction coefficient for direct visible light [-].
%                   Default is 0.90 and range is [0.00, 2.00].
V.kdir              = 0.90;

% * Part 3: Leaf area index or soil cover fraction
% islai:            Choice between:
%                     *1: Leaf Area Index       [implemented]
%                     *0: Soil Cover Fraction   [not implemented]
V.islai             = 1;

% LAI:              Leaf area index (if islai=1) or soil cover fraction (if
%                   islai=0).
%                   *If islai=1, leaf area index [0..12 ha/ha, R]
% * If islai = 0, list soil cover fraction [0..1 m2/m2, R], as function of dev. stage [0..2 -,R]:
% *                     DVS     LAI/SCF
V.LAI               = [
                        0.00    0.05
                        0.30    0.14
                        0.50    0.61
                        0.70    4.10
                        1.00    5.00
                        1.40    5.80
                        2.00    5.20
                       ];

%% IRRIGATION

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

% interval:         Irrigation scheduling interval.
%                       irSTART     irEND
%                       [mm-dd]     [mm-dd]
V.interval          = [ '04-01',    '08-15' ];

% ------------------
% irri:             Table for the definition of GROSS irrigation provided
%                   by farmer [cm d-1]. 
%                   This parameter can assume a value even if the crop is
%                   absent while in SWAP you can have an irrigation
%                   schedule only if the crop is present.
%                   Five parameters of seven must be defined according to
%                   the value assumed by V.flirri (see the '*' below to
%                   understand which one to configure).
%                   It is equal to SWAP-32 "gird", which is gross
%                   irrigation depth.
%                   The following is the list of parameters given in
%                   V.irri:
%                     *Time:    format is 'YYYY-MM-DD,HH'. Valorization
%                     changes accoring to the value of flirri:
%                        1: discrete irrigation dates with iDEPTH volume
%                        2: discrete irrigation dates with h_to threshold
%                        3: continuous dates defining h intervals [from,to]
%                     *iDEPTH:  amount of gross irrigation [cm].
%                     *cirr:    concentration of irrigation water ==> not implemented yet!
%                               Set to ZERO to let the program work
%                     *Type:    switch for type of irrigation
%                                 *0:   surface
%                                 *1:   sprinkler
%                                 *2:   drip        |\
%                                 *3:   fertigation |--> they require depth
%                     *rDEPTH:  depth at which irrigation water is given ==> ask A.Coppola 
%                     *h_from:  pressure head threshold below which
%                               irrigation is performed by the program.
%                     *h_to:    the value of pressure head to which each
%                               soil compartment is restored.
% 
%                   -----------------------------------
V.irri              = {
% flirri-defined
% 1 *               *       *           *       *
% 2 *                       *           *       *                   * 
% 3 #                       *           *       *       *           * 
%   Time            iDEPTH  cirr        Type    rDEPTH  h_from      h_to 
%   [day]           [cm]    [mg cm-3]   [-]     [cm]    [cm]        [cm] 
%   [sdate,edate]   [0,100] [0,1d03]    [0,3]   [0,200] [-1000,0]   [-500,0]  
%   '01-01,12'      NaN     0           NaN     NaN     NaN         NaN
    '01-01,00'      NaN     0           1       0       -500        -20
                      };
%                   -----------------------------------


% * Specify pressure head at field capacity
% * required for timing options  TCS = 2, 3, or 4 and depth option DCS = 1, else dummy 
V.hFC = -100.0;   % soil hydraulic pressure head [-1000.0 .. 0.0,cm, R]







