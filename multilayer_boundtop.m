%% READ VARIABLES

% *CLIMATIC VARIABLES [meteoinput.for]
%   -reference evapotranspiration:
etr             = B.top.eto(    P.tidx )*10;%   [mm d-1]

%   -rainfall [cm d-1]:
%       "grai" [mm] in SWAP-32 is my gauged rain measured by sensor, line:130
%       "grai" is converted from [mm] to [cm], line:161
%       hence I use my "rain" var already expressed in [cm]
rain            = B.top.rain(   P.tidx );%      [cm d-1]

% *CROP VARIABLES
%   -leaf Area Index at current day [m3 m-3]:
lai             = P.LAI(        P.tidx );% need of crop module to compute it
%   -see: http://www.fao.org/docrep/x0490e/x0490e0c.htm#evaporation%20component%20(ke%20eto)
% Kc              = P.Kc(         P.tidx ); % V.Kc should be set as V.Ke!!
% Ke              = P.Ke(         P.tidx );
Ke              = B.top.cfbs;

% *MANAGEMENT VARIABLES
irri            = P.irri(       P.tidx );%      [cm d-1]
% isua:         Switch for sprinkler irrigation:
%                   *1 = sprinkling irrigation,
%                   *0 = surface irrigation.
isSprikler      = P.irTYPE(     P.tidx )==1;%   [-]
%% to be parameterised in conf!!!!!!
% The runon term should be solved using a 2-D [X,Y] grid by means of which
% runoff balance from all pixels of a regional application scheme should be
% calculated as coming out from the simulation step in which IN (rain,
% irrigation, surface soil exfiltration, other) and OUT (soil infiltration,
% evaporation) are both accounted for.
% At now it is set to zero.
runon           = 0;%           [cm]
% as P.runon

% develop!!
melt            = 0;%           [cm]

% ArMpSs:       Area fraction of macropores at soil surface. I set it to
%               zero so that it is not accounted (untill I will implement
%               macropore module.
% ArMpSs          = 0.d0;%        [-]

% SC:           Soil cover fraction. SWAP-32 uses it in place of the the
%               exponential model by Goudriaan (1977) and Belmans (1983)
%               when no LAI is available. I decided not to implement it,
%               because if I don't have LAI I cannot compute wfrac either
%               (and hence I cannot use the equation using SC in which
%               wfrac is also required).

% isflux:       Is flux (and not pressure head) prescribed at the soil
%               surface?
%% init
fl_isRunoff     = false;
%% *WATER at TOP BOUNDARY
% --- determine whether the precipitation falls as rain or snow
%       -fPrecNoSnow:   Ratio rain (excl. snow and rain on snow) / gross
%                       rain flux in case of detailed rainfall data [-].
if W.SwSnow
    % ...to be developed!
else
    fPrecNoSnow = 1;
end

% Amount of rainfall interception during current day [cm] [meteoinput.for]:
%   -check that unit of measurements of aintc is fine!
%   -aintc   --> requires {kdir,kdif,rain,irri,lai}
if lai < 1.d-3 || rain+(irri*isSprikler) < 1.0d-05
% --- neither vegetation nor rainfall/irrigation
        aintc       = 0.d0;
else
% --- canopy interception for agricultural crops:
%     Von Hoyningen-Hune (1983) and Braden (1985), with coefficients:
%     	*a  --> V.avhhb;
%       *b  --> bvhhb
    if (V.avhhb <= 0.000001d0)% [cm]
        aintc       = 0.0d0;
    else
        rpd         = rain*10;
        % if sprinkler, it accounts for interception
        if (isSprikler==1), rpd = (rain + irri)*10; end
        % b coefficient of soil cover fraction:
        %   -exponential relation between soil cover fraction and lai;
        %   -see Eq. 2.53, page, 54, SWAP-32 manual.
        bvhhb       = 1.0d0 - exp(-V.kdif*V.kdir*lai);
        % The calculation min(bvhhb,1.0d0) was found in meteoinput.for, but
        % it is undocumented in SWAP-32 manual:
        bvhhb       = min(bvhhb,1.0d0);
        % See Eq. 2.52, page 54, SWAP-32 manual:
        aintc       = V.avhhb*lai*(1.0d0-(1/(1.0d0+(bvhhb*rpd)/(V.avhhb*lai))))/10;
    end
end

% --- Divide interception into rain part and irrigation part [meteoinput.for]
%       -nraida:    Daily average net precipitation flux [cm d-1]
%       -P.nird:    Net irrigation depth [cm]
if aintc < 0.001d0
    nraida          = rain;% - gsnow - snrai;
    P.nird          = irri;
else
    if isSprikler
        nraida      = rain - aintc*(rain/(rain+irri));
        P.nird      = irri - aintc*(irri/(rain+irri)) ;
    else
        nraida      = rain - aintc;
        P.nird      = irri;
    end
end

% *BARE SOIL? (I defined this "if" condition, check with Antonio!!!)
if ~W.iveg || isnan(V.currcrop)
    flBareSoil      = true;%  soil is bare because no vegetation exists
elseif W.iveg && lai < 1.d-3 % parameterise this threshold?
    flBareSoil      = true;%  soil is bare because vegetation does not exist yet
else
    flBareSoil      = false;% soil is covered by vegetation
end

% Evaporation & Transpiration terms {et0,ew0,es0} [meteoinput.for]:
%   -requires {fao Kc(=swap cf),fao Ke(swap cfbs)}
if flBareSoil
% --- NO CROP
    % Potential transpiration rate, dry crop completely covering the soil [mm d-1]
    et0             = 0;% crop factor is assumed! (and not crop height)
    % Potential transpiration rate, wet crop completely covering the soil [mm d-1]
    ew0             = 0;% crop factor is assumed! (and not crop height)
    % Potential evaporation rate, wet and bare soil [mm d-1]
    es0             = Ke*etr;% convert potential ET into potential E
else
% --- CROP IS PRESENT
    % Potential transpiration rate, dry crop completely covering the soil [mm d-1]
    et0             = P.Kc( P.tidx )*etr;% crop factor is assumed! (and not crop height)
    % Potential transpiration rate, wet crop completely covering the soil [mm d-1]
    ew0             = P.Kc( P.tidx )*etr;% crop factor is assumed! (and not crop height)
    % Potential evaporation rate, wet and bare soil [mm d-1]
    es0             = Ke*etr;% convert potential ET into potential E
end

% Fraction of day the crop is wet [meteoinput.for]:
%   -requires {ew0,aintc}
%   -Eq. 2.62, page 60, SWAP-32 manual
if ew0 < 0.0001d0
    wfrac           = 0.0d0;
else
    % nested max(.,0) and min(.,1) functions force wfrac in range [0, 1]:
    wfrac           = max( min(aintc*10/ew0, 1.0d0), 0.0d0 );
end

% Potential soil evaporation [cm] [meteoinput.for]:
% 	-exponential function by Goudriaan (1977) and Belmans (1983)
%   -requires {es0,lai,wfrac}
%   -Eq. 2.63, page 60, SWAP-32 manual
%   -code in meteoinput.for:332 divides by 10, but I cannot understand why
peva                = max( 0.0d0,  es0*(1.0d0-wfrac)*exp(-V.kdir*V.kdif*lai)/10 );
% --- adapt peva in case of ponding
%   -SWAP-32 uses cfevappond = 1.25d0 [readswap.for:403]
if (pond > 1.0d-10) && (es0 > 1.0d-08)
    peva            = 1.25d0 * peva;
end
% Potential crop transpiration [cm] [meteoinput.for]:
%   -requires {wfrac,et0,peva}
%   -ptra enters sink definition!
P.ptra              = max( 1.0d-9, et0*(1.0d0-wfrac)-peva*10 ) / 10.0;
%                                                    |_mm__|   |_cm_|

% --- ACTUAL RAINFLUX [meteoinput.for:363]
%       -graidt:    Gross precipitation flux during iteration timesteps [cm d-1]
%       -P.nraidt:  Net precipitation flux during iteration timesteps   [cm d-1]
%       -aintcdt:   Interception flux during iteration timesteps        [cm d-1]
if nraida > 1.0d-05
    fInterception   = nraida / rain;
    if aintc < 1.0d-05, fInterception = 1; end
    rainflux        = fPrecNoSnow * rain;
    netrainflux     = fInterception * rainflux;
% --- Set actual gross and net rainflux and interception on TIMESTEP basis
    graidt          = rainflux;
    P.nraidt        = netrainflux;
    aintcdt         = rainflux-netrainflux;
else
    fInterception   = 1;
    graidt          = 0;
    P.nraidt        = 0;
    aintcdt         = 0;
end

% Reduction of potential soil evaporation on daily basis [reduceva.for]
if pond >= 1.0d-10 || B.top.revameth==0
% ---   in case of ponding no reduction ---
	peva            = peva; %#ok<ASGSL>
    P.t_DRY         = 0.0d0;% t_DRY :: is ldwet of SWAP
elseif B.top.revameth==1
% ---   reduction with Black model ---
    % reset t_DRY to zero:
    if (nraida+P.nird) > B.top.Rsigni
        P.t_DRY     = 0.0d0;
    end
    P.t_DRY         = P.t_DRY + 1.0d0;
    peva            = min(peva, B.top.revacoef*(sqrt(P.t_DRY)-sqrt(P.t_DRY-1.0d0)) );
elseif B.top.revameth==2
% ---   reduction with Boesten & Strootsnijder model ---
%   ...to be implemented!!
    error('Reduction of soil evaporation with Boesten & Strootsnijder not implemented yet!!')
else
    error('Parameter "%s=%1d" not properly configured!','B.top.revameth',B.top.revameth)
end

% S O I L   E V A P O R A T I O N  -- [boundtop.for]
%   -moisture fraction (volumetric water content)
% tetatm              = fnteta( hatm, P.sh, 1 ); % NOTE add hysteresis in fnteta
tetatm              = multilayer_fnteta_vgm( hatm, P.sh, 1 ); % NOTE add hysteresis in fnteta
%   -hydraulic conductivity of top boundary
Ksurf               = multilayer_conductivity_node(tetatm,P.sh,1);
if W.SwMacro
    Ksurf           = P.FrArMtrx(1) * Ksurf;
end
%   -hydraulic internodal conductivity at [(top boundary) vs (first
%    compartment)]:
K1atm               = multilayer_conductivity_internode(Ksurf,P.K(1),W.Kmeth,P.nodes.dz(1),P.nodes.dz(1));
%   -Darcian maximum evaporation rate:
%      °see also: Eq. 73, FAO paper 56, Chapter 7, page 144 of pdf
%       ATTENTION: check that the sign of Emax is in accordance with what
%       is required after when you use it (and remember that the code was
%       taken from SWAP).
%      °this multiplies the h gradient by internodal conductivity, which on
%       a general basis computes a flux (Emax is a flux).
Emax                = -K1atm * ( (hatm-P.h(1))/P.nodes.disnod(1) +1 );

% [boundtop.for]:
% reduced soil evaporation rate:
%   -peva is reduced to B.top.redu{
%                                   :1=maximum Darcy flux;
%                                   :2=maximum Darcy+Black;
%                                   :3=maximum Darcy+Boesten/Stroosnijder;
%                                 }
%   -peva used here already accounts for the B.top.redu switch and
%    therefore peva is already reduced accoriding to the selected model.
reva                = min( peva, max(0.0d0, Emax) );

% H I G H   A T M O S P H E R I C   D E M A N D  -- [boundtop.for]
%       -melt:      Melting rate.
%                   [cm d-1]
% I expect that reva is greater than all other terms (q0 is negative) when
% no rain or irrigation are given:
q0                  = (P.nraidt+P.nird+melt) * (1-ArMpSs) + runon -reva;
% the negative q0 is here set as positive, and pond from previous step is
% accounted for in the balance of flux throw ground surface (q1):
q1                  = -q0 - pond_jm1/P.dt;
% NOTE: this of q0 and q1 as the unique following equation:
% q1 = reva - (rain+irri+melt)*(1-ArMpSs) + runon + pond_jm1/dt;

% Check whether the atmospheric demand condition applies:
if q1>=0 && q1>Emax
    isflux          = false;%   a pressure head is prescribed at soil surface!
    hsurf           = hatm;%    pressure head of air is applied to top boundary
    P.Kim2(1)       = K1atm;%   internodal conductivity based on hatm
    pond            = 0;%       ponding is null
    runoff          = 0;%       runoff is null

else
% --- maximum conductivity assuming saturation at ground surface
    % rfcp    : Reduction factor for frozen conditions in each model
    %           compartment [-].
    % hcvsmall: Hydraulic conductivity for complete frozen soils.
    % ksatfit : Saturated hydraulic conductivity (L/T) for each soil layer
    %           (fitted on VG based on lab data) --> P.sh.k0
    % I do not consider frozen conditions in multilayer, therefore I skip
    % the following equation:
    %       Ks = rfcp(1)*ksatfit(1) + (1.0d0-rfcp(1))*hcvsmall;
    % and I use directly my saturated hydraulic conductivity (P.sh.k0) in
    % the calculation of the maximum internodal conductivity:
    K1max           = multilayer_conductivity_internode( P.sh.k0(1), P.K(1), W.Kmeth, P.nodes.dz(1), P.nodes.dz(1) );
    
    % Check whether application of flux=q1 will yield a pressure head >0 at
    % ground surface. If not: flux boundary condition is valid!
    h0              = P.h(1) - P.nodes.disnod(1)*(q1/K1max+1.0d0);
    if h0<1.0d-6%   FLUX controlled
        isflux      = true;
        hsurf       = 0.0d0;
        P.Kim2(1)   = 0.0d0;
        pond        = 0.0d0;
        runoff      = 0.0d0;
        P.qtop      = q1;
    else%           HEAD controlled (ponding occurs!)
        isflux      = false;
        P.Kim2(1)   = K1max;
        fl_isRunoff = true; % runoff potential possible
% --- calculate max value of pond(=h0max) without runoff
        p1          = K1max/P.nodes.disnod(1) * P.dt;
        p2          = 1.0d0/(p1+1.0d0);
        h0max       = p2 * ( pond_jm1 + q0*P.dt - K1max*P.dt + p1*P.h(1) );
% --- in case of macropores, calc. potential overland flow into macrop.: QMpLatSs
        if W.SwMacro
            if h0max > W.pondmax% PndmxMp is used in SWAP-32
                RsRoMp = (h0max + (P.nraidt+P.nird+melt)*ArMpSs*P.dt) / KsMpSs;
                p2Mp = 1 / (p1 + 1 + P.dt*RsRoMp);
                pond = (h0max - W.pondmax) * p2Mp/p2;
                QMpLatSs = min( QMpLatSs, h0max );
                if QMpLatSs < 1.0d-07, QMpLatSs = 0; end
            else
                QMpLatSs = 0;
            end
        end
    end
end
%% clean workspace:
% clear qtop hsurf