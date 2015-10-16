%% DOCUMENTATION
% 
% DESCRIPTION
%   This is the module for the analytical modelling of two forms of N,
%   namely nitrate (NO3-) and ammonia (NH4+).
%   This module accounts for the following processes [see the paper by
%   Liang et al. (2007)]:
%       UREA
%               --> IDROLISI            (pedice 'h')    --> NH4
%               --> VOLATILIZZAZIONE    (pedice 'v')    --> NH3
%       MATERIA ORGANICA
%               --> MINERALIZZAZIONE    (pedice 'm')    --> NH4
% 
%   For nitrification, immobilization and denitrification see Shi et al.
%   (2007).
% 
%   Remember that sl accounts for the type of solute:
%       sl=1    --> NH4+
%       sl=2    --> NO3-
%   while P.tidx accounts for time-dependent inputs on a integer-time
%   length (e.g. day).
% 
% References
%   Gusman & Marino, 1999. Analytical Modeling of Nitrogen Dynamics in soil
%     and groundwater. J.of Irrigation and Drainage Engineering.
%   Cabon et al., 1991. Modelling of the nitrogen cycle in farm land areas.
%     Fertilizer Research 27: 161–169.
%   Sparks, 2003. Environmental Soil Chemistry, 147–186, Second Edition.
%   Shi et al., 2007. An Inverse Method to Estimate the Source-Sink Term in
%     the Nitrate Transport Equation. SSSAJ Soil Sci. Soc. Am. J. 71:26-34.
%   Maggi F. et al., 2008. A mechanistic treatment of the dominant soil
%     nitrogen cycling processes: Model development, testing, and
%     application. Journal of Geophysical Research, VOL. 113, G02016,
%     doi:10.1029/2007JG000578.
%   Liang et al., 2007. Modeling transport and fate of nitrogen from urea
%     applied to a near-trench paddy field. Environmental Pollution 150,
%     313–320.
%   Zhu et al., 2014. Optimizing First-order Rate Coefficients for Soil
%     Nitrate Transformation Processes Applying an Inverse Method. J. Agr.
%     Sci. Tech., Vol. 16: 1173–1185.

%% Variables
% P.cml         :   Solute concentration in mobile region [g cm-3 water]
% csurf         :   Total amount of solutes in ponding layer on soil surface [g cm-2] 
% cfluxb        :   Convective and dispersive flux at bottom border of node [g cm-2]
% cfluxt        :   Convective and dispersive flux at top border of node [g cm-2]
% 
% 
% P.compounds________________________________________________
%            1   2     3     4      5      6      7
% rows --> [ UR; O_rp; O_sw; SD_nh; SD_no; FR_nh; FR_no ]
% cols --> P.tidx
% ___________________________________________________________

%% IMPORTANT NOTES

% LINKS         :   http://ascelibrary.org/doi/abs/10.1061/(ASCE)0733-9437(1999)125:6(330)
%                   
% 
% %****we need a kinetic for NH4+ solid-->liquid solubilization*****
% 
% flZeroCumu    :   Understand when it switches and apply to my code. Also
%                   adjust the initialization of all variables within it!
% cml           :   Add to "cml" the solute quantities coming from
%                   fertilizer applications. 
% cfluxb        :   I think that cfluxb should be changed according to the
%                   value assumed by cml applying the equations for
%                   adsorption (see "iteration procedure for calculation of
%                   cml").
%                   Indeed, when you pass to the next compartment you
%                   cannot have a flux at the top border of node (after the
%                   program set cfluxt=cfluxb) larger than the quantity
%                   effectively available in previous compartment.
% cpond         :   In the code see "solute flux at soil surface". There
%                   you can find an important note on why I think it is not
%                   correct to set as zero cpond and cfluxt in case of
%                   evaporation.
% S.CDE.Cinput  :   I should delete this variable from the whole program
%                   and configuration because solute input is now
%                   implemented in S.compounds & S.fertilizers.
% P.cirr        :   It is now indexed on sl. Is there a sense to its time
%                   dependent valorization?
% CNH3_UR       :   Error in defining this variable, indeed see Liang et
%                   al. (2007), Eq. (8).
% C(t0)         :   Initial mineralizing N concentration in upper-root zone
%                   (0-30 cm) ==> I have to add this parameter to the
%                   config file! Gusman & Marino (1999) used Cr(t0) and
%                   Cs(t0) in Eqs. 13 and 15 to account for slow and rapid
%                   organic-N
% S.CDE.Topt    :   Optimum temperature is usually taken as 35 [°C], while
%                   we set it as 25 [°C]: we should correct it (according
%                   to Cabon et al. 1991 and Shi et al. 2007).
% *print*       :   Select which variables need to be "printed" and save
%                   them in O structure array.

%% DEFs
vsmall      = 1.d-15;
rer         = 1.d-3;
% ---------------------
% TO BE PARAMETERIZED
% REMEMBER TO CONVERT [mg] to [g] BEFORE ENTERING SOLUTE TRANSPORT!!!!!!
DDIF = 0.0;     % Solute molecolar diffusion coefficient in free water [cm2 d-1]
cpre = [0 0];   % Solute concentration in precipitation [mg cm-3]
GAMPAR = 0.0;   % Factor reduction decomposition due to temperature, [0..0.5 °C-1]
RTHETA = 0.3;   % Minimum water content for potential decomposition, [0..0.4 cm3 cm-3]
BEXP   = 0.7;   % Exponent in reduction decomposition due to dryness, [0..2 -]
DECPOT = 0.0;   % Potential decomposition rate, [0..10 d-1, R]
FDEPTH = 1.0;   % List the reduction of pot. decomposition for each soil type [0..1 -]
CREF   = 1.0;   % Reference solute concentration for Freundlich adsorption, [0..1000 mg cm-3]
SWBR   = 0;     % Switch, consider mixed reservoir of saturated zone [Y=1, N=0]
% without mixed reservoir (SWBR = 0), specify:
cdrain = [0 0]; % Mean solute concentration in aquifer or drainage system [g cm-3 water], if SWBR==0 
% In case of mixed reservoir (SWBR = 1), specify:
DAQUIF = NaN;   % Thickness saturated part of aquifer, [0..10000 cm, R]
POROS  = NaN;   % Porosity of aquifer, [0..0.6 -, R]
KFSAT  = NaN;   % Linear adsorption coefficient in aquifer, [0..100 cm3/mg, R]
DECSAT = NaN;   % Decomposition rate in aquifer, [0..10 /d, R]
CDRAINI = NaN;  % Initial solute concentration in groundwater, [0..100 mg/cm3, R]
% ---------------------

% ---------------------
% TO BE INIT in top of the main program:
cseep = cdrain; % Mean solute concentration in upward seepage water at bottom of profile [g cm-3 water] 
% ---------------------

% ---------------------
% TO BE IMPLEMENTED
qdrtot = 0.0d0; % drainage.for:328
% ---------------------

%% init

CNH3_UR                 = 0; %--> andrebbe messa in print (in "O" structure array)
Fmw                     = 0;
Kmr                     = zeros(P.nz,1);
Kms                     = zeros(P.nz,1);
CNH4_UR                 = zeros(P.nz,1);
CNH4_ORG_rp             = zeros(P.nz,1);
CNH4_ORG_sw             = zeros(P.nz,1);
CNH4_SD                 = zeros(P.nz,1);
CNO3_SD                 = zeros(P.nz,1);
CNH4_pn                 = zeros(P.nz,1);
CNO3_pn                 = zeros(P.nz,1);
S1                      = NaN(P.nz,2);
Cfert                   = zeros(2,1);
solbal                  = NaN(1,2);

% ______________________________________
%                      'NH'   'NO'
%       -dispersion length [cm]
LDIS                    = P.CDElambda;
% ______________________________________


% --- reset cumulative solute fluxes
% I'm not sure when this flag turned on/off!
% I feel that when it is false all variables that are equal to zero are
% used anyway without initialization (e.g. when accumulating sqprec!).
% On the contrary, when it is true it seems that I do not have sampro!!
if (flZeroCumu) % **** I have to set this flag in timecontrol
    % See NOTE below on csurf --> I shouldn't put csurf = 0 here!!!!!!
    csurf       = zeros(1,2);
    creactot    = zeros(1,2); % it is the corrispective of SWAP dectot
    sqbot       = zeros(1,2);
    sqsur       = zeros(1,2);
    P.samini    = P.sampro;
    sqprec      = zeros(1,2);
    sqirrig     = zeros(1,2);
    rottot      = zeros(1,2); % to be reviewed when activated!!

    % (NOT IMPLEMENTED YET)
%     sqdra       = 0.0d0;
end

isqbot = zeros(1,2);
isqtop = zeros(1,2);
% isqdra = 0.0d0;         % (NOT IMPLEMENTED YET) 

%% parte adsorbita ***ANTONIO::use cmsy/cml ratio instead!!!
% Dear Antonio, I think we should use the cmsy/cml ratio from previous
% dtwater in order to compute the adsorbate as done in the previous step.
% S1(:,1) = P.cmsy(:,1) ./ P.cml(:,1);
% S1(:,2) = P.cmsy(:,2) ./ P.cml(:,2);

% adsorbate-N by exponential Freundlich:
S1(:,1)             = KF(1)*P.cml(:,1).^FREXP(1);
S1(:,2)             = KF(2)*P.cml(:,2).^FREXP(2);
%% fertigation ***ANTONIO::add this to cpond/cml(1?)***
% Dear Antonio:
% **se lo aggiungiamo a cpond**
%   Dobbiamo necessariamente includere uno nuovo stato al bordo superiore
%   del profilo di suolo, ossia la "deposizione secca" del soluto a formare
%   uno strato in superficie (ad altezza diciamo nulla per semplicità, e
%   con nessun tipo di interazione con l'infiltrazione di acqua, ecc.
%   sempre per semplicità). 
%   Questo perchè per come è fatto adesso, se si verifica evaporazione il
%   soluto non entra e se non lo fa quando lo somministriamo risulta come
%   non pervenuto al suolo.
%   Invece inserendo la deposizione secca, al primo apporto di acqua (irri
%   AND/OR rain) al quale segue un'infiltrazione, assumiamo che la
%   concentrazione de(i/l) solut(i/o) fornita come fertirrigazione fluisce
%   nel primo nodo.
%   QUESTO MECCANISMO DELLA "DEPOSIZIONE SECCA" ANDREBBE COMUNQUE
%   IMPLEMENTATO A PRESCINDERE A CHI ASSEGNAMO IL "Cfert" PROPRIO PER NON
%   FAR SCOMPARIRE IL "cpond" IN UN BUCO NERO.
% **se lo aggiungiamo direttamente al primo nodo**
%   Lo si aggiunge e non vedo al momento prescrizioni particolari. Per il
%   momento opto per questa opzione a minore impatto sul codice.

% AMMONIA (NH4+):
Cfert(1,1)      = P.compounds(6,P.tidx);
% NITRATE (NO3-):
Cfert(2,1)      = P.compounds(7,P.tidx);
%% UREA(hydr,vol} + SOM{min} ***Antonio :: read notes!!***
% NOTES
%   Se parliamo di SOM, dobbiamo necessariamente introdurre un valore di OM
%   come condizione iniziale. Poi si dovrebbe modelizzare come "manure" e
%   "residuals" mineralizzano per formare sia AMMONIO (NH4+, già
%   modellizato) ma credo pure la stessa SOM.
%   Poi dalla mineralizzazione della SOM si produce ancora NH4+.
%   E' come se dovessimo considerare 3 tipo di SOM (vd immagine
%   "SOM-cycle.png" in Dropbox):
%       -POC (1-5y),        --> labile      (P=particulte)
%       -HOC (20-40y),      --> slow        (H=humus)
%       -ROC (500-1000y),   --> stable      (R=resistant)
%****we need a kinetic for NH4+ solid-->liquid solubilization*****

% -------------------------------------------------------------------------
% SOM mineralization :: START
% -------------------------------------------------------------------------

% Nodes receiving fertilizer application:
idL                 = sum(P.nodes.z < B.Ctop.dL);
% Compartmentation of fertilizer application:
Compart             = P.nodes.dz(1:idL) / sum(P.nodes.dz(1:idL));%      [-]

for k=1:idL
    
% --- Calculation of water content factor for mineralization:
%       -It should be computed for each dt-water and used in equation for
%        Km below.
%       -See both Cabon et al. (1991) and Gusman & Marino (1999)
    if P.teta(k)<=P.sh.tetafc(1)
        Fmw         = P.teta(k)/P.sh.tetafc(1);
    else
        Fmw         = P.sh.tetafc(1)/P.teta(k);
    end
    
% --- Rate coefficient for net mineralization :: ƒ(T_abs,teta)
%     ERROR     :   The coefficients Km are valid to compute the daily
%                   rates of mineralization: this means that – as suggested
%                   for UREA-derived N-compounds – we should divide the
%                   daily amounts in all the dt-s in current day. I should
%                   divide the total amount considering the current
%                   dt-water and each compartment height (by means of
%                   Compart).
    % Gusman & Marino (1999), Eq. (6):
    % -rate coefficient for net mineralization of rapidly mineralizing
    %  organic N [dt-1]
    Kmr(k)          = 5.6*10^12 * exp(-9800/(B.Ctop.Tstar(P.tidx)+273))*Fmw;% [dt-1]
    % Gusman & Marino (1999), Eq. (7):
    % -rate coefficient for net mineralization of slowly mineralizing
    %  organic N [dt-1]
    Kms(k)          = 4.0*10^09 * exp(-8400/(B.Ctop.Tstar(P.tidx)+273))*Fmw;% [dt-1]
end

% --- Asynthotic decays of {Urea,ORG-N{rp,sw}} compounds:
%       -ERRORS: 1. Sotto andremo a calcolare l'ammontare dei diversi
%                   composti azotati provenienti dal decadimento
%                   giornaliero di UREA-N e ORG-N. Adesso non è corretto
%                   ripetere questo calcolo per ogni dt-acqua: bisognerebbe
%                   farlo una volta per die (utlizzando il flag
%                   P.flStartOfDay) e frazionare questo ammontare nei
%                   dt-acqua (o dt-soluto?).
%                2. When I apply the Liang et al. Eqs. I should multiply by
%                   "Compart" to divide the amount (which is the total
%                   amount decayed in the upper-root zone).
for i=1:P.tidx% P.kk
    
    % Skip days in which no fertilizers are given:
    if sum(P.compounds(1:3,i)) == 0, continue, end
    
    % Relative time, from today to day of fertilizer application
    tm              = P.time(P.j) -i+1;%floor(P.time(i)); % check with Antonio

% --- (1) UREA :: Liang et al. (2007)
    % Liang et al. (2007), derivative of Eq. (7):
    %   -Daily amount of ammonium produced during urea hydrolysis process
    %    in the upper-root zone (from top of soil profile to B.Ctop.dL cm).
    %   -[g cm-2 day-1] = [g cm-2 day-1] + [g cm-2]*[day-1]*[-]
    CNH4_UR(1:idL)  = CNH4_UR(1:idL) + P.compounds(1,i)*B.Ctop.KhUR*exp(-B.Ctop.KhUR*tm);% [g cm-2 day-1 soil volume]
    % CNH4_UR_CUM(i) = P.compounds(1,i)*(1-exp(-B.Ctop.KhUR*tm));

    % Liang et al. (2007), derivative of Eq. (8):
    %   -Daily amount of ammonia produced during volatilization process
    %    in the upper-root zone (from top of soil profile to B.Ctop.dL cm).
    %   -ERROR: I should use CNH4_UR-partial from previous step instead of
    %           P.compounds(1,i).
    %   -[g cm-2 day-1 soil volume] = [g cm-2 day-1] + [g cm-2]*[day-1]*[-]
    CNH3_UR         = CNH3_UR       + P.compounds(1,i)*B.Ctop.KvUR*exp(-B.Ctop.KvUR*tm);% [g cm-2 day-1 soil volume]
    % CNH3_UR_CUM(i) = P.compounds(1,i)*(1-exp(-B.Ctop.KvUR*tm));      

% --- (2) ORG-N :: Gusman & Marino (1999)
%   -This is computed using the Km coefficient as a:
%     ƒ(daily Temperature, teta in current dt)
%    and the amount of ORG-N derived ammonium is relative to the whole day
%    assuming B.Ctop.Tstar(P.tidx) average temperature and as the current
%    dt-water moisture teta was applied at the whole day.
%    In order to obtain the amount for current dt-water and for each
%    compartment I have to multiply by dt-water and "Compart".

    % Gusman & Marino (1999), Eq. (13):
    % [g cm-2 day-1] = [g cm-2 day-1] + [g cm-2] * [day-1] * [-],     exp(.) has [day]*[day-1] = [-]
    CNH4_ORG_rp         = CNH4_ORG_rp + P.compounds(2,i).*Kmr.*exp(-Kmr.*tm);% [g cm-2 day-1 soil volume]
    
    % Gusman & Marino (1999), Eq. (14):
    % [g cm-2 day-1] = [g cm-2 day-1] + [g cm-2] * [day-1] * [-],     exp(.) has [day]*[day-1] = [-]
    CNH4_ORG_sw         = CNH4_ORG_sw + P.compounds(3,i).*Kms.*exp(-Kms.*tm);% [g cm-2 day-1 soil volume]
end

% --- (3) INORG-N :: solid NH4+ & NO3-
% L'apporto della forma solida NH4 ed NO3 dura per l'intero periodo
% P.tidx:P.tidx+1. Essendo in [g cm-2 dt-1], l'input è già nella forma
% richiesta dall'ADE.
% Questo viene poi moltiplicato per P.dt nella stessa equazione.
% Una volta costituito, il POOL si assume distribuito per l'intero spessore
% di suolo B.Ctop.dL ==> lo si divide per B.Ctop.dL e si ottiene una
% concentrazione in g/cm3 di suolo.
% ***ERRORS:    -1- La ripartizione deve essere proporzionale allo spessore 
%                   di ciascun nodo!! ==> FATTO!!
%               -2- Manca il meccanismo di solubilizzazione in acqua!
%                   Mentre adesso stiamo assumendo che tutto l'N fornito
%                   come ammonio e nitrato solido è istantaneamente
%                   disponibile in acqua!!
%               -3- Si potrebbe implementare un discioglimento in acqua
%                   frazionato e graduale nel tempo
%                   (passando per una variabile di stato).
%               -4- A me sembra che la quantità P.compounds(4:5,P.tidx)
%                   viene somministrata N volte per tutti i dt che si
%                   vengono a creare nel giorno di concimazione! O
%                   ripartiamo nella giornata sui dt, oppure tutto viene
%                   dato al primo dt (assunzione della somministrazione
%                   istantanea) e non più a tutti gli altri dt del giorno
%                   corrente!!! Sarei per moltiplicare i composti azotati
%                   solidi ed inorganici somministrati per (dt/day), che
%                   sarebbe proprio uguale a dt/1, ossia adimensionale
%                   [dt * dt-1] = [-].
%                   D'altronde serve anche ai fini di una correttezza
%                   formale delle unità di misura per entrare in CDE.
%               -5- L'apporto solido è definito in [g cm-2] nel file
%                   config, mentre qui stando accorti si potrebbe optare
%                   per un [g cm-2 day-1], cioè si assume che l'apporto
%                   solido non è un tasso istantaneo ma giornaliero
%                   (significa che il solido è addittivo alle frazioni da
%                   UREA e ORG-N) e quindi va ripartito per i dt-water che
%                   si vengono a creare moltiplicando per il fattore 
%                   dt/day = dt/1 = dt con units [-].
% NH4+ from solid fertilization:
CNH4_SD(1:idL)  = P.compounds(4,P.tidx) * Compart;%                 [g cm-2 day-1 soil volume]
% NO3- from solid fertilization:
CNO3_SD(1:idL)  = P.compounds(5,P.tidx) * Compart;%                 [g cm-2 day-1 soil volume]

% --- (4) ORG-N :: Build the pool, it should be [g cm-2 day-1]
% [g cm-2 day-1] = [g cm-2 day-1] + [g cm-2 day-1] + [g cm-2 day-1] + [g cm-2 day-1]
CNH4_pn         = CNH4_UR + CNH4_ORG_rp + CNH4_ORG_sw + CNH4_SD;%   [g cm-2 day-1 soil volume]
CNO3_pn         = CNO3_SD;%                                         [g cm-2 day-1 soil volume]

% The pools are in [g cm-2 day-1] and enter the SST which is expressed in
% units [g cm-3 dt-1]. This means that when we multiply by dt in the
% conservation equation for the substance with get [g cm-3 dt-1] as we
% need!
% We also need to divide by theta in order to get (soil volume) ––> (soil
% water), and we have to decide where to apply this transformation (or when
% computing creact or in place in the conservation equation).
CNH4_pn         = CNH4_pn ./ P.nodes.dz(1:P.nz);% * (P.dt/1);%        [g cm-3 day-1]
CNO3_pn         = CNO3_pn ./ P.nodes.dz(1:P.nz);% * (P.dt/1);%        [g cm-3 day-1]

% -----------------------------------------------------------------------
% SOM mineralization :: END
% -----------------------------------------------------------------------
%% COUPLED TRANSPORT/REACTIVE ANALYTICAL MODELLING (i.e. CDE + SST)
% See Zhu et al. (2014):
%   CDE   :   Convection-dispersion equation
%   SST   :   Source-sink term

% -------------------------------------------------------------------------
% (CDE + SST) :: START
% -------------------------------------------------------------------------

% ___________________________________________________________________________________________ 
%  D T – S O L U T E
% ___________________________________________________________________________________________         
% RATIO:
% Compared to an implicit, iterative scheme, the explicit scheme in Eq.
% 8.15 has the advantage that incorporation of non-linear adsorption,
% mobile/immobile concepts, and other non-linear processes is relatively
% easy.
% In order to ensure stability of the explicit scheme, the time step
% dtsolu_j should meet the criterium (Van Genuchten and Wierenga, 1974):
%   dtsolu(j) <= ( dz(i)^2 * theta(i)(j) ) / 2*D(i)(j)      (Eq. 8.16)
% 
% % % % IMPLEMENTATION:
% % % % --- determine maximum timestep [computed on nodes (1:nz-1)]:
% % % thetav = P.nodes.inpola(2:P.nz).*P.teta(1:P.nz-1) + ...
% % %          P.nodes.inpolb(1:P.nz-1).*P.teta(2:P.nz);
% % % %       -diffusion  coefficient [cm2 d-1], Eq. 8.2, SWAP manual
% % % %       -tourtuosity factor in the water phase according to Millington and
% % % %        Quirk [1961], Eq. 3.56, HYDRUS1D manual
% % % diffus = DDIF * (thetav.^2.33)./(P.sh.tetas(1:P.nz-1).^2);
% % % %       -dispersion coefficient [cm2 d-1], Eq. 8.5
% % % disper = LDIS(1:P.nz-1).*abs(P.q(1:P.nz-1))./P.teta(1:P.nz-1);
% % % %       -overall (diffusion + dispersion) coefficient [cm2 d-1]
% % % Dfs    = diffus + disper;
% % % Dfs(Dfs<1.0d-8) = 1.0d-8;
% % % %       -criterium on solute dt, Eq. 8.16
% % % dummy  = P.nodes.dz(1:P.nz-1).^2 .* P.teta(1:P.nz-1) /2.0 ./Dfs;
% % % dtsolu = min([P.dt;dummy]);

dtsolu = P.dt;

% === List of Variables accumulated in while ==============
%    {  tcumsum, csurf(sl), dectot(sl), rottot, cdrtot, ...
%       isqdra, sqdra, P.cmsy(inode,sl), cdrain(), sqsur(), ...
%       sqbot(sl), 
%    }
% =========================================================
tcumsol = 0;
while P.dt-tcumsol >= 1.0d-08% I should consider a threshold accounting dtmin
    % Counts the iterations number required by the solute transport module:
    P.Ndtsolute(P.tidx) = P.Ndtsolute(P.tidx) +1;
    
% ---    time step and cumulative time
    dtsolu      = min(dtsolu,(P.dt-tcumsol)); % take dt if dtsolu is greater
    dtsolu      = max(dtsolu,W.dtmin);        % take dtmin if dtsolu is lesser
    tcumsol     = tcumsol + dtsolu;
        
    for sl=1:2
% ___________________________________________________________________________________________ 
%  C D E   F L U X   A T   T O P   B O R D E R   O F   F I R S T   N O D E
% ___________________________________________________________________________________________              
% --- solute flux at soil surface
        % NOTE: Here there should be a term accounting for existing ponding
        %       even though the flZeroCumu switch zeroed csurf. As a matter
        %       of fact you can have a new day with ponding coming from
        %       previous day, and therein have a concentration which would
        %       be accounted for by csurf together with upcoming rain and
        %       irrigation volumes in dt-solute. Furthermore, cpond cannot
        %       become zero in case of evaporation, because I don't know if
        %       current evaporation consumes the whole ponding, therefore
        %       the only place in which I can set cpond as zero is in top
        %       boundary module in which I know when pond is present or
        %       not and its height!
        % csurf : Total amount of solutes in ponding layer on soil surface [g cm-2] 
        csurf(sl)   = csurf(sl) + (P.nird*P.cirr(sl,P.tidx) + ...
                         P.nraidt*cpre(sl))*dtsolu;         % [g cm-2]
        if (P.qtop < -1.d-6)% infiltration
            % cpond : Mean solute concentration in ponding layer on soil surface [g cm-3] 
            cpond   = csurf(sl) / (pond-P.qtop*dtsolu);     % [g cm-3]
            % convective and dispersive flux at top border of first node
            cfluxt  = P.qtop*(1.0d0-ArMpSs)*cpond*dtsolu;   % [g cm-2]
            csurf(sl)= csurf(sl) + cfluxt;                  % [g cm-2]
            % isqtop : Solute flux through the soil top surface [g cm-2 d-1] 
            isqtop(sl) = P.qtop*(1.0d0-ArMpSs)*cpond;       % [g cm-2 d-1]
        else % evaporation
            cpond   = 0.0d0;% to be deleted!?               % [g cm-3]
            cfluxt  = 0.0d0;                                % [g cm-2]
        end
        
% --- calculate mass balance for each compartment
        for inode = 1:P.nz

% ___________________________________________________________________________________________ 
%  C D E   F L U X   A T   B O T T O M   B O R D E R   O F   N O D E S
% ___________________________________________________________________________________________         
% --- convective and dispersive fluxes [g cm-2]
            if (inode < P.nz)% all nodes except the bottom one
                cmlav   = P.nodes.inpola(inode+1) * P.cml(inode,sl) + P.nodes.inpolb(inode) * P.cml(inode+1);% [g cm-3 water]
                thetav  = P.nodes.inpola(inode+1) * P.teta(inode) + P.nodes.inpolb(inode) * P.teta(inode+1);%  [-]
                %       -pore water velocity v = q/theta (Bolt, 1979)
                vpore   = abs(P.q(inode+1))/thetav;% [cm dt-1]
                %       -diffusion  coefficient [cm2 d-1], Eq. 8.2
                diffus  = DDIF*(thetav^2.33d0)/(P.sh.tetas(inode)^2.0d0);% [cm2 dt-1]
                %       -dispersion coefficient [cm2 d-1], Eq. 8.5
                disper  = LDIS(inode)*vpore; %+ 0.5d0*dtsolu*vpore^2;% [cm2 dt-1]
                %       -overall (diffusion + dispersion) coefficient [cm2 dt-1]
                Dfs     = diffus + disper;% [cm2 dt-1]
                % convective and dispersive flux at bottom border of inode:
                %  This is the piece with all i+1/2 terms of Eq. 8.15:
                %    (q*c)_ip1 + (theta*D*deltaC/disnod)_ip1
                % [g cm-2] = ([cm dt-1] * [g cm-3 water] + [-]*[cm2 dt-1]*[g cm-3 water]/[cm]) * [dt] 
                cfluxb  = ( P.q(inode+1) * cmlav + ... % convective
                            thetav * Dfs * (P.cml(inode+1,sl)-P.cml(inode,sl)) / P.nodes.disnod(inode+1) ... % dispersive
                          )*dtsolu;% [g cm-2]
            else% bottom node
                if (P.q(inode+1) > 0.0d0)
                    cfluxb = P.q(inode+1)*cseep*dtsolu;%           [g cm-2]
                else
                    cfluxb = P.q(inode+1)*P.cml(inode,sl)*dtsolu;% [g cm-2]
                end
            end

% ___________________________________________________________________________________________ 
%  S O L U T E   R E A C T I O N   B Y   B I O L O G I C A L   T R A N S F O R M A T I O N S
% ___________________________________________________________________________________________ 
% 
%     Scheme of nitrogen transformation processes:
%  _____________________________________________________
% |                                                     |
% |      volatilization    plant uptake      volat.     |
% |                 ^       ^        ^         ^        |
% |                 |____   |        |         |        |
% |                      \  |        |         |        |
% |          (NH4+) <–––> NH4+ –––> NO3- –––> N2O       |
% |           ads          ^                            |
% |                        |                            |
% |                        |                            |
% |                       SOM                           |
% |_____________________________________________________|
% 
%     See Cabon et al. (1991), Shi et al. (2007) and Zhu et al. (2014) for
%     equations and parameters:
%       -nitrification      :   NH4+ ––> NO3-   (Cabon et al., 1991)
%           -lq             :   //
%           -sd             :   ??
%       -immobilization     :   NO3- ––> N-org  (Cabon et al., 1991)
%       -denitrification    :   NO3- ––> N2O    (Lafolie, 1991; MeGechan and Wu, 2001)
%       -root uptake        :   Zhu et al. (2014) which refers to Schoups
%                               and Hopmans (2002).
            if P.teta(inode)<P.sh.tetafc(inode)
                teta_ratio  = P.teta(inode)/P.sh.tetafc(inode);
            else
                teta_ratio	= P.sh.tetafc(inode)/P.teta(inode);
            end
% --- nitrification     [g cm-3 dt-1 soil volume]
            % Sn(z,t), Eq. 15 in Shi et al. 2007 [g cm-3 dt-1 soil volume]
            Cntf_lq         = P.CDEKnitr(inode)*1.071^ ...
              (B.Ctop.Tstar(P.tidx)-S.CDE.Topt)*teta_ratio*P.cml(inode,1)*P.teta(inode);% [g cm-3 dt-1 soil volume]
            % Reference ???
            %   -Check that this amount calculation is correct assuming S1
            %    adsorbed.
            Cntf_sd         = P.CDEKnitr(inode)*1.071^ ...
              (B.Ctop.Tstar(P.tidx)-S.CDE.Topt)*teta_ratio*P.sh.dap(inode)*S1(inode,1);% [g cm-3 dt-1 soil volume]
            if sl==1
% --- SST(NH4+)         [g cm-3 dt-1 soil volume]
                % solute conc. from equilibrium reaction [g cm-3 dt-1 soil volume]
                creact      = -Cntf_lq -Cntf_sd +CNH4_pn(inode);% [g cm-3 dt-1 soil volume]
            elseif sl==2
% --- immobilization    [g cm-3 dt-1 soil volume]
                % Sm(z,t), Eq. 16 in Shi et al. 2007 [g cm-3 soil volume]
                Cimb        = P.CDEKimmob(inode)*1.050^ ...
                  (B.Ctop.Tstar(P.tidx)-S.CDE.Topt)*teta_ratio*P.cml(inode,2)*P.teta(inode);% [g cm-3 soil volume]
% --- denitrification   [g cm-3 dt-1 soil volume]
                % teta_tsh: the threshold water content for de-nitrification:
                teta_tsh    = 0.627*P.sh.tetafc(inode)-0.0267*(P.sh.tetas(inode) ...
                             -P.teta(inode))/P.sh.tetas(inode)*P.sh.tetafc(inode);
                if P.teta(inode)>teta_tsh
                    % Sd(z,t), Eq. 17 in Shi et al. 2007, [g cm-3 dt-1 soil volume]
                    Cden	= P.CDEKdenitr(inode)*1.07^(B.Ctop.Tstar(P.tidx)- ...
                              S.CDE.Topt)*(P.teta(inode)-teta_tsh) / ...
                              (P.sh.tetafc(inode)-teta_tsh)*P.cml(inode,2)*P.teta(inode);% [g cm-3 dt-1 soil volume]
                else
                    Cden    = 0;
                end
% --- SST(NO3-)         [g cm-3 soil volume]
                % solute conc. from equilibrium reaction [g cm-3 soil volume]
                %   fprintf('Cntf_lq=%g,  Cntf_sd=%g,  Cimb=%g,  Cden=%g,  CNO3_pn=%g\n', Cntf_lq,Cntf_sd,Cimb,Cden,CNO3_pn(inode) )
                creact      = +Cntf_lq +Cntf_sd -Cimb -Cden +CNO3_pn(inode);% [g cm-3 dt-1 soil volume]
%             C2_chem         = ( Creact*P.dt + P.cml(inode,sl)*P.teta(inode) )/P.teta(inode);
            end
% --- cumulative SST [g cm-2 soil volume]
            % creactot  :  Cumulative amount of solute after chemical and
            %              biological reactions (same as SWAP dectot)
            creactot(sl) = creactot(sl) + creact*dtsolu*P.nodes.dz(inode);% [g cm-2] %% ANTONIO??
%            dectot(sl) = dectot(sl) + ctrans*dtsolu*P.nodes.dz(inode);

% ___________________________________________________________________________________________ 
%  S E L E C T I V E   S O L U T E   U P T A K E   B Y   P L A N T   R O O T S
% ___________________________________________________________________________________________         
% --- selective solute uptake by plant roots
            % crot : selective (NH4+ or NO3- only) solute uptake by plant roots [g cm-3 dt-1] 
            % tscf : Relative uptake of solutes by roots [-] 
            % qrot : Array with root water extraction flux for each compartment [cm dt-1] 
            % cml  : Solute concentration in mobile region [g cm-3 water]
            % [g cm-3 dt-1] = [-] * [cm dt-1] * [g cm-3 water] / [cm]
            crot        = S.CDE.NX.Kr(sl)*P.sink(inode)*P.cml(inode,sl)/P.nodes.dz(inode);% ANTONIO: add "/dz" in your code!!
            %crot = tscf*qrot(inode)*P.cml(inode)/dz(inode);
            % [g cm-2] = [-] * [cm dt-1] * [g cm-3 water] * [dt]
            rottot(sl)  = rottot(sl) + S.CDE.NX.Kr(sl)*P.sink(inode)*P.cml(inode,sl)*dtsolu;
            %rottot = rottot + tscf*qrot(inode)*P.cml(inode)*dtsolu;

% ___________________________________________________________________________________________ 
%  S O L U T E   L A T E R A L   D R A I N A G E   F L U X
% ___________________________________________________________________________________________                     
% --- lateral drainage (NOT IMPLEMENTED YET)
            % qdra      : Array with lateral drainage flux [g d-1] for each drainage level and compartment
            % nrlevs    : Number of drainage levels
            cdrtot(sl)  = 0.0d0;
%             for level = 1 : nrlevs
%                if (qdra(level,inode) > 0.0d0)
%                   cdrtot = cdrtot+qdra(level,inode)*P.cml(inode)/dz(inode);
%                else
%                   cdrtot = cdrtot+qdra(level,inode)*cdrain/dz(inode);
%                end
%             end
% --- cumulative amount of solutes to lateral drainage (NOT IMPLEMENTED YET) 
%             isqdra = isqdra + cdrtot*dz(inode)*dtsolu;
%             sqdra = sqdra + cdrtot*dz(inode)*dtsolu;

% ___________________________________________________________________________________________ 
%  C O N S E R V A T I O N   E Q U A T I O N   F O R   T H E   S U B S T A N C E 
% ___________________________________________________________________________________________         
% --- (CDE+SST) :: conservation equation for the substance
            % This considers all the terms with i+1/2 in Eq. 8.15, SWAP manual: 
            %   -I substituted ctrans with creact that accounts for more
            %    transformations in multilayer than in SWAP.
            %   -cmsy  :  dissolved + adsorbed solute concentration in
            %             mobile region [g cm-3 soil volume]
            %   -cflux*:  convective and dispersive flux borders of inode [g cm-2]
            % [g cm-3 soil volume] = [g cm-3 soil volume] + ([g cm-2]) / [cm] + ([g cm-3 dt-1 soil volume])*[dt-1] 
            P.cmsy(inode,sl) = P.cmsy(inode,sl) + (cfluxb-cfluxt) / P.nodes.dz(inode) + ...
                    (+creact-crot-cdrtot(sl)) * dtsolu;% [g cm-3 soil volume]

% ___________________________________________________________________________________________ 
%  S O L U T E   C O N C E N T R A T I O N   I N   M O B I L E   R E G I O N 
% ___________________________________________________________________________________________         
% --- iteration procedure for calculation of cml
%  ***ERROR:
%     I think that cfluxb should be changed according to the value assumed
%     by cml applying the following equations. Indeed, when you pass to the
%     next compartment you cannot have a flux at the top border of node
%     (after the program set cfluxt=cfluxb) larger than the quantity
%     effectively available in previous compartment. Maybe this reasoning
%     is valid only if cmsy < vsmall.
%     In any case I think that we use a smaller dtsolute than dtwater in
%     order to manage this inconsistency.
            if (P.cmsy(inode,sl) < vsmall)
               P.cmsy(inode,sl) = 0.0d0;
               P.cml(inode,sl)  = 0.0d0;
               % ...I should apply an update of both the "flux" (i.e. CDE)
               % and the "react" (i.e. SST) terms in order to consider
               % physically possible amounts.
               % -SST term- At first step I assume that the substance
               %            undergoes to trasformations (and I should check
               %            whether this is possible or it would consume
               %            more than what is available)
               % ...               
               % -CDE term- At second step, I subtract the transformations
               %            from the substance amount and then I leave the
               %            substance residual to the cfluxb, which is
               %            suddenly used below to update the cfluxt of the
               %            next node.
               % ...
            else%    linear      Freundlich
                if (abs(FREXP-1.0d0) < 0.001d0)
                    % [g cm-3 water] = [g cm-3 soil volume] / ([cm3 water cm-3 soil volume]+[g cm-3]*[cm3 g-1]) 
                    P.cml(inode,sl) = P.cmsy(inode,sl) / (P.teta(inode)+P.sh.dap(inode)*KF(sl));
                else% exponential Freundlich
                    if (P.cml(inode,sl) < vsmall), P.cml(inode,sl) = vsmall; end
                    differ = true;
                    while (differ)
                        old = P.cml(inode,sl);
                        % [??] = [g cm-3]*[cm3 g-1] * ([g cm-3 water]/[g cm-3])^[–]
                        dummy = P.sh.dap(inode)*KF(sl)*(P.cml(inode,sl)/CREF)^(FREXP-1.0d0);
                        % [g cm-3 water] = [g cm-3 soil volume] / [cm3 water cm-3 soil volume] 
                        P.cml(inode,sl) = P.cmsy(inode,sl)/(P.teta(inode)+dummy);
                        if (abs(P.cml(inode,sl)-old) < rer*P.cml(inode,sl)), differ = false; end
                    end
                end
            end

% ___________________________________________________________________________________________ 
%  C D E   F L U X   A T   T O P   B O R D E R   O F   N O D E S
% ___________________________________________________________________________________________              
% --- make top flux next compartment equal to current bottom flux
            cfluxt = cfluxb;

% --- next 'compartment'
        end% node

% ___________________________________________________________________________________________ 
%  T O   B E   I M P L E M E N T E D
% ___________________________________________________________________________________________              
% --- solute balance in aquifer for breakthrough curve :: TO BE IMPLEMENTED YET 
         if (SWBR == 1) 
            if (qdrtot > 0.0d0)
               cdrain(sl) = cdrain(sl) + dtsolu/(POROS + P.sh.dap(inode)*KFSAT)*...
                       ( (isqdra - qdrtot*cdrain(sl))/DAQUIF - ...
                         DECSAT*cdrain(sl)*(POROS + P.sh.dap(inode)*KFSAT) );
            else
               cdrain(sl) = cdrain(sl) + dtsolu/(POROS + P.sh.dap(inode)*KFSAT)*...
                       ( isqdra/DAQUIF - DECSAT*cdrain(sl)* ...
                         (POROS + P.sh.dap(inode)*KFSAT) );
            end
            cseep(sl) = cdrain(sl);
         end
% --- flux to surface water from aquifer :: TO BE IMPLEMENTED YET
         if (SWBR == 1)
            sqsur(sl) = sqsur(sl) + qdrtot*cdrain(sl)*dtsolu;
         end

% --- flux through bottom of soil profile
         if (P.qbot > 0.0d0)
            sqbot(sl) = sqbot(sl) + P.qbot*cseep(sl)*dtsolu;
         else 
            sqbot(sl) = sqbot(sl) + P.qbot*P.cml(P.nz,sl)*dtsolu;
         end
         
% --- next 'solute'
    end% solute
    
% --- continue with next solute time step 'dtsolute'
end% while

% ___________________________________________________________________________________________ 
%  S O L U T E   B A L A N C E   C O M P O N E N T S 
% ___________________________________________________________________________________________              
for sl = 1:2
% --- current solute flux at bottom of soil column
    if (P.q(P.nz+1) > 0.0d0)
        isqbot(sl) = P.q(P.nz+1) * cseep(sl);
    else
        isqbot(sl) = P.q(P.nz+1) * P.cml(P.nz,sl);
    end
    
% === calculate solute balance components ========================

% --- total amount in soil profile
    P.sampro(sl) = csurf(sl) + sum( P.cmsy(:,sl) .* P.nodes.dz(1:P.nz) );% [g cm-2]

% --- add time step fluxes to total cumulative values
    sqprec(sl) = sqprec(sl) + P.nraidt * cpre(sl) * P.dt;
    sqirrig(sl) = sqirrig(sl) + P.nird * P.cirr(sl,P.tidx) * P.dt;

% --- cumulative solute balance
    % solbal : Cumulative solute balance for current balance period [g cm-2]
    %   -unsure about the sign of creactot, because in the original it is
    %    total decay while here it is total substance formation.
    %   -unable to understand how the signs of all term must be composed.
    solbal(sl) = P.sampro(sl) -sqprec(sl) -sqirrig(sl) -sqbot(sl) +creactot(sl) -P.samini(sl) +rottot(sl);% +sqdra
end

% -------------------------------------------------------------------------
% (CDE + SST) :: END
% -------------------------------------------------------------------------
%% end
return