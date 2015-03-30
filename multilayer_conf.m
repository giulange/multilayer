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
W.MTCL              = 1;

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
% nz:                   Number of nodes defining the geometry of the soil
%                       column during simulation of water flow. It includes
%                       the top and the botom boundaries.
W.nz                = 50;
% nlay:                 Number of soil layers.
W.nlay              = 3;
% zint:                 Soil layer bottom boundaries.
W.zint              = [25, 60, 300];

% SOIL GRID NODES WITH HYDRAULIC CHARACTERISTICS:
% W.crc.? --> curva ritenzione/conducibilit�: es. W.crc.dap, W.crc.tetas, ecc. 
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
W.hin               = repmat(-100,W.nz,1);
% W.inhin     = 1;
% W.hin       = load( fullfile(proj.ipath,'initial_inp.txt') )

W.tmax              = 120;
W.ntp               = 131;
% **DATA**
W.tp                = [0.01, 1:130];

%%   INITIAL INPUT [DELETE]
% ----------------------------------

I.description       = 'MARWA CORRECTED';    % DELETE
% **DATA**
I.hin               = repmat(-1000,W.nz,1); % --> as W.Ihin

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
% **DATA**
% times of flux/potential at bottom boundary.
B.bot.thqstar       = [ 0.00, 1.00, 1.50, 3.00, 4.00, 4.50, 5.00, 6.00, 7.00, 7.50 ];
% flux/potential at bottom boundary.
B.bot.hqstar        = [ 0.01, 0.05, 0.00, 0.01, 0.02, 0.00, 0.01, 0.02, 0.01, 0.00 ];

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
% **DATA**
%   [days]
B.Ctop.tqstar       = 0:1:B.Ctop.nCtop-1;
%   [�C]
B.Ctop.Tstar        = [ 22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	22.6	23.2	22.2	23.5	23.4	23.4	23.5	22.6	22.8	22	23.2	22.5	22.5	23.2	24.1	23.6	22.1	22.1	21.5	21.5	21.5	20.3	21.4	21.4	21.9	22	22	22.5	23.2	23.2	23.1	22.7	21.4	21.4	21.1	20.9	21.4	21	21.2	21.1	21.1	21	20.4	20.4	20.7	20.5	20.5	20.5	20.2	20.2	20.6	20	19.7	19	19	18.9	19.5	19.5	21.1	18.6	18.6	18.6	18.1	18.1	18.3	18.3	17.3	17.3	17.2	17.2	16.5	16.5	16	14.6	14.6	14.1	14.1	15.5	15.5	15.6	14.9	14.9	14.4	14.4	12.9	12.9	13.3	13.8	13.8	13.5	13.5	13.4	13.4	13.4	12.3	12.3	12.5	12.5	12.8	12.8	12.8	12.8	12.5	12.5	12.5	14	14	13.7	13.6 ];
%   [g/cm3H2O] -- FR: fertirrigazione on Top
B.Ctop.Cstar.NH.FR  = [ 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ];
%   [g/cm3H2O]
B.Ctop.Cstar.NO.FR  = [ 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ];
%   [g/cm2] -- SD: solid on dL
B.Ctop.Cstar.NH.SD  = [ 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ];
%   [g/cm2]
B.Ctop.Cstar.NO.SD  = [ 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ];
%   [g/cm2]
B.Ctop.Cstar.UR     = [ 0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0.00001	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ];
%   [g/cm2] -- rp: rapid
B.Ctop.Cstar.ORG.rp = [ 0.00432	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0.00001	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ];
%   [g/cm2] -- sw: slow
B.Ctop.Cstar.ORG.sw = [ 0.02565	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00001	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ];

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
% **DATA** [we can also use "load" from external file]
V.tqstar            = 0:1:V.nET-1;
V.ETr               = [ 0.961	0.361	0.361	1.325	0.314	0.489	0.489	0.472	0.489	0.489	0.911	0.489	0.489	0.492	0.389	0.417	0.417	0.407	0.407	0.647	0.647	0.833	0.833	0.8	0.69	0.69	0.685	0.685	0.518	0.518	0.5	0.858	0.858	0.661	0.661	0.638	0.638	0.6	0.578	0.578	0.623	0.623	0.516	0.516	0.6	0.52	0.52	0.451	0.451	0.443	0.443	0.4	0.4	0.4	0.526	0.526	0.463	0.463	0.4	0.526	0.526	0.611	0.611	0.6	0.6	0.6	0.402	0.402	0.598	0.598	0.593	0.593	0.6	0.637	0.637	0	0.527	0.527	0.5	0.5	0.575	0.575	0.55	0.499	0.499	0.5	0.5	0.5	0.5	0.499	0.499	0.5	0.5	0.508	0.508	0.5	0.462	0.462	0.5	0.5	0.518	0.518	0.5	0.426	0.426	0.5	0.5	0.463	0.463	0.5	0.437	0.437	0.4	0.4	0.401	0.401	0.4	0.4	0.4	0.4	0.4	0.568	0.568	0.4	0.405	0.405 ];
V.Kc                = [ 1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6	0.6 ];
V.LAI               = [ 0.13	0.19	0.23	0.26	0.38	0.42	0.47	0.598	0.651	0.712	0.78	0.856	0.938	1.026	1.12	1.22	1.324	1.432	1.545	1.661	1.78	1.902	2.026	2.152	2.28	2.409	2.539	2.669	2.8	2.93	3.061	3.19	3.319	3.446	3.572	3.696	3.818	3.938	4.055	4.169	4.281	4.389	4.494	4.596	4.694	4.788	4.878	4.963	5.045	5.122	5.195	5.262	5.326	5.384	5.437	5.486	5.529	5.568	5.601	5.629	5.652	5.67	5.682	5.69	5.692	5.689	5.681	5.668	5.65	5.627	5.599	5.567	5.529	5.488	5.441	5.391	5.336	5.277	5.215	5.148	5.078	5.005	4.929	4.849	4.767	4.682	4.595	4.505	4.414	4.322	4.228	4.133	4.037	3.941	3.844	3.748	3.652	3.557	3.463	3.371	3.281	3.193	3.107	3.025	2.946	2.871	2.8	2.734	2.673	2.618	2.569	2.526	2.49	2.465	2.443	2.422	2.403	2.386	2.37	2.354	2.34	2.327	2.315	2.303	2.293	2.283 ];
V.Droot             = [ 0.1	0.3	0.5	0.7	0.9	1.1	1.3	1.5	1.7	1.9	2.1	2.3	2.5	2.7	2.9	3.1	3.3	3.5	3.7	3.9	4.1	4.3	4.5	4.7	4.9	5.1	5.3	5.5	5.7	5.9	6.1	6.3	6.5	6.7	6.9	7.1	7.3	7.5	7.7	7.9	8.1	8.3	8.5	8.7	8.9	9.1	9.3	9.5	9.7	9.9	10.1	10.3	10.5	10.7	10.9	11.1	11.3	11.5	11.7	11.9	12.1	12.3	12.5	12.7	12.9	13.1	13.3	13.5	13.7	13.9	14.1	14.3	14.5	14.7	14.9	15.1	15.3	15.5	15.7	15.9	16.1	16.3	16.5	16.7	16.9	17.1	17.3	17.5	17.7	17.9	18.1	18.3	18.5	18.7	18.9	19.1	19.3	19.5	19.7	19.9	20.1	20.3	20.5	20.7	20.9	21.1	21.3	21.5	21.7	21.9	22.1	22.3	22.5	22.7	22.9	23.1	23.3	23.5	23.7	23.9	24.1	24.3	24.5	24.7	24.9	25.1 ];

%%   SOLUTE Jury INPUT
% ----------------------------------
if W.isol==1
S.J.description     = 'prova puglianello solute transport';
S.J.decad           = 0;
S.J.retard          = 0;
% **DATA**
S.J.tetasst         = repmat(0.500,W.nz,1);
S.J.sigma           = repmat(0.200,W.nz,1);
S.J.mu              = [ 0	0	0	0	0	0	0	0	0.114	0.226	0.326	0.417	0.5	0.577	0.648	0.715	0.778	0.836	0.892	0.945	0.995	1.042	1.088	1.131	1.173	1.213	1.251	1.288	1.324	1.359	1.392	1.424	1.456	1.486	1.515	1.544	1.572	1.599	1.625	1.651	1.676	1.7	1.724	1.747	1.77	1.792	1.814	1.835	1.856	1.876	1.896	1.916	1.935	1.954	1.972	1.991	2.009	2.026	2.043	2.06	2.077	2.093	2.109	2.125	2.141	2.156	2.171	2.186	2.201	2.216	2.23	2.244	2.258	2.272	2.285	2.298	2.312	2.325	2.337	2.35	2.363	2.375	2.387	2.399	2.411	2.423	2.434	2.446	2.457	2.469	2.48	2.491	2.501	2.512	2.523	2.533	2.544	2.554	2.564	2.574 ]';
S.J.Rcoeff          = repmat(1.000,W.nz,1);
S.J.Dcoeff          = repmat(0.000,W.nz,1);
end

%%   SOLUTE CDE INPUT
% ----------------------------------
if W.isol==2
S.CDE.description   = 'MARWA CORRECTED';
S.CDE.tCinput       = 0.000;
S.CDE.tCinput_end   = 0.500;
S.CDE.Cinput        = 0.040;
% Topt:                 The optimum temperature (�C) for the XXX process
S.CDE.Topt          = 25;
% NX:                   Where X = {'H':NH4, 'O':NO3}
S.CDE.NX.Kf1        = [0.100, 1.000]; % ['NH','NO']
S.CDE.NX.Kf2        = [1.000, 1.000];
S.CDE.NX.Kr         = [0.000, 1.000];
% **DATA**
S.CDE.Cin.NH        = repmat(0.0000,W.nz,1);
S.CDE.Cin.NO        = [ repmat(0.0002,11,1); repmat(0.0000,W.nz-11,1) ]; % you can parameterize "11" if it is useful!
S.CDE.lambda        = repmat(3.0,W.nz,1);

% Knitr:                Coeff. nitrificazione
S.CDE.Knitr         = [ repmat(1.0000,11,1); repmat(0.0100,W.nz-11,1) ];

% Kimmb:                Coeff. immobilizzazione
S.CDE.Kimmob        = [ repmat(0.0400,11,1); repmat(0.0300,W.nz-11,1) ];

% Kdntr:                Coeff. denitrificazione
S.CDE.Kdenitr       = [ repmat(0.0400,11,1); repmat(0.0200,W.nz-11,1) ];
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
