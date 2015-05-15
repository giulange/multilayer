%% ****DATA****
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

%% **B.top**
% times of flux/potential at top boundary.
% B.top.thqstar       = 0:1:B.top.ntbc-1;
% To update this with the new definition approach delineated above, I
% firstly have to solve the P.kk variable: that is I really need it how it
% is now? I think that if we define all time-dependent variables on a daily
% basis the program needs to update the values of this variables once every
% time a new day starts.

% flux/potential at top boundary.
% B.top.hqstar        = [ -0.961	-0.361	0	-1.325	-0.314	-0.489	-0.489	-0.472	-0.489	0	-0.911	-0.489	-0.489	-0.492	-0.389	-0.417	0	-0.407	0	-0.647	0	-0.833	0	-0.6	0	-0.69	0	-0.685	0	-0.518	0	0	-0.68	0	-0.661	0	-0.638	0	0	-0.578	0	-0.623	0	-0.516	-0.6	0	-0.52	0	-0.451	0	-0.643	0	-0.487	0	0	-0.526	0	-0.463	-0.5	0	-0.726	0	-0.611	0	-0.549	-0.7	0	-0.402	0	-0.598	0	-0.793	0	0	-0.537	0	0	-0.527	0	-0.5	0	-0.575	0	0	-0.499	0	-0.512	0	0	-0.499	0	-0.48	0	-0.508	0	0	-0.562	0	0	0	-0.518	0	0	-0.526	0	0	0	-0.463	0	0	-0.537	0	-0.5	0	-0.401	0	0	-0.544	0	0	0	-0.568	0	0	-0.405	0 ];

%% **B.bot**
%% **B.Ctop**
%% **V**
V.tqstar            = 0:1:V.nET-1;
% We might provide a meteo file:
V.ETr               = [ 0.961	0.361	0.361	1.325	0.314	0.489	0.489	0.472	0.489	0.489	0.911	0.489	0.489	0.492	0.389	0.417	0.417	0.407	0.407	0.647	0.647	0.833	0.833	0.8	0.69	0.69	0.685	0.685	0.518	0.518	0.5	0.858	0.858	0.661	0.661	0.638	0.638	0.6	0.578	0.578	0.623	0.623	0.516	0.516	0.6	0.52	0.52	0.451	0.451	0.443	0.443	0.4	0.4	0.4	0.526	0.526	0.463	0.463	0.4	0.526	0.526	0.611	0.611	0.6	0.6	0.6	0.402	0.402	0.598	0.598	0.593	0.593	0.6	0.637	0.637	0	0.527	0.527	0.5	0.5	0.575	0.575	0.55	0.499	0.499	0.5	0.5	0.5	0.5	0.499	0.499	0.5	0.5	0.508	0.508	0.5	0.462	0.462	0.5	0.5	0.518	0.518	0.5	0.426	0.426	0.5	0.5	0.463	0.463	0.5	0.437	0.437	0.4	0.4	0.401	0.401	0.4	0.4	0.4	0.4	0.4	0.568	0.568	0.4	0.405	0.405 ];
V.ETr               = V.ETr(1:W.tmax);
% V.ETr               = B.top.eto;

% We should develop a full crop module from which a LAI is simulated:
V.LAI               = [ 0.13	0.19	0.23	0.26	0.38	0.42	0.47	0.598	0.651	0.712	0.78	0.856	0.938	1.026	1.12	1.22	1.324	1.432	1.545	1.661	1.78	1.902	2.026	2.152	2.28	2.409	2.539	2.669	2.8	2.93	3.061	3.19	3.319	3.446	3.572	3.696	3.818	3.938	4.055	4.169	4.281	4.389	4.494	4.596	4.694	4.788	4.878	4.963	5.045	5.122	5.195	5.262	5.326	5.384	5.437	5.486	5.529	5.568	5.601	5.629	5.652	5.67	5.682	5.69	5.692	5.689	5.681	5.668	5.65	5.627	5.599	5.567	5.529	5.488	5.441	5.391	5.336	5.277	5.215	5.148	5.078	5.005	4.929	4.849	4.767	4.682	4.595	4.505	4.414	4.322	4.228	4.133	4.037	3.941	3.844	3.748	3.652	3.557	3.463	3.371	3.281	3.193	3.107	3.025	2.946	2.871	2.8	2.734	2.673	2.618	2.569	2.526	2.49	2.465	2.443	2.422	2.403	2.386	2.37	2.354	2.34	2.327	2.315	2.303	2.293	2.283 ];
V.LAI = V.LAI(1:W.tmax);

% We should develop a full crop module from which roots are simulated:
V.Droot             = [ 0.1	0.3	0.5	0.7	0.9	1.1	1.3	1.5	1.7	1.9	2.1	2.3	2.5	2.7	2.9	3.1	3.3	3.5	3.7	3.9	4.1	4.3	4.5	4.7	4.9	5.1	5.3	5.5	5.7	5.9	6.1	6.3	6.5	6.7	6.9	7.1	7.3	7.5	7.7	7.9	8.1	8.3	8.5	8.7	8.9	9.1	9.3	9.5	9.7	9.9	10.1	10.3	10.5	10.7	10.9	11.1	11.3	11.5	11.7	11.9	12.1	12.3	12.5	12.7	12.9	13.1	13.3	13.5	13.7	13.9	14.1	14.3	14.5	14.7	14.9	15.1	15.3	15.5	15.7	15.9	16.1	16.3	16.5	16.7	16.9	17.1	17.3	17.5	17.7	17.9	18.1	18.3	18.5	18.7	18.9	19.1	19.3	19.5	19.7	19.9	20.1	20.3	20.5	20.7	20.9	21.1	21.3	21.5	21.7	21.9	22.1	22.3	22.5	22.7	22.9	23.1	23.3	23.5	23.7	23.9	24.1	24.3	24.5	24.7	24.9	25.1 ];
V.Droot = V.Droot(1:W.tmax);
% V.Droot(90:end) = 0;
%% **S.J**
% How did you measured/calculated this parameter?
% We should define it using the same approache used for others
% depth-dependent variables.
S.J.mu              = [ 0	0	0	0	0	0	0	0	0.114	0.226	0.326	0.417	0.5	0.577	0.648	0.715	0.778	0.836	0.892	0.945	0.995	1.042	1.088	1.131	1.173	1.213	1.251	1.288	1.324	1.359	1.392	1.424	1.456	1.486	1.515	1.544	1.572	1.599	1.625	1.651	1.676	1.7	1.724	1.747	1.77	1.792	1.814	1.835	1.856	1.876	1.896	1.916	1.935	1.954	1.972	1.991	2.009	2.026	2.043	2.06	2.077	2.093	2.109	2.125	2.141	2.156	2.171	2.186	2.201	2.216	2.23	2.244	2.258	2.272	2.285	2.298	2.312	2.325	2.337	2.35	2.363	2.375	2.387	2.399	2.411	2.423	2.434	2.446	2.457	2.469	2.48	2.491	2.501	2.512	2.523	2.533	2.544	2.554	2.564	2.574 ]';