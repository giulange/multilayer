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
%% **S.J**
% How did you measured/calculated this parameter?
% We should define it using the same approache used for others
% depth-dependent variables.
% S.J.mu              = [ 0	0	0	0	0	0	0	0	0.114	0.226	0.326	0.417	0.5	0.577	0.648	0.715	0.778	0.836	0.892	0.945	0.995	1.042	1.088	1.131	1.173	1.213	1.251	1.288	1.324	1.359	1.392	1.424	1.456	1.486	1.515	1.544	1.572	1.599	1.625	1.651	1.676	1.7	1.724	1.747	1.77	1.792	1.814	1.835	1.856	1.876	1.896	1.916	1.935	1.954	1.972	1.991	2.009	2.026	2.043	2.06	2.077	2.093	2.109	2.125	2.141	2.156	2.171	2.186	2.201	2.216	2.23	2.244	2.258	2.272	2.285	2.298	2.312	2.325	2.337	2.35	2.363	2.375	2.387	2.399	2.411	2.423	2.434	2.446	2.457	2.469	2.48	2.491	2.501	2.512	2.523	2.533	2.544	2.554	2.564	2.574 ]';