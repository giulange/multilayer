%% on general structure of the model
% Potrebbe essere un'alternativa risolvere prima la parte idrologica poi
% quella dei soluti!! Pensiamoci (non saprei in termini di performance...)
%% on geometry & time of multilayer
% Dear Antonio, we should think of a new kind of implementation about time
% & geometry of the simulated system:
%   > GEOMETRY
%     We have to implement different kinds of geometry that can be
%     activated according to different conditions (also set by the program
%     user). For instance if we define 100 nodes, they must comprise all
%     the nodes of the system (also for the two boundary conditions)
%     without adding misterious nodes in the meanwhile.
%     hsurf --> Qtop, 
% 
%   > TIME
%     We have to simulate each time step, but we should think to a
%     different way in which variables are stored in runtime. This means
%     that I only have to store all the j-s for current print-timestep
%     (e.g. [day]), while for each past timestep (e.g. days) I have only to
%     store the values of parameters for THE day (maybe it's useless to
%     hold the information for all the j-s).
%     The user must define the start/end dates of simulation, that we will
%     use to print the final results in output files!!
% 
%   > GENERAL COMMENT
%     The same logic used to simplify the time domain (regarding what
%     information to store in runtime to be smarter) can be applied to the
%     GEOMETRY of the system, at least for those variables that must be
%     printed: do we really need a print of nitrate concentration for each
%     node of the system (and for each time j)?
%     Maybe we can aggregate the information (and only if strictly required
%     by the problem at hand, hold all the information, but this can be set
%     by the user as input).

%% pedici usati per iterare sulla griglia di suolo -- check with Antonio
% Esistono 2 modalità di accesso alla griglia di suolo nel programma:
%   -la prima usa tutti e solo in nodi interni 1:P.nz
%   -la seconda include nella griglia anche il primo nodo esterno sotto la
%    griglia. Questo nodo non è definibile in termini idrologici quindi
%    bisognerebbe eliminare il nodo P.nz+1 dalle variabili matlab definite
%    su P.nz e magari creare nuove variabili per definire lo stato del nodo
%    P.nz+1. Viceversa se bisogna definire tutto il programma su P.nz+1
%    allora diventa direttamente P.nz=P.nz+1 il che vuol dire includere
%    il nodo P.nz+1 nella griglia di suolo. DECIDIAMO rispettando la
%    componente fisica del programma!!!
%% _ADE_ function for nitrogen :: check with Antonio
% idL   --> I changed its definition, is it ok?
% O.C2  --> grows on "nodes" dimension from P.nz to P.nz+1. We have to
%           decide the size of it (and of all other variables like it).
%% W.hsurf / W.qsurf
%   I cannot change the two as hqsurf, because they both coexist in
%   multilayer_prog.m (lines 325-326):
%                     W.qsurf         = P.Ep;
%                     W.hsurf         = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
% Either we change that, or we leave all as it is now.
%% P.Cinput
% It is preallocated, then it is defined in _prog:
%   P.Cinput(1)=...  &  P.Cinput(2)=...
% but only used in _ADE_
% Try to adjust things!!!!!
%% P.h1star [remove?]
%% fluxin
% In _ADE_N you are reading at P.nz+1 of P.flux, which is impossible (it's
% defined on 1:P.nz !!).
%% stats on montecarlo: cumulative std during Monte Carlo Simulation ??
% We should compute on the mm repetitions the following statistics for the
% Ti output of the i simulation, where i=[1, ..., M.nvp]:
%   ++mean++    Tave = sum( Ti ) / M.nvp
%   ++bias++    NULL
%   ++ std++    Tstd = sqrt( sum( (Ti - Tave)^2 ) / (M.nvp-1) )
%   ++ mse++    NULL
%   other?
% For instance, for O.C2 we have:
%   Tave = sum( Ti ) / M.nvp
%       where Ti = squeeze( O.C2(:,:,:,mm) )
%   Tstd = sqrt( sum( (Ti - Tave)^2 ) / (M.nvp-1) )
%       where Tave = Ti / M.nvp at current i, and finally the equation is:
%   (Ti - Ti/M.nvp)
%     O.C2 = O.C2 + squeeze( O.C2(:,:,mm,:) );
%     
%     % 
%     squeeze( O.C2(:,:,mm,:) ) - squeeze( O.C2(:,:,mm,:) )/M.nvp

% An example follows:
Ti          = rand(10,1);

Tave        = sum( Ti ) / length( Ti );

Tave_i = 0;
for ii=1:length(Ti)
    Tave_i  = Tave_i + Ti(ii) / length( Ti );
end
fprintf('mean:\n reference=%.4f, manual=%.4f, iterative=%.4f\n',mean(Ti), Tave, Tave_i);

Tstd        = sqrt( sum( (Ti - Tave).^2 ) / (length(Ti)-1) );
% Tstd_i      = zeros(size(Ti));
% for ii=1:length(Ti)
%     Tstd_i  = Tstd_i - Ti ./ (length( Ti )-1);
% end
% iterative is not possible due to the power 2 in expression!!
fprintf('std:\n reference=%.4f, manual=%.4f, iterative=%.4f\n',std(Ti), Tstd, NaN);



