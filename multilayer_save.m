function multilayer_save( O, proj, Pj )
% NOTES:
% 
% Decide with Antonio: do we print time in rows and nodes in columns? At
% now I decided like this!
% 
% Implement: Set the check for Monte Carlo Simulation in multilayer_prog in
% order to split for multilayer_save_mcs.m or multilayer_save_norm.m.
% 
% Le matrici importanti sono le seguenti:
%   > O.C2          --> solutes concentrations      [nodes x times x montecarlo x solutes]
%   > O.h22         --> flussi ai nodi intermedi    [nodes x times x montecarlo]
%   > O.fluxsurf    --> flusso al contorno superiore[1 x times x montecarlo]
%   > O.fluxbot     --> flusso al contorno inferiore[1 x times x montecarlo]
%   > O.runoff      --> runoff                      [1 x times x montecarlo]
% 
%   > O.runon       --> in una implementazione futura.
%
% I farei in questo modo: conserverei solo le matrici utili
% elencate sopra, tal quali e senza rimaneggiamenti. inoltre in
% Initialization prevederei una preallocazione considerando
% anche il numero di simulazioni montecarlo cos da mantenere
% tutto quanto - ridotto ai minimi termini - in memoria RAM, e
% solo alla fine di tutto il run si salva quello che il
% programma DEVE salvare (le info da salvare possono anche
% essere configurate dall'utente, per cui si prevede una
% sezione in 'conf' apposita).
% 
% Remember also to account for M.nnc : this can be avoided when a
% simulation will never be abondoned!

%% flux at intermediate nodes
%  to test saving function you can do this:
%  O.h22 = rand(size(O.h22));

% CREATE filename
FILE = [proj.opath, O.files.h22];
% OVERWRITE existing ??
if exist( FILE, 'file' ) == 2
    warning( 'The %s file already exists! It will be re-written!!',FILE )
end
% OPEN file
fid = fopen( FILE, 'w' );
% WRITE file
for r = 1:Pj%size(O.h22,2)      % r is the row of file (not of array)
    for c = 1:size(O.h22,1)  % c is the col of file (not of array)
        fprintf(fid,' %7.3f',O.h22(c,r));
    end
    fprintf(fid,'\n');
end
% CLOSE file
fclose( fid );
%% concentration of solutes
%  to test saving function you can do this:
%  O.C2 = rand(size(O.C2));
%% --- ammonia NH4
% CREATE filename
FILE = [proj.opath,'NH4_',O.files.C2];
% OVERWRITE existing ??
if exist( FILE, 'file' ) == 2
    warning( 'The %s file already exists! It will be re-written!!',FILE )
end
% OPEN file
fid = fopen( FILE, 'w' );
% WRITE file
% TMP = squeeze( O.C2(:,:,1,1) );
for r = 1:Pj%size(O.C2,2)      % r is the row of file (not of array)
    for c = 1:size(O.C2,1)  % c is the col of file (not of array)
        fprintf(fid,' %7.3e',O.C2(c,r,1,1));
    end
    fprintf(fid,'\n');
end
% CLOSE file
fclose( fid );
%% --- nitrate NO3
% CREATE filename
FILE = [proj.opath, 'NO3_', O.files.C2];
% OVERWRITE existing ??
if exist( FILE, 'file' ) == 2
    warning( 'The %s file already exists! It will be re-written!!',FILE )
end
% OPEN file
fid = fopen( FILE, 'w' );
% WRITE file
for r = 1:Pj%size(O.C2,2)      % r is the row of file (not of array)
    for c = 1:size(O.C2,1)  % c is the col of file (not of array)
        fprintf(fid,' %7.3e',O.C2(c,r,1,2));
    end
    fprintf(fid,'\n');
end
% CLOSE file
fclose( fid );
%% flux
%  to test saving function you can do this:
%  O.fluxsurf   = rand(size(O.fluxsurf));
%  O.fluxbot    = rand(size(O.fluxbot));
%  O.runoff     = rand(size(O.runoff));
flux    = [ O.fluxsurf; O.fluxbot; O.runoff; ];
% CREATE filename
FILE = [proj.opath,O.files.flux];
% OVERWRITE existing ??
if exist( FILE, 'file' ) == 2
    warning( 'The %s file already exists! It will be re-written!!',FILE )
end
% OPEN file
fid = fopen( FILE, 'w' );
% WRITE file
for r = 1:Pj%size(flux,2)      % r is the row of file (not of array)
    for c = 1:size(flux,1)  % c is the col of file (not of array)
        fprintf(fid,' %7.3f',flux(c,r));
    end
    fprintf(fid,'\n');
end
% CLOSE file
fclose( fid );
%% end
return