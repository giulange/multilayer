function multilayer_save_mcs( O, proj )
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
mean    = nanmean( O.h22, 3 );
std     =  nanstd( O.h22, 1, 3 );

% CREATE filename
FILE{1} = [proj.opath, 'mean_',O.files.h22];
FILE{2} = [proj.opath, 'std_',O.files.h22];
% OVERWRITE existing ??
for ii = 1:length(FILE)
    if exist( FILE{ii}, 'file' ) == 2
        warning( 'The %s file already exists! It will be re-written!!',FILE{ii} )
    end
end
% OPEN file
fid(1) = fopen( FILE{1}, 'w' );
fid(2) = fopen( FILE{2}, 'w' );
% WRITE file
for r = 1:size(mean,2)      % r is the row of file (not of array)
    for c = 1:size(mean,1)  % c is the col of file (not of array)
        fprintf(fid(1),' %7.3f',mean(c,r));
        fprintf(fid(2),' %7.3f',std (c,r));
    end
    fprintf(fid(1),'\n');
    fprintf(fid(2),'\n');
end
% CLOSE file
fclose( fid(1) );
fclose( fid(2) );
%% concentration of solutes
%  to test saving function you can do this:
%  O.C2 = rand(size(O.C2));
%% --- ammonia NH4
mean    = nanmean( O.C2(:,:,:,1), 3 );
std     =  nanstd( O.C2(:,:,:,1), 1, 3 );

% CREATE filename
FILE{1} = [proj.opath, 'mean_NH4_',O.files.C2];
FILE{2} = [proj.opath, 'std_NH4_',O.files.C2];
% OVERWRITE existing ??
for ii = 1:length(FILE)
    if exist( FILE{ii}, 'file' ) == 2
        warning( 'The %s file already exists! It will be re-written!!',FILE{ii} )
    end
end
% OPEN file
fid(1) = fopen( FILE{1}, 'w' );
fid(2) = fopen( FILE{2}, 'w' );
% WRITE file
for r = 1:size(mean,2)      % r is the row of file (not of array)
    for c = 1:size(mean,1)  % c is the col of file (not of array)
        fprintf(fid(1),' %7.3f',mean(c,r));
        fprintf(fid(2),' %7.3f',std (c,r));
    end
    fprintf(fid(1),'\n');
    fprintf(fid(2),'\n');
end
% CLOSE file
fclose( fid(1) );
fclose( fid(2) );
%% --- nitrate NO3
mean    = nanmean( O.C2(:,:,:,2), 3 );
std     =  nanstd( O.C2(:,:,:,2), 1, 3 );

% CREATE filename
FILE{1} = [proj.opath, 'mean_NO3_',O.files.C2];
FILE{2} = [proj.opath, 'std_NO3_',O.files.C2];
% OVERWRITE existing ??
for ii = 1:length(FILE)
    if exist( FILE{ii}, 'file' ) == 2
        warning( 'The %s file already exists! It will be re-written!!',FILE{ii} )
    end
end
% OPEN file
fid(1) = fopen( FILE{1}, 'w' );
fid(2) = fopen( FILE{2}, 'w' );
% WRITE file
for r = 1:size(mean,2)      % r is the row of file (not of array)
    for c = 1:size(mean,1)  % c is the col of file (not of array)
        fprintf(fid(1),' %7.3f',mean(c,r));
        fprintf(fid(2),' %7.3f',std (c,r));
    end
    fprintf(fid(1),'\n');
    fprintf(fid(2),'\n');
end
% CLOSE file
fclose( fid(1) );
fclose( fid(2) );
%% flux
%  to test saving function you can do this:
%  O.fluxsurf   = rand(size(O.fluxsurf));
%  O.fluxbot    = rand(size(O.fluxbot));
%  O.runoff     = rand(size(O.runoff));
mean    = [ nanmean( O.fluxsurf, 3 );
            nanmean( O.fluxbot, 3 );
            nanmean( O.runoff, 3 ); ];
std     = [ nanmean( O.fluxsurf, 3 );
            nanmean( O.fluxbot, 3 );
            nanmean( O.runoff, 3 ); ];
% CREATE filename
FILE{1} = [proj.opath, 'mean_',O.files.flux];
FILE{2} = [proj.opath, 'std_',O.files.flux];
% OVERWRITE existing ??
for ii = 1:length(FILE)
    if exist( FILE{ii}, 'file' ) == 2
        warning( 'The %s file already exists! It will be re-written!!',FILE{ii} )
    end
end
% OPEN file
fid(1) = fopen( FILE{1}, 'w' );
fid(2) = fopen( FILE{2}, 'w' );
% WRITE file
for r = 1:size(mean,2)      % r is the row of file (not of array)
    for c = 1:size(mean,1)  % c is the col of file (not of array)
        fprintf(fid(1),' %7.3f',mean(c,r));
        fprintf(fid(2),' %7.3f',std (c,r));
    end
    fprintf(fid(1),'\n');
    fprintf(fid(2),'\n');
end
% CLOSE file
fclose( fid(1) );
fclose( fid(2) );
%% end
return