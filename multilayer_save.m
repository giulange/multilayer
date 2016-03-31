function multilayer_save( O, proj, Pj, Wisol )
% multilayer_save( O, proj, Pj, Wisol )
%
% DESCRIPTION
%   The function saves the most important results from current simulation:
%       (1) pressure head [always]
%       (2) ammonia [if Wisol~=0]
%       (3) nitrate [if Wisol~=0]
% 
% INPUTS
%   O       :   A structure array storing all the most important simulation
%               information about water and solutes.
%   proj    :   General information about the current project.
%   Pj      :   The latest simulation step before the run safely stopped.
%   Wisol   :   The user defined setting about the if and/or type of solute
%               transport.
%                   *Wisol == 0, no solute transport simulation.
%                   *Wisol ~= 0, solute transport by a pre-defined method.
% 
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

%% printing time-step
fprintf('Defining printing time-step...\n')

% This is the array that define which j's the user want to print:
iPrint          = NaN(size(O.tprint));

print_issue     = zeros(size(O.tprint));
for ii = 1:numel(O.tprint)
    Diff        = abs(O.time - O.tprint(ii));
    [mnD,idD]   = min(Diff);
    iPrint(ii)  = idD;
    if mnD > O.tptolle
        print_issue(ii) = 1;
    end
end
sumI = sum(print_issue);
if ~sumI
    fprintf('   All printing time-steps are correctly defined [tollerance is %5.5f]\n',O.tptolle)
else
    print_issue = find(print_issue);
    el_iss      = iPrint(print_issue~=0);% locations in "time" var
    fprintf('   Found %d issues.\n',sumI)
    fprintf('     %4s, %10s, %10s\n','idx','tprint','time')
    for ii = 1:sumI
        fprintf('     %4d, %10.4f, %10.4f\n',ii, O.tprint(print_issue(ii)), O.time(el_iss(ii)) )
    end
end
fprintf('...done!\n')
%% Soil Water Pressure Head (h) & Volumetric Soil Water Content (teta) [mandatory]
%  to test saving function you can do this:
%  O.h22 = rand(size(O.h22));
fprintf('Saving pressure head at user-defined time-step [W.tp]... ')
% CREATE filename
FILE  = [proj.opath, O.files.h22];
FILE2 = [proj.opath, O.files.theta];
% OVERWRITE existing ??
if exist( FILE, 'file' ) == 2
    warning( 'The %s file already exists! It will be re-written!!',FILE )
end
if exist( FILE2, 'file' ) == 2
    warning( 'The %s file already exists! It will be re-written!!',FILE2 )
end
% OPEN file
fid  = fopen( FILE,  'w' );
fid2 = fopen( FILE2, 'w' );
% WRITE file
%  -- col header is DEPTH [cm]
% fprintf(fid,' %10d',NaN);% I should put this: '[day] / [cm]'
fprintf(fid,' %10s','');% I should put this: '[day] / [cm]'
fprintf(fid2,' %10s','');% I should put this: '[day] / [cm]'
for c = 1:size(O.h22,1)
    fprintf(fid, ' %10.3f',O.z(c));
    fprintf(fid2,' %10.3f',O.z(c));
end
fprintf(fid, '\n');
fprintf(fid2,'\n');
for r = 1:numel(iPrint)      % r is the row of file (not of array)
    fprintf(fid, ' %10.3f',O.time(iPrint(r))); % row header is TIME [day]
    fprintf(fid2,' %10.3f',O.time(iPrint(r))); % row header is TIME [day]
    for c = 1:size(O.h22,1)  % c is the col of file (not of array)
        fprintf(fid, ' %10.3f',O.h22 (c,iPrint(r)));
        fprintf(fid2,' %10.3f',O.teta(c,iPrint(r)));
    end
    fprintf(fid, '\n');
    fprintf(fid2,'\n');
end
% CLOSE file
fclose( fid  );
fclose( fid2 );
fprintf(' done!\n')
%% concentration of solutes [optional]
%% --- ammonia NH4
if Wisol
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
else
    fprintf('Ammonia simulation not saved! [solute transport unselected ––> W.isol=%d]\n',Wisol)
end
%% --- nitrate NO3
if Wisol
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
else
    fprintf('Nitrate simulation not saved! [solute transport unselected ––> W.isol=%d]\n',Wisol)
end
%% flux at boundaries [mandatory]
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