function multilayer_save( O )
%   > O.C2          --> concentrazioni soluti [nodi x tempi x soluti x montecarlo]
%   > P.h22         --> flussi ai nodi intermedi [nodi x tempi x montecarlo]
%   > O.fluxsurf    --> flusso al contorno superiore; [tempi x montecarlo]
%   > P.fluxbot     --> flusso al contorno inferiore; [tempi x montecarlo]
%   > P.runoff      --> runoff; [tempi x montecarlo]
%% da scommentare in caso di simulazione MONTECARLO
[r,c] = size(C);

%OPEN file
fid = fopen( strcat(proj.path,'teta_print',num2str(mm),'.dat'), 'w' );
%cicli loop
for riga=1:r
    for colonna = 1:c
        v = C(riga,colonna);
        if  colonna == c
            fprintf(fid,'%f\n',v);
        else
            fprintf(fid,'%f ',v);
        end
    end
end
state = fclose(fid);
close;

M.tetasum=C+M.tetasum;
M.tetasumSQ=(C.^2)+M.tetasumSQ;


%%%%MONTECARLO CONCENTRAZIONI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rr,cc] = size(D);

%OPEN file
fid = fopen( strcat(proj.path,'conc_print',num2str(mm),'.dat'), 'ww' );
%cicli loop
for riga=1:rr
    for colonna = 1:cc
        vv = D(riga,colonna);
        if  colonna == cc
            fprintf(fid,'%f\n',vv);
        else
            fprintf(fid,'%f ',vv);
        end
    end
end
state = fclose(fid);
close;

M.concsum=D+M.concsum;
M.concsumSQ=(D.^2)+M.concsumSQ;

%%%%MONTECARLO FLUSSI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rrr,ccc] = size(G);

%OPEN file
fid = fopen([strcat(proj.path,'flux_i_j',num2str(mm),'.dat')],'ww');
%cicli loop
for riga=1:rr
    for colonna = 1:ccc
        vvv = G(riga,colonna);
        if  colonna == ccc
            fprintf(fid,'%f\n',vvv);
        else
            fprintf(fid,'%f ',vvv);
        end
    end
end
state = fclose(fid);
close;

M.fluxsum=G+M.fluxsum;
M.fluxsumSQ=(G.^2)+M.fluxsumSQ;


tetamed=M.tetasum/(mm-M.nnc);
tetavar=(mm*M.tetasumSQ-M.tetasum.^2)/((mm-M.nnc)*(mm-M.nnc-1));

concmed=M.concsum/(mm-M.nnc);
concvar=(mm*M.concsumSQ-M.concsum.^2)/((mm-M.nnc)*(mm-M.nnc-1));

fluxmed=M.fluxsum/(mm-M.nnc);
fluxvar=(mm*M.fluxsumSQ-M.fluxsum.^2)/((mm-M.nnc)*(mm-M.nnc-1));

eval(['save ''' proj.path 'tetamed.dat'' tetamed -ascii'])
eval(['save ''' proj.path 'tetavar.dat'' tetavar -ascii'])

eval(['save ''' proj.path 'concmed.dat'' concmed -ascii'])
eval(['save ''' proj.path 'concvar.dat'' concvar -ascii'])

eval(['save ''' proj.path 'fluxmed.dat'' fluxmed -ascii'])
eval(['save ''' proj.path 'fluxvar.dat'' fluxvar -ascii'])

% salva matrice di valori
eval(['save ''' proj.path 'pote_tot.dat'' P.A -ascii'])
eval(['save ''' proj.path 'pote_print.dat'' B -ascii'])
eval(['save ''' proj.path 'teta_print.dat'' C -ascii'])
eval(['save ''' proj.path 'conc_NH4print.dat'' DNH -ascii'])
eval(['save ''' proj.path 'conc_NOprint.dat'' DNO -ascii'])
eval(['save ''' proj.path 'sink_print.dat'' E -ascii'])
eval(['save ''' proj.path 'cap_print.dat'' F -ascii'])
eval(['save ''' proj.path 'flux_i_j.dat'' G -ascii'])
eval(['save ''' proj.path 'runoff.dat'' RN -ascii'])
eval(['save ''' proj.path 'flux_surf_bot.dat'' Qtb -ascii'])
