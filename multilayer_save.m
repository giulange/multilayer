function multilayer_save(P)

if W.MTCL==1
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
end

% salva matrice di valori
%save pote_tot.dat P.A -ascii
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