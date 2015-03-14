function nitrati_IO_gui_soilconsweb( Itime, Otime )
% This function is called by the php engine installed in
% /var/www/cross-server/remote.php

% Read:
% {'ADD1'}						soil
% {'2015-01-01','2015-02-02'}	time
% {'Olivo'}						crop
% {'True'}						irrigazione di soccorso ( true/false/NaN )
% {'Domanda'}					tipo di irrigazione ( domanda/turno fisso/none ) *
% [5]							e la sua quantita' *
% {'Organica-Letame'}			tipo concimazione ( vari valori ) *
% [10]							e la sua quantita' *
fid = fopen(Itime, 'r');
for ii = 1:8
    eval( ['I{ii,1} = ',fgetl( fid ),';'] )
end
fclose(fid);

% Launch:


% Write:
fid = fopen(Otime, 'w');
fprintf(fid,'');
fclose(fid);

return