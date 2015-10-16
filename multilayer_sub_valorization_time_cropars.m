function PPAR = multilayer_sub_valorization_time_cropars( PAR, rot_edate, Wsdate, ts, PPAR )
% 
% DESCRIPTION
%   This function sets values in time for crop parameters according to the
%   rotation dates (rot_edate).
%   Values are set in PPAR using a step function where each step is given
%   by a row defined in PAR
% 
% INPUTS
%   PAR:            The crop parameter to be valorised. It is assumed that
%                   in the first column are dates and in the second column
%                   are values. Valorization takes place using a step
%                   function.
%   rot_edate:      End date for current crop in rotation. This is
%                   generally equal to:
%                       V.rotation(V.currcrop,2)
%   Wsdate:         Starting date of simulation.
%   ts:             Timestep as defined in W.timestep.
%   PPAR:           The PAR defined in P structure array, which is defined
%                   in time from W.sdate to W.edate (current parameter is
%                   generally defined in a time-interval inside this
%                   extremes).
%                   
% OUTPUT
%   PPAR:           The same PPAR is returned with added the definition of
%                   PAR in the time-interval in which the crop exists.

%% pre
len_units_time      = @(ed,sd) (datenum(ed)-datenum(sd)+1)/ts;% +1 to deal with both extremes of interval 
%% main
if size(PAR,1)==1
    dat             = len_units_time( PAR(1,1),     Wsdate ) : ...
                      len_units_time( rot_edate,    Wsdate );
    PPAR(dat)       = PAR{1,2};
else
    for ii = 1:size(PAR,1)-1
        dat         = len_units_time( PAR(ii,1),    Wsdate ) : ...
                      len_units_time( PAR(ii+1,1),  Wsdate );
        PPAR(dat)   = PAR{ii,2};
    end
    dat             = len_units_time( PAR(end,1),   Wsdate ) : ...
                      len_units_time( rot_edate,    Wsdate );
    PPAR(dat)       = PAR{end,2};
end
%% end
return