% -1-
if nr_breaked%  convergence has been reached!
% --- RISE   TIMESTEP
    if (p<=3)        ,  P.dt = min(P.dt*W.dtfactor,W.dtmax); end
% --- REDUCE TIMESTEP
    if (p>=W.maxit-2),  P.dt = max(P.dt*W.dtfactor,W.dtmin); end    
end

% -2-
if P.flEndOfDay% at previous dt I closed a DAY!
% --- I have to deal with the START-OF-DAY
    P.dt                = W.dtin;
    P.flEndOfDay        = true;

else
% --- I have to deal with the END-OF-DAY
    P.tidx              = floor(P.time(P.j-1)+P.dt)+1;
    P.L                 = P.tidx > P.tidx_jm1;
    if P.L% if END-OF-DAY
    % --- set P.dt to simulate the last piece of current day
        P.flEndOfDay    = true;
        P.time(P.j)     = P.tidx;
        P.dt            = P.tidx-P.time(P.j-1);
        % otherwise last dt of current DAY would use tidx of next DAY!
        P.tidx          = P.tidx - 1;
    end
end

% -3-
P.j                     = P.j+1;
P.time(P.j)             = P.time(P.j-1)+P.dt;
P.tidx_jm1              = P.tidx;
            
%% Code to be revised for inclusion

% -4-
% --- check maximum number of time steps during this day
% isteps = isteps + 1
% if (isteps .gt. msteps) then
% write(tmp,'(i11)') daynr
% tmp = adjustl(tmp)
% messag ='The maximum number of time steps for a day is exceeded'
% &    //' at daynumber '//trim(tmp)//'. Check input for numerical'
% &    //' solution of Richards equation'
% call fatalerr ('timer',messag)
% endif







% --- determine next time step, based on numerical performance
%       if (flDayStart) then
%         if(flprevious.eq.2)then
%            dt = max(dt, sqrt(dtmin*dtmax),dtprevious)
%         else
%            dt = max(dt, sqrt(dtmin*dtmax))
%         end if
%         dtprevious = dt
%       else
%         if(flprevious .eq. 2)then
%            dt = dtprevious
%         else
%            if (numbit.le.3)     dt = min(dt*2.0d0,DtMax)
%            if (numbit.ge.MaxIt) dt = max(dt*0.5d0,DtMin)
%            dtprevious = dt
%          endif
%       endif
%       flprevious = 1
% 
%       if ( tcum + dt - tEvent .gt. dtCrit) then
%          dt = tEvent - tcum
%          flprevious = 2
%          flTnext = .true.
%       endif
%       dt = max(dt,dtmin)
%       dt = min(dt,dtmax)

% add a piece found in headcalc.for at lines 667-713 !!!
    
    % --- save state variables of time = j
    %hm1         = h;
    %thetm1      = theta;
    %gwlm1       = gwl;
    %pondm1      = pond;
