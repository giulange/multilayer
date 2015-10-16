%% update: { j, dt, time, tidx, flEndOfDay }

% -0-
P.tidx_jm1              = P.tidx;

% -1-
if P.iter(2,P.j)%  convergence has been reached!
% --- RISE   TIMESTEP
    if (P.iter(1,P.j)<=3),      P.dt = min(P.dt*W.dtfactor,W.dtmax); end
% --- REDUCE TIMESTEP
    if (P.iter(1,P.j)>=W.maxit),P.dt = max(P.dt/W.dtfactor,W.dtmin); end
%     if (P.iter(1,P.j)>=10),     P.dt = max(P.dt/W.dtfactor,W.dtmin); end
end

% -2-
P.j                     = P.j+1;

% -3-
if P.flEndOfDay% at previous j I closed a DAY!
% --- I deal with the START-OF-DAY
    P.flEndOfDay        = false;
    P.flStartOfDay      = true;
    flZeroCumu          = true;% reset solute balance vars
    if B.top.rain(P.tidx_jm1+1)<=4.0d-01
        P.dt            = W.dtmax;
    elseif B.top.rain(P.tidx_jm1+1)<=4.0d-01 && P.iter(1,P.j-1)<=25
        P.dt            = min(P.dt*W.dtfactor*2,W.dtmax);
    elseif diff(B.top.rain(P.tidx_jm1:P.tidx_jm1+1)) <= 1.5d-00 && P.iter(1,P.j-1)<=20
        P.dt            = min(P.dt*W.dtfactor*1.5,W.dtmax);
    elseif diff(B.top.rain(P.tidx_jm1:P.tidx_jm1+1)) > 8.0d-00% 2.0d-01
%     if diff(B.top.rain(P.tidx_jm1:P.tidx_jm1+1)) > 8.0d-00% 2.0d-01
        P.dt            = W.dtin;%max(P.dt,sqrt(W.dtmin*W.dtmax));%sqrt(W.dtmin*W.dtmax);% sqrt(W.dtmin*W.dtmax) % W.dtin
    elseif P.iter(1,P.j-1)>3 && P.iter(1,P.j-1)<=10
        P.dt            = min(P.dt*W.dtfactor,W.dtmax);
    else
        P.dt            = max(P.dt,sqrt(W.dtmin*W.dtmax));
    end
    P.time(P.j)         = P.time(P.j-1)+P.dt;
    P.tidx              = floor(P.time(P.j))+1;

else
    P.flStartOfDay      = false;
    flZeroCumu          = false;% maintain solute balance vars
% --- is the now-defined dt to much near to the next tidx?
    if (P.tidx>(P.time(P.j-1)+P.dt)) &&  (P.tidx-(P.time(P.j-1)+P.dt) < 1.0d-02)
        P.dt            = P.dt + (P.tidx-(P.time(P.j-1)+P.dt));
    end
% --- set the current DAY
    P.time(P.j)         = P.time(P.j-1)+P.dt;
    P.tidx              = floor(P.time(P.j))+1;
% --- am I at the END-OF-DAY?
    P.L                 = P.tidx > P.tidx_jm1;
% --- I deal with the END-OF-DAY:
    if P.L% if END-OF-DAY
    % --- set P.dt to simulate the last piece of current day
        P.flEndOfDay    = true;
        P.time(P.j)     = P.tidx_jm1;
        P.dt            = P.time(P.j)-P.time(P.j-1);
        % otherwise last dt of current DAY would use tidx of next DAY!
        P.tidx          = P.tidx - 1;
    end
end          
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
