%% multilayer_wtus_newtonraphson
% This module solve the Richard's equation with finite differences by means
% of a Newton-Raphson iteration scheme. This approach is similar to that
% implemented in SWAP-32.
%% Issues to solve
%   (1) k for implicit calculation of conductivities
%       a.  It is missing in the calculation of P.K (and maybe on other
%           related variables).
%       b.  Check that all terms depending on k=0/1 are correctly defined
%           in order to properly run each of the two flags (0 or 1).
%   (2) P.nodes.z
%       a.  Check that our definition does agree with SWAP equations in
%           which they consider an opposite definition of z.
%       b.  I have to understand the exact position of nodes with respect
%           to compartment limits!!
%   (6) develop sinks, such as: root, drain and macroporosity.
%% **IMPLICIT LINEARIZATION OF HYDRAULIC CONDUCTIVITIES
if W.SwkImpl==0
% --- 0: internodal conductivities refer to the old time level j-1
    k	= 0;
elseif W.SwkImpl==1
% --- 1: internodal conductivities refer to the new time level j
    k	= 1;
end
%% Nodal/Internodal Hydraulic Conductivity [P.K,P.Kim2,P.Kip2] <-- h,teta
multilayer_Kis2
%% top      boundary condition <-- h,pond_jm1,K,
% --- you should decrease time step to dtmin if required by boundtop
%     see SWAP-32, TimeControl(3):629
multilayer_boundtop
%% root, drain, macropores [sink,macr] <-- h_jm1
% ** sink **
% %      |--root--|   |--drain-|
% sink = zeros(P.nz,1) + zeros(P.nz,1);
% **TEMPORARY SINK:
for ii=1:P.nz
    if W.iveg==1
        if P.Droot(P.tidx)>0
            if W.iosm==1 && V.ifs>3
                P.op        = -P.ECstar(ii)*360;
            else
                P.op        = 0;
            end
            if P.nodes.z(ii) < P.Droot(P.tidx)                 
                P.sink(ii)   = fnsink( P.h_jm1(ii), P.nodes.z(ii),P.op,P.ptra,P.Droot(P.tidx),V.ifg,V.ifs,V);
                %P.sink(ii)   = fnsink_OLD(P.h_jm1(ii),P.nodes.z(ii),P.op,P.ptra,V.hI,V.hII,V.hIIIH,V.hIIIL,V.hIV,V.hw50,V.pw1,V.hs50,V.ps1,V.aMH,V.bMH,P.Droot(P.tidx),V.ifs,V.rda,V.rdb,V.rdc,P.nz,V.ifg,V.zc,V.g0,V.gzc,V.Drf);
            else
                P.sink(ii)   = 0;
            end
        else
            P.sink(ii)       = 0;
        end
    end
end
%      |--MaPo--| i.e. macroporosity
if W.SwMacro
    multilayer_macropore% --> output{ macr, ... }
end
%% pondrunoff   <-- h,pond_jm1
multilayer_pondrunoff
%% F-function [Fi #1] <---{P.Kim2,P.Kip2, P.h, P.teta,P.teta_jm1}
multilayer_Fi
% Estimate Fi inner product at init (sumFi2_pm1):
sumFi2_pm1 = sum(Fi.^2)*0.5;
%% Maximum No of iterations, according to P.dt [pmax]
if P.dt == W.dtmin
    pmax    = W.maxit*2;
else
    pmax    = W.maxit;
end
%% **Newton-Raphson iteration scheme [loop on "p"]
for p=1:pmax
%% save values of h
    % hold:
    P.h_hold = P.h;% same as hold(i), headcalc.for:278
%% dF_dh
    % **UPDATE conditions at the current iteration-step (at "p")
    %   considering the actual lambda.
    %   -[dKi_dhi]  (y) --> Conductivity derivative to the pressure head
    %   -[Ci]       (y) --> Capillary Capacity
    %   -[teta_pm1] (n) --> Moisture fraction
    [dKi_dhi,Ci,~] = multilayer_dKi_dhi(P.h, P.sh, 1:P.nz, P.dt);
    multilayer_dFdh% requires: {dKi_dhi,P.K, Ci, P.h, P.Kim2,P.Kip2}
%% dh_p --- tridiagonal system solution
    % Tri-diagonal system of equations (SWAP-32, Eq. 2.29, page 36):
    dh_p    = (dF_dh \ Fi); % note that inv(A)*b = A\b
    % I should implement a more general solution of the Thomas algorithm
    % that is a simplified form of Gaussian elimination (used for
    % tri-diagonal systems of equations).

    % Add a check on the solution provided by MatLab (e.g. singularity), so
    % that a different method can be accomplished!
    
    % ...develop the alternative SOLVER to provide solution at current
    %    p-step.
    % dh_p = ...
%% exit point
    if sum(isnan(dh_p))
        error('%s\n\tP.time=%.3f, P.j=%d\n', ...
              'Unable to properly solve the tri-diagonal system!',...
              P.time(P.j), P.j)
    end
%% BackTrack :: lambda [lambda=lambda/3]
    lambda      = 1;
    % backtrack block was broken:
    bt_breaked  = false;
    for iL = 1:W.maxbacktrack
%% factmax
        if P.dt==W.dtmin && p>pmax
            warning('%s\n\tP.time=%.3f, P.j=%d, p=%d\n', ...
                'factmax used:',P.time(P.j), P.j, p)
            % use the factmax to update the Newton step:
            factmax = 0.0d0;
            for ii = 1:P.nz
                if(abs( P.h_hold(ii) ) < 1.0d0 )
                    factmax = max( factmax, abs( dh_p(ii) ) ) ;
                else
                    factmax = max( factmax, abs( dh_p(ii) / P.h_hold(ii) ) );
                end
            end
            for ii = 1:P.nz
                P.h(ii) = P.h_hold(ii) - dh_p(ii) * min(1.0d0, 1.0d0 / factmax);
            end
        else
%% Newton-Raphson step :: at current lambda
            % regular case
            P.h = P.h_hold - lambda*dh_p;
        end
%% theta [!! teta update !!]
        % Moisture at "p" iteration:
        %   -check that teta is updated outside headcalc in SWAP-32, so
        %    that I have to do the same!! <-- VERY IMPORTANT
        %P.teta = fnteta( P.h, P.sh, 1:P.nz );
        P.teta = multilayer_fnteta_vgm( P.h, P.sh, 1:P.nz );
%% implicit [k=1]
        if k==1
%% root extraction
            %multilayer_rootextraction
%% Nodal/Internodal Hydraulic Conductivity [P.K,P.Kim2,P.Kip2] <-- h,teta
            multilayer_Kis2
        end
%% top boundary condition
        if W.SwMacro, QMpLatSsSav = QMpLatSs; end
        multilayer_boundtop
%% macropore
        fl_unsat_ok(3) = 0;
        if W.SwMacro && ~fl_unsat_ok(3)
            multilayer_macropore
        elseif W.SwMacro && fl_unsat_ok(3)
            QMpLatSs = QMpLatSsSav;
        end
%% pondrunoff
        multilayer_pondrunoff
%% F-function [Fi #2]
        % The original headcalc.for/SWAP-32 implementation uses the same
        % piece of code for computing the F-function here (Fi #2) and before
        % the Newton-Raphson iteration scheme (which was called Fi #1).
        % The only differences (what Fi #2 has and not Fi #1) are at lines:
        %   -(193,195)
        %   +(445,450,452,453,454)
%         figure(6),subplot(211),plot(1:100,Fi),subplot(212),plot(1:100,diag(dF_dh))%,axis([-1,101,-1,1])
        multilayer_Fi% requires: {P.Kim2,P.Kip2, P.h, P.teta,P.teta_jm1}
        % Estimate Fi inner product at current "p" (sumFi_p):
        sumFi2_p    = sum(Fi.^2)*0.5;
        % *IMPORTANT
        %   -BackTrack block exit criterion:
        maxFi_p     = max(abs(Fi));
        %   -Newton-Raphson block exit criterion:
        sumFi_p     = sum(Fi); % [cm d-1]
%% break of BackTrack block
        if sumFi2_p<sumFi2_pm1 || maxFi_p<W.CritDevBalCp
            bt_breaked = true;
            break
        end
%% update iteration progress
        % If Newton step is too large, repeat using reduced lambda:
        lambda      = lambda/3;
    end
%% check
    % Do something if you exit backtrack block with bt_breaked=false (it
    % means that you were not able to apply a reduced Newton-step as
    % expected)!!
    if ~bt_breaked
        % do something
    end
%% break Newton-Raphson iteration?
    nr_breaked  = true;
    isAdj_flNonConv1 = false;
% --- APPLY PERFORMANCE CRITERIA PER COMPARTMENT
% --- Test for water balance deviation of soil compartments
    flNonConv1 = abs(Fi)>W.CritDevBalCp;
% --- ESCAMOTAGE*, while waiting for finding the bug on bottom node
    % if only the bottom node is non-convergent
    if ~sum(flNonConv1(1:end-1)) && flNonConv1(end)
        % then apply a different threshold to bottom node and continue
        if abs(Fi(end))<1.0d-01
            flNonConv1(end) = 0;
            isAdj_flNonConv1 = true;
        end
    end
    if sum(flNonConv1), nr_breaked=false; end
% --- Test for change of pressure head
    flNonConv2 = abs(P.h_hold)<1.0d0 & abs(P.h-P.h_hold)>W.CritDevh2Cp;
    if sum(flNonConv2), nr_breaked=false; end
    flNonConv3 = abs(P.h-P.h_hold)./abs(P.h_hold)>W.CritDevh1Cp;
    if sum(flNonConv3), nr_breaked=false; end
    
% --- Test for water balance of ponding layer
    if ~isflux
        P.qtop = -P.Kim2(1)*((hsurf - P.h(1))/P.nodes.disnod(1)+1.0d0);
        if nr_breaked && ~W.SwMacro
            deviat = pond - pond_jm1 + reva*P.dt - (P.nraidt+P.nird+melt)*P.dt ...
                    - runon*P.dt  +  runoff  - P.qtop * P.dt;
            if abs(deviat) > W.CritDevPondDt
                flnonconv3 = true;
                nr_breaked = false;
           end
        end
    end
    if W.SwMacro% always FALSE untill I develop it
        deviat = pond - pond_jm1 + reva*P.dt - (P.nraidt+P.nird+melt)*P.dt ...
                    - runon*P.dt  +  runoff  - P.qtop * P.dt; %...
                  % + ArMpSs*(P.nraidt+P.nird+melt)*P.dt + QMpLatSs;
        if abs(deviat) > W.CritDevPondDt
            flnonconv3 = true;
            nr_breaked = false;
        end            
        error('To be completed!')
    end

% --- Critical water balance deviation for the whole soil profile
    if abs(sumFi_p) > W.CritDevBalTot
        nr_breaked = false;
        if isAdj_flNonConv1% can be removed! --> ESCAMOTAGE*
            if sum(Fi(1:end-1)) < W.CritDevBalTot% can be removed! --> ESCAMOTAGE*
                nr_breaked = true;% can be removed! --> ESCAMOTAGE*
            end
        end
    end

% --- EXIT WITH SUCCESS !!
    if nr_breaked
        break
    else
%         plot_performance_criteria
    end
%% what to do with this??
    % %     % *****maybe this must be solved after solving h at "p+1"******
    % %     % First order approximation of the Taylor expansion to calculate the
    % %     % new moisture fraction (SWAP 32 manual, Eq. 2.27, page 35):
    % %     % --not sure about the epoch of "p" and how to calculte hp at this
    % %     %   stage!
    % %     teta_p  = teta_pm1 + (h_p-h_pm1) .* Ci;% + ... + ...
    % %     %************
%% update loop
    % sum of Fi used in the backtrack block at step "p+1":
    sumFi2_pm1  = sumFi2_p;% sumold of SWAP
end% p
%% non-convergent
% --- see TimeControl(3):612
if ~nr_breaked%     convergence could not be reached!
    my_dtfactor         = 2;
    n_noconv            = n_noconv +1;
    % *REDUCE TIMESTEP:
    % Pay attention to prevent very small dt values
    dtprevious          = P.dt;
    if P.dt > my_dtfactor*W.dtmin
        P.dt            = P.dt / my_dtfactor;
        P.time(P.j)     = P.time(P.j) -dtprevious +P.dt;
        % If we go back to the previous integer timestep, we must update:
        P.tidx          = floor(P.time(P.j))+1;
        % The dt before reduction––>(P.dt/my_dtfactor) could close a day
        % if the P.L=1! After reduction the current dt cannot close the day
        % and the program have to know this by setting the flEndOfDay flag
        % to false!
        if P.L
            P.flEndOfDay = false;
        end
        
    elseif P.dt > W.dtmin
        P.dt            = W.dtmin;
        P.time(P.j)     = P.time(P.j) -dtprevious +P.dt;
        % If we go back to the previous integer timestep, we must update:
        P.tidx          = floor(P.time(P.j))+1;
        % The dt before reduction––>(P.dt/my_dtfactor) could close a day
        % if the P.L=1! After reduction the current dt cannot close the day
        % and the program have to know this by setting the flEndOfDay flag
        % to false!
        if P.L
            P.flEndOfDay = false;
        end

    elseif P.dt==W.dtmin
        % *CONTINUE WITHOUT CONVERGENCE:
        % ...what to do?
        fl_noconv       = true;% continue without convergence?
    else
        error('Something wrong in setting P.dt!')
    end
    
    % RESTORE:
    if ~fl_noconv%  if dt can be reduced
% --- reset soil state variables of time = j-1
        P.h             = P.h_jm1;
        P.teta          = P.teta_jm1;
        %kmean(numnod+1) = k(numnod);
        %gwl         = gwlm1;
        pond            = pond_jm1;
    else%           if dt cannot be reduced
        % do nothing to hold current state as good one!
        warning('\n%s\n%s\n %7s\t%-11.6f\n %7s\t%-4d\n %7s\t%-g\n %7s\t%-2d\n%s',...
            repmat('_',1,35), 'non-convergence of Newton-Raphson:',...
            'time:',P.time(P.j),'j:',P.j,'dt:',P.dt,'p:',P.iter(1,P.j),...
            repmat('_',1,35) )
    end
end
%% clean
clear ii my_dtfactor