%% multilayer_wtus_newtonraphson
% This module solve the Richard's equation with finite differences by means
% of a Newton-Raphson iteration scheme. This approach is similar to that
% implemented in SWAP-32.
%% Issues to solve
%   (1) k for implicit calculation of conductivities
%       a.  It is missing in the calculation of K (and maybe on other
%           related variables).
%       b.  Check that all terms depending on k=0/1 are correctly defined
%           in order to properly run each of the two flags (0 or 1).
%   (2) P.nodes.z
%       a.  Check that our definition does agree with SWAP equations in
%           which they consider an opposite definition of z.
%       b.  I have to understand the exact position of nodes with respect
%           to compartment limits!!
%       
%   (3) define initial conditions!!
%   (4) define top boundary conditions in different settings!!
%   (5) define bottom boundary conditions in different settings!!
%   (6) develop sinks, such as: root, drain and macroporosity.
%% NOTES
%   ??
%% initialize
% partial derivative of the macro-pore exchange to the pressure head:
% dMacrPor    = zeros(P.nz,1); % must be developed!!!
% Sink term must be developed:
% sink        = zeros(P.nz,1); % must be developed!!!
%% pre-elaboration (in/out this function?)
% compute here or outside the function:
% I know that Antonio uses P.dt as current dt, but we should avoid to write
% user defined var. We simply update P.time(P.j) according to current
% simulation conditions and then calculate dt using the line below (it's
% advisable to write once in P.dt before to start new P.j):
i       = 1:P.nz;
i_i     = 2:P.nz-1;
i_im1   = 1:P.nz-2;
i_ip1   = 3:P.nz-0;
dz      = P.nodes.dz;
dz_im1  = P.nodes.dz(i_im1);
dz_i    = P.nodes.dz(i_i);
dz_ip1  = P.nodes.dz(i_ip1);
%% **IMPLICIT LINEARIZATION OF HYDRAULIC CONDUCTIVITIES
%   0   --> internodal conductivities refer to the old time level j-1
%   1   --> internodal conductivities refer to the new time level j
k       = 1;
%% **INITIAL CONDITIONS on h and teta to solve first iteration on p=1

% Initial pressure head at current j timestep
if P.j==1           % initial condition
    h_pm1       = P.hin;            %check that it is what you expect here!!
    pond        = 0;%   -I reset pond state variable in case of non-convergence
else
    h_pm1       = O.h22(:,P.j-1,mm);	%check that it is what you expect here!!
    pond        = O.pond(1,P.j-1,mm);%  -I reset pond state variable in case of non-convergence
end
% Volumic soil water content at former time level [-]
%   -this term is the same during the whole Newton-Raphson iteration
%    scheme (check if it is true!!)
teta_jm1= fnteta( h_pm1, P.sh, 1:P.nz ); % check with Antonio!!

% teta at first iteration p:
%   -I have to differentiate teta_jm1 from teta_p !!!!! How?
h_p             = h_pm1;
teta_p          = fnteta( h_p, P.sh, 1:P.nz ); % check with Antonio!!
% -------------------------------------------------------------------------
%% Nodal Hydraulic Conductivity [K] <-- h_pm1
% hydraulic conductivity, treated implictly (SWAP-32, Eq. 2.6, page
% 28):
% [Is it correct to treat the K in Eq. 2.6 as given below, for all "p"?]
% NOTES:
%   -if you transform ifc to a scalar, transform the
%    multilayer_conductivity_node to be parametric on it (so you can any
%    time decide whether to pass h o teta).
x               = teta_jm1; % unsure !!!!!
notTeta         = P.sh.ifc~=1 & P.sh.ifc~=3;
x(notTeta)      = h_pm1(notTeta);
K               = multilayer_conductivity_node( x, P.sh, i );
clear x notTeta;
%% root, drain, macropores [sink,macr]
% ** sink **
%      |--root--|   |--drain-|
sink = zeros(P.nz,1) + zeros(P.nz,1);
% **TEMPORARY SINK:
for i=1:P.nz
    if W.iveg==1 && W.itopvar==1
        if P.Droot>0
            P.dpt           = P.nodes.z(i);
            if W.iosm==1 && V.ifs>3
                P.op        = -P.ECstar(i)*360;
            else
                P.op        = 0;
            end
            if P.nodes.z(i) < P.Droot                    
                %P.sink(i)   = fnsink( P.h1(i), P, W, V );
                P.sink(i)   = fnsink_OLD(P.h1(i),P.dpt,P.op,P.Tp,V.hI,V.hII,V.hIIIH,V.hIIIL,V.hIV,V.hw50,V.pw1,V.hs50,V.ps1,V.aMH,V.bMH,P.Droot,V.ifs,V.rda,V.rdb,V.rdc,P.nz,V.ifg,V.zc,V.g0,V.gzc,V.Drf);
            else
                P.sink(i)   = 0;
            end
        else
            P.sink(i)       = 0;
        end
    end
end
%      |--MaPo--| i.e. macroporosity
macr = zeros(P.nz,1);
%% Internodal Hydraulic Conductivity [Kim2,Kip2] <-- K
multilayer_Kis2% requires: {K,dz}
%% top      boundary condition
% --- you should decrease time step to dtmin if required by boundtop
%     see SWAP-32, TimeControl(3):629
multilayer_boundtop
%% bottom   boundary condition -- put before while in wtus!!!
multilayer_boundbot
%% F-function [Fi #1]
multilayer_Fi% requires: {Kim2,Kip2, h_p, teta_p,teta_jm1}
% Fi(1)=Fi(2)-1.0d-2;
% Fi(end)=Fi(end-1)+1.0d-2;
% Estimate Fi inner product at init (sumFi2_0) & at "p"=0 (sumFi2_p):
sumFi2_0 = sum(Fi.^2)*0.5;
% To be able to enter while loop and use the full Newton step:
% sumFi2_p = sumFi2_0+1;
%% Maximum No of iterations, according to P.dt [pmax]
if P.dt == W.dtmin
    pmax = W.maxit*2;
else
    pmax = W.maxit;
end
%% **Newton-Raphson iteration scheme [loop on "p"]
for p=1:pmax
%% dF_dh
    % **UPDATE conditions at the current iteration-step (at "p")
    %   considering the actual lambda.
    %   -[dKi_dhi]  (y) --> Conductivity derivative to the pressure head
    %   -[Ci]       (y) --> Capillary Capacity
    %   -[teta_pm1] (n) --> Moisture fraction
    [dKi_dhi,Ci,~] = multilayer_dKi_dhi(h_pm1, P.sh, 1:P.nz, P.sh.k0);
    multilayer_dFdh% requires: {dKi_dhi,K, Ci, h_pm1, Kim2,Kip2}
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
%% BackTrack :: lambda [lambda=lambda/3]
    lambda      = 1;
    % backtrack block was broken:
    bt_breaked  = false;
    for iL = 1:W.maxbacktrack
%% factmax
        if P.dt == W.dtmin && p>W.maxit
            % use the factmax to update the Newton step:
            % ...to be developed
            error('Develop this piece of code!!')
            % h_p = ...
        else
%% Newton-Raphson step :: at current lambda
            % regular case
            h_p = h_pm1 - lambda*dh_p;
        end
%% theta
%         [~,Ci,teta_p] = multilayer_dKi_dhi(h_p, P.sh, 1:P.nz, P.sh.k0);
        % Moisture at "p" iteration:
        teta_p = fnteta( h_p, P.sh, 1:P.nz );
%% implicit [k=1]
        if k==1
%% Nodal Hydraulic conductivity [K]
            % IMPORTANT: I have to change ifc being a scalar and the
            % multilayer_conductivity_node must be parametric on it!!!!
            K = multilayer_conductivity_node( teta_p, P.sh, i );
%% Internodal Hydraulic Conductivity [Kim2,Kip2]
            multilayer_Kis2% requires: {K}
        end
%% top      boundary condition -- to be activated?
        % the actual pond here is the one calculated at p-1 (or outside
        % the p-loop if p=1)
        multilayer_boundtop
%% F-function [Fi #2]
        % The original headcalc.for/SWAP-32 implementation uses the same
        % piece of code for computing the F-function here (Fi*2) and before
        % the Newton-Raphson iteration scheme (which was called Fi*1).
        % The only differences (what Fi*2 has and not Fi*1) are at lines:
        %   -(193,195)
        %   +(445,450,452,453,454)
        multilayer_Fi% requires: {Kim2,Kip2, h_p, teta_p,teta_jm1}
        % Estimate Fi inner product at current "p" (sumFi_p):
        sumFi2_p    = sum(Fi.^2)*0.5;
        % *IMPORTANT
        %   -BackTrack block exit criterion:
        maxFi_p     = abs(max(Fi));
        %   -Newton-Raphson block exit criterion:
        sumFi_p     = sum(Fi); % [cm d-1]
%% break of BackTrack block
        if sumFi2_p<sumFi2_0 || maxFi_p<W.CritDevBalCp
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
%% break of Newton-Raphson block
    if abs(sumFi_p) <= W.CritDevBalTot
        nr_breaked = true;
        break
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
    % h at previous iteration-step "p-1" when I go to next iteration-step
    % "p":
    h_pm1       = h_p;
    % sum of Fi used in the backtrack block at step "p+1":
    sumFi2_0    = sumFi2_p;
end% p


if nr_breaked%  convergence has been reached!
    % *STORE current pressure head:
    O.h22(:,P.j,mm) = h_p;
    
    % *RISE TIMESTEP:   ***to be corrected***
    % To speedup simulation convergence
   if (p<=3)      ,     P.dt = min(P.dt*W.dtfactor,W.dtmax); end
   if (p>=W.maxit-2),   P.dt = max(P.dt*W.dtfactor,W.dtmin); end

% --- see TimeControl(3):612
else%           convergence could not be reached!
    my_dtfactor         = 3;
    n_noconv            = n_noconv +1;
    % *REDUCE TIMESTEP:
    % Pay attention to prevent very small dt values
    dtprevious          = P.dt;
    if P.dt > my_dtfactor*W.dtmin
        P.dt            = P.dt / my_dtfactor;
        P.time(P.j)     = P.time(P.j) -dtprevious +P.dt;
        % If we go back to the previous integer timestep, we must update:
        P.tidx          = floor(P.time(P.j))+1;

    elseif P.dt > W.dtmin
        P.dt            = W.dtmin;
        P.time(P.j)     = P.time(P.j) -dtprevious +P.dt;
        % If we go back to the previous integer timestep, we must update:
        P.tidx          = floor(P.time(P.j))+1;

    elseif P.dt==W.dtmin
        % *CONTINUE WITHOUT CONVERGENCE:
        % ...what to do?
        fl_noconv       = true;% is non-convergent?
    else
        error('Something wrong in setting P.dt!')
    end
%     P.iter(:,P.j) = [p;nr_breaked;fl_noconv;n_noconv;iL;bt_breaked];
    % --- reset soil state variables of time = j-1
    %h           = hm1;
    %theta       = thetm1;
    %kmean(numnod+1) = k(numnod);
    %gwl         = gwlm1;
end
%% clean
clear i i_i i_im1 i_ip1 dz_i dz dzp_ddg dz_im1 dz_ip1