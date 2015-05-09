%% define internal pars:
err_msg_wrong_par_set   = @(PAR) sprintf('This parameter was set in a wrong way: %s', PAR);
len_units_time          = @(ed,sd) (datenum(ed)-datenum(sd)+1)/W.timestep;
len_units_time_str1     = @(ed,sd) (datenum(ed,'yyyy-mm-dd,hh') - datenum(sd)+1)/W.timestep;
% len_units_time_str2     = @(ed,sd) (datenum(ed) - datenum(sd,'yyyy-mm-dd,hh')+1)/W.timestep;
%% check dirs
if ~exist( proj.ipath, 'dir' )
    error('The Project directory does not exist!')
end
if ~exist( proj.opath, 'dir' )
    mkdir( proj.opath )
elseif length(dir(proj.opath))-2>0
    ListOfFiles     = dir(proj.opath);
    ListOfFiles     = {ListOfFiles(3:end).name};
    warning('The output folder\n\t%s\nalready has this %d files:\n', ...
            proj.opath,length(ListOfFiles) )
    for ii=1:length(ListOfFiles)
        fprintf('\t> %s\n',ListOfFiles{ii})
    end
    
    STOP            = input('Type 0 to continue OR 1 to stop execution!\n');
    if STOP, error('You decided to stop the program'), end
    clear ListOfFiles ii STOP
end
%% W.tmax -- Simulation duration
% tmax is adapted according to size of W.timestep !!
if W.timestep<1/24 || W.timestep>10
    err_msg_wrong_par_set('W.timestep')
end
% if W.timestep~=1.00
%     error('A timestep different from "day" for any time-dependent data is not yet implemented!')
% end
W.tmax                  = len_units_time( W.edate, W.sdate );
%% n-iterations
% Nj:                        Number of iterations as a function 
%                               f(W.tmax,W.dtin,W.dtmax) according to the
%                               W.Nj_shp shape parameter.
P.Nj                    = round( W.tmax / (10^(W.Nj_shp + ...
                                  mean([log10(W.dtin),log10(W.dtmax)]) ...
                               )          )   );
%% adapt all time-related parameters affected by W.timestep
W.dtin          = W.dtin/W.timestep;
W.dtmin         = W.dtmin/W.timestep;
W.dtmax         = W.dtmax/W.timestep;
%% Read I_depth.txt -- to be implemented
%% Read I_time.txt -- to be implemented
%% adapt climatic data according to W.timestep
B.top.eto       = reshape(repmat(B.top.eto,1,1/W.timestep)',W.tmax,1);% [mm]
B.top.radi      = reshape(repmat(B.top.radi,1,1/W.timestep)',W.tmax,1);%[]
B.top.rain      = reshape(repmat(B.top.rain,1,1/W.timestep)',W.tmax,1);%[mm]
B.top.tmax      = reshape(repmat(B.top.tmax,1,1/W.timestep)',W.tmax,1);%[C]
B.top.tmin      = reshape(repmat(B.top.tmin,1,1/W.timestep)',W.tmax,1);%[C]
%% P.nz & SOIL GRID GEOMETRY
W.sg.sublayers_names        = {'SoilLay','SubLay','hSubLay','hNode','nNodes'};
switch W.sg.type
    case 1 % regular soil grid
        P.nz                = W.sg.regular(1);
        P.nodes             = multilayer_soilgrid(P.nz,W.zint);
    case 2 % sl --> sub-layers
        P.nz                = sum(W.sg.sublayers(:,5));
        P.nodes             = multilayer_soilgrid_sl(W.sg.sublayers,W.zint);
    otherwise
        % Not yet implemented
        err_msg_wrong_par_set( 'W.sg.type' )
end

% PLOT the soil grid geometry configuration:
if W.sg.plotme, multilayer_soilgrid_graph(P.nodes,W.zint); end
%% W.dtin<W.dtmin
if W.dtin<W.dtmin
    warning('W.dtin=%f < W.dtmin=%f',W.dtin,W.dtmin)
    W.dtin = W.dtmin;
end
%% P.tprint (print steps)
% This is the only time-dependent parameter whose LENGTH is not affected by
% the valorization of W.timestep, but ony the VALUES are.
if isempty(W.tp{1})
    % read from I_time.txt file.
    % P.tprint = VAR(:,col);
else
    if size(W.tp,1)==1
        P.tprint        = 1:W.tp{1,2}:W.tmax;
    else
        P.tprint        = [];
        for ii = 1:size(W.tp,1)-1
            P.tprint    = [P.tprint, datenum( W.tp{ii,1},'yyyy-mm-dd,hh' )+W.tp{ii,2}:W.tp{ii,2}:datenum( W.tp{ii+1,1},'yyyy-mm-dd,hh' )];
        end
        P.tprint        = [P.tprint, datenum( W.tp{end,1},'yyyy-mm-dd,hh' )+W.tp{end,2}:W.tp{end,2}:datenum( W.edate )];
        if P.tprint(end)~=datenum(W.edate)
            P.tprint(end+1) = datenum(W.edate);
        end
        P.tprint        = (P.tprint - datenum(W.sdate))/W.timestep;
    end
end
%% V.Kc (reduction coeff. of potential evapotraspiration)
if isempty(V.Kc{1})
    err_msg_wrong_par_set('V.Kc')
else
    if size(V.Kc,1)==1
        P.Kc            = repmat(V.Kc{1,2},1,W.tmax);
    else
        P.Kc            = [];
        for ii = 1:size(V.Kc,1)-1
            LEN         = len_units_time( datenum( V.Kc{ii+1,1},'yyyy-mm-dd,hh' ),datenum( V.Kc{ii,1},'yyyy-mm-dd,hh' ) );
            P.Kc        = [P.Kc, repmat(V.Kc{ii,2},1, LEN )];
        end
        LEN             = len_units_time( datenum( W.edate ), datenum( V.Kc{end,1},'yyyy-mm-dd,hh' ) );
        P.Kc            = [P.Kc, repmat(V.Kc{end,2},1, LEN )];
    end
end
%% V.Ke (reduction coeff. of potential soil evaporation)
if isempty(V.Ke{1})
    % calculate following FAO paper 56, Chapter 7
    err_msg_wrong_par_set('V.Ke')
else
    if size(V.Ke,1)==1
        P.Ke            = repmat(V.Ke{1,2},1,W.tmax);
    else
        P.Ke            = [];
        for ii = 1:size(V.Ke,1)-1
            LEN         = len_units_time( datenum( V.Ke{ii+1,1},'yyyy-mm-dd,hh' ),datenum( V.Ke{ii,1},'yyyy-mm-dd,hh' ) );
            P.Ke        = [P.Ke, repmat(V.Ke{ii,2},1, LEN )];
        end
        LEN             = len_units_time( datenum( W.edate ), datenum( V.Ke{end,1},'yyyy-mm-dd,hh' ) );
        P.Ke            = [P.Ke, repmat(V.Ke{end,2},1, LEN )];
    end
end
%% B.bot.hqstar (flux/potential at bottom boundary)
if isempty(V.Ke{1})
    % calculate following FAO paper 56, Chapter 7
    err_msg_wrong_par_set('B.bot.hqstars')
else
    if size(B.bot.hqstar,1)==1
        P.bothq         = repmat(B.bot.hqstar{1,2},1,W.tmax);
    else
        P.bothq         = [];
        for ii = 1:size(B.bot.hqstar,1)-1
            LEN         = len_units_time( datenum( B.bot.hqstar{ii+1,1},'yyyy-mm-dd,hh' ),datenum( B.bot.hqstar{ii,1},'yyyy-mm-dd,hh' ) );
            P.bothq     = [P.bothq, repmat(B.bot.hqstar{ii,2},1, LEN )];
        end
        LEN             = len_units_time( datenum( W.edate ), datenum( B.bot.hqstar{end,1},'yyyy-mm-dd,hh' ) );
        P.bothq         = [P.bothq, repmat(B.bot.hqstar{end,2},1, LEN )];
    end
end
%% S.compounds & S.fertilizers
if S.isfert
    P.compounds = NaN;
    
    % Then convert to something similar to P.compounds readly used by
    %    multilayer:
    % ...to be developed!!
    P.fertilizers = NaN;
    
else
    P.fertilizers = NaN;
    
    P.compounds = zeros( size(S.compounds,2)-1, W.tmax );
    for ii = 1:size(S.compounds,1)
        t_elem = len_units_time_str1( S.compounds{ii,1}, W.sdate );
        P.compounds(:,t_elem) = cell2mat(S.compounds(ii,2:end));
    end
    % Convert [mg cm-2] to [g cm-2] used in multilayer:
    P.compounds = P.compounds / 1000;
    
end
%% P.hin
P.hin       = multilayer_sub_valorization_depth( W.hin, P.nz, P.nodes.z(1:end-1) );
%% P.CDECinNH
P.CDECinNH  = multilayer_sub_valorization_depth( S.CDE.Cin.NH, P.nz, P.nodes.z(1:end-1) );
%% P.CDECinNO
P.CDECinNO  = multilayer_sub_valorization_depth( S.CDE.Cin.NO, P.nz, P.nodes.z(1:end-1) );
%% P.CDElambda
P.CDElambda = multilayer_sub_valorization_depth( S.CDE.lambda, P.nz, P.nodes.z(1:end-1) );
%% P.CDEKnitr
P.CDEKnitr = multilayer_sub_valorization_depth( S.CDE.Knitr, P.nz, P.nodes.z(1:end-1) );
%% P.CDEKimmob
P.CDEKimmob = multilayer_sub_valorization_depth( S.CDE.Kimmob, P.nz, P.nodes.z(1:end-1) );
%% P.CDEKdenitr
P.CDEKdenitr = multilayer_sub_valorization_depth( S.CDE.Kdenitr, P.nz, P.nodes.z(1:end-1) );
%% condizioni al contorno superiore e inferiore
% if (W.itopvar==0 && W.ibotvar==0) || (W.itopvar==1 && W.ibotvar==1)
%     error('W.itopvar=%d cannot be equal to W.ibotvar=%d',W.itopvar,W.ibotvar)
% end
%% hqstar -- flux/potential at top boundary
if W.itbc==0
    if ~B.top.isirri
        B.top.irri  = 0;
    end

    B.top.hqstar    = -(B.top.rain) -(B.top.irri);
end
%% itopvar -- it's set here and not in conf anymore
% itopvar:          Flag to set top boudary conditions:
%                       *0  --> fixed
%                               The value used is that of hsurf/qsurf
%                               (B.top.hqstar is ignored).
%                       *1  --> variable
%                               B.top.hqstar is used (hsurf/qsurf ignored).
if sum(abs(diff(B.top.hqstar))) > 0
    W.itopvar       = 1;% variable  top boudary condition
else
    W.itopvar       = 0;% fixed     top boudary condition
end
%% trasporto soluti
switch W.isol
    case 0
        % do nothing!

    case 1
        if ~W.iCtopvar==0, error('Something wrong!'), end

    case 2
        if ~W.iCtopvar==1, error('Something wrong!'), end
        %Con questo controllo si assume che il primo strato del profilo
        %corrisponda allo strato di interramento del concime azotato
        if and(B.Ctop.dL<W.zint(2),B.Ctop.dL>W.zint(2)) %# --> never TRUE!!!
            error('check B.Ctop.dL-W.zint(2)')
        end

end
%% vegetazione
if ~W.iveg==1, V = NaN; end
%% TEMP???????
B.Ctop.Tstar        = (B.top.tmax+B.top.tmin)/2;
%% montecarlo

if W.MTCL==0
    M.nvp               = 1;
elseif W.MTCL==1
    if M.nvp<2
        warning('\n\t%s (W.MTCL=%d),\n\t%s M.nvp=%d times!',  ...
            'You set stochastic simulation',W.MTCL, ...
            'but multilayer will run',M.nvp )
    end
    if sum(M.nlay) ~= size(M.data,1)
        error('%s\n  (M.nlay=%d, while M.data size is %d)', ...
            'The number of stochastic simulations is not properly defined',...
            sum(M.nlay), size(M.data,1))
    end
    if M.combinatorial==0
        if sum( min( M.nlay , M.nvp )~=M.nvp )
            M.bound     = min( M.nlay , M.nvp )~=M.nvp;
            fprintf('\nM.nvp = %d\n', M.nvp)
            warning('\nlayer={%d} is limited according to the value of M.nvp!\n', ...
                    find(M.bound)' )
            M.nvp_old   = M.nvp;
            M.nvp       = min( M.nlay );
            warning('A new value is set for M.nvp=%d (instead of %d)', ...
                    M.nvp,M.nvp_old)
        end

        % -----------------------------------------------------------------
        % compute the combinations using the first M.nvp stochastic
        % simulations from each layer in sequence:
        % -----------------------------------------------------------------
        M.combinations = NaN(M.nvp,W.nlay);
        for ii = 1:W.nlay
            if ii==1
                M.combinations(:,1) = 1:M.nvp;
            else
                s=sum(M.nlay(1:ii-1))+1;
                M.combinations(:,ii) = s:s+M.nvp-1;
            end
            
        end
        clear ii s
        % -----------------------------------------------------------------
        
    elseif M.combinatorial==1
        M.nvp_old       = M.nvp;
        M.nvp           = prod( M.nlay );
        
        % -----------------------------------------------------------------
        % compute all possible combinations from W.nlay sets of stochastic
        % simulations:
        % -----------------------------------------------------------------
        varin = cell( 1,W.nlay );
        for ii = 1:W.nlay
            if ii==1, start=1; else start=sum(M.nlay(1:ii-1))+1; end
            varin{ii}=start:start+M.nlay(ii)-1;
        end
        % replicate the grid vectors to compute the full grid:
        [tmp{1:W.nlay}]     = ndgrid( varin{1:3} );
        % concatenate
        M.combinations      = reshape(cat(W.nlay+1,tmp{:}),[],W.nlay);
        clear varin ii start tmp
        % -----------------------------------------------------------------

    else error( err_msg_wrong_par_set( 'W.combinatorial' ) );
    end
else
    error( err_msg_wrong_par_set( 'W.MTCL' ) );
end
if (W.MTCL==0 && M.nvp>1) || (W.MTCL==1 && M.nvp<2)
    error( err_msg_wrong_par_set( 'M.nvp' ) );
end
%% cutted code:

% %% condizioni al contorno superiore e inferiore ==> mettere in _conf.m
% if W.itopvar==1 && W.ibotvar==0
%     B.sel.name          = 'top';
%     B.sel.tqstar        = B.top.thqstar;
%     B.sel.qstar         = B.top.hqstar;
% elseif W.itopvar==0 && W.ibotvar==1
%     B.sel.name          = 'bot';
%     B.sel.tqstar        = B.bot.thqstar;
%     B.sel.qstar         = B.bot.hqstar;
% else % if they are [0,0] OR [1,1]
%     % see multilayer_checkpars.m
% end
% 
% %% trasporto soluti ==> mettere in _conf.m
% switch W.isol
%     case 0
%         % do nothing!
%     case 1
%         % do nothing!
%     case 2 % ==> check if that is what you were thinking about!!
%         B.sel.name      = 'Ctop';
%         B.sel.tqstar    = B.Ctop.tqstar;
%         B.sel.qstar     = NaN; % B.Ctop.qstar does not exist!
% end
% 
% %% vegetazione ==> mettere in _conf.m, solo se non � inutile!!
% if W.iveg==1 % ==> inutile, meglio definire un unico tqstar per tutto il programma!! 
%     B.sel.name          = 'veg';
%     B.sel.tqstar        = V.tqstar;
% %     B.sel.qstar         = B.bot.hqstar; % ??non c'�??
% end
%% ...something else?
%% include checks on EC (in particular EC.matrix!!)
%% end
clear err_msg_wrong_par_set len_units_time ii LEN len_units_time_str1 len_units_time_str2