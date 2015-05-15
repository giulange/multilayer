% --- update hydraulic conductivities to time level t+1
%       do i = 1,numnod
%          k(i) = hconduc(i,theta(i),cofgen(1,i),swfrost,rfcp(i),
%      &                  swsophy,numtab,sptab)
%          if(swmacro.eq.1)  k(i) = FrArMtrx(i) * k(i)
%          if(i.gt.1)then
%             kmean(i) = hcomean(swkmean,k(i-1),k(i),dz(i-1),dz(i))
%          end if
%       enddo
%       kmean(numnod+1) = k(numnod)
%% Nodal/Internodal Hydraulic Conductivity [P.K,P.Kim2,P.Kip2] <-- h,teta
multilayer_Kis2

% --- calculate actual water content of profile
%       call watstor (volm1,volact,numnod,theta,dz,FrArMtrx)
% multilayer_waterstorage

% --- calculate water fluxes between soil compartments
%       call fluxes (q,qbot,dt,inq,numnod,thetm1,theta,dz,qrot,qdra,
%      &  qimmob,qtop,qrosum,qdrtot,volact,volm1,swbotb,FrArMtrx,
%      &  QExcMpMtx,QMaPo,nrlevs,fllowgwl)
multilayer_fluxes

% --- calculation of states macropores and intermediate & cumulative values
%       if (flMacroPore) call macropore(4)
if W.SwMacro multilayer_macropore, end

% --- calculate cumulative fluxes
%       call integral 
% multilayer_integral

% --- update parameters for soil water hystereses
%       if (swhyst.ne.0) 
%      &  call hysteresis (numnod,layer,h,hm1,indeks,tau,paramvg,cofgen,
%      &                      dimoca,theta,disnod,dt,swsophy,numtab,sptab)
% multilayer_hysteresis
