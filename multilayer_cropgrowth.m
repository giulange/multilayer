%% develop a crop growth module
%% Roots depth
% P.Droot = [];
%% cropfixed ~SWAP
% % %     purpose            : simple crop growth routine for swap 
% % % ----------------------------------------------------------------------
% % 
% % %       parameter (nihil=1.0d-7)
% % 
% % 
% % 
% % % 1000  continue
% % 
% % % === initialization ===================================================
% % 
% % % --- read crop data
% %       gctb = 0.0d0
% %       cftb = 0.0d0
% %       chtb = 0.0d0
% %       rdtb = 0.0d0
% % %       call readcropfixed (cropfil(icrop),pathcrop,idev,lcc,tsumea,
% % %      & tsumam,tbase,kdif,kdir,gctb,swgc,cftb,swcf,rdtb,rdctb,hlim1,
% % %      & hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,kytb,
% % %      & cofab,logf,schedule,swinter,pfreetb,pstemtb,scanopytb,
% % %      & avprectb,avevaptb,cumdens,chtb,albedo,swetr,
% % %      & flsolute,ecmax,ecslop,c2eca,c2ecb,c2ecf,numlay,
% % %      & alphacrit,swroottyp,wiltpoint,rootradius,rootcoefa,rsw)           % NwRootExtr
% % 
% % % --- development stage
% %       dvs = 0.0d0
% % 
% % % --- initial lai or sc
% %       lai = afgen (gctb,(2*magrs),dvs)
% %       if (swgc==2) 
% %         gc = lai
% %         lai = lai*3.0d0
% %       end
% % 
% % % --- initial crop factor or crop height
% %       cf = afgen (cftb,(2*magrs),dvs)
% %       ch = afgen (chtb,(2*magrs),dvs)
% % 
% % % --- actual rooting depth [cm]
% %       rd = min (rds,afgen (rdtb,(2*magrs),dvs))
% % 
% % % --- initial summation variables of the crop
% %       cptr0 = 0.
% %       ctr0 = 0.
% % 
% % % --- init arrays with cum. pot. and act. transpiration 
% %       cptr = 0.0d0
% %       ctr = 0.0d0
% % 
% % % --- initialize matric flux potential                                  % NwRootExtr
% %       phead = 0.d0   % dummy                                            %
% %       if (swroottyp == 2)                                         %
% %         node = 1                                                        % dummy node nr
% % %         call MatricFlux(1,phead,node)                                   %
% % % ---   output of matric flux potential                                 %
% % %        if (swmfp==1) call outmatricflux(2,mfp,numnod,tcum,          %
% % %     &   mflux,z,outfil,pathwork,project,ptra,h)                       %
% %       end                                                             % NwRootExtr
% % 
% %       return          
% % 
% % % 2000  continue
% % 
% % % === calculate potential rate and state variables ======================
% % 
% % % 3000  continue
% % 
% % % === calculate actual rate and state variables ======================
% % 
% % % --- increase in temperature sum
% %       dtsum = max (0.0d0,tav-tbase)
% % 
% % % --- development rate
% %       if (idev==1) 
% %         dvr = 2.0/lcc
% %       elseif (idev==2) 
% %         if (dvs<1.0d0) 
% %           dvr = dtsum/tsumea
% %         else
% %           dvr = dtsum/tsumam
% %         end
% %       end
% % 
% % % --- determination of current growing stage
% %       for i = 1:magrs
% %         help(i) = kytb(2*i-1)
% %       end
% %       icgs = stepnr(help,magrs,dvs)
% % 
% % % --- water stress
% %       if(abs(ptra)<nihil) 
% %         reltr = 1.0d0
% %       else
% %         reltr = max(min(tra/ptra,1.0d0),0.0d0)
% %       end
% % 
% % % ----integrals of the crop --------------------------------------------
% % 
% % % --- phenological development stage
% %       dvs = dvs+dvr
% % 
% % % --- leaf area index or soil cover fraction    
% %       lai = afgen (gctb,(2*magrs),dvs)
% %       if (swgc==2) 
% %         gc = lai
% %         lai = lai*3.0d0
% %       end
% % 
% % % --- crop factor or crop height
% %       cf = afgen (cftb,(2*magrs),dvs)
% %       ch = afgen (chtb,(2*magrs),dvs)
% % 
% % % --- rooting depth [cm]
% %       rd = min (rds,afgen (rdtb,(2*magrs),dvs))
% % 
% % % --- cumulative relative transpiration, total growing season 
% %       cptr0 = cptr0 + ptra  
% %       ctr0 = ctr0  + tra
% %       if (cptr0.le.nihil) 
% %         crt0 = 0.0d0
% %       else
% %         crt0 = max(min(ctr0/cptr0,1.0d0),0.0d0)
% %       end
% % 
% % % --- cumulative relative transpiration, current growing stage
% %       cptr(icgs) = cptr(icgs) + ptra
% %       ctr(icgs) = ctr(icgs)  + tra
% %       if (cptr(icgs).le.nihil) 
% %         crt(icgs) = 0.0d0
% %       else
% %         crt(icgs) = max(min(ctr(icgs)/cptr(icgs),1.0d0),0.0d0)
% %       end
% % 
% % % --- relative yield per growing stage and cumulated
% %       crely = 1.0d0 
% %       for i = 1,icgs
% %         rely(i) = 1.0d0-((1.0d0-crt(i))*kytb(2*i))
% %         crely = crely*rely(i)
% %       end
% % 
% %       return