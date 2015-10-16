%% INITIALIZE SOIL WATER RATE/STATE VARIABLES -- Swap::SoilWater(1)
%% ! --- initialize miscellaneous
%   -pressure head of air near the soil surface [see soilwater.for]:
%   - you can implement it's calculation as done in multilayer by Coppola!
hatm            = -2.75d+05;%   [cm??]
P.nraidt        = 0.0d0; % net precipitation flux during iteration timesteps [cm d-1] 
P.nird          = 0.0d0; % net irrigation depth [cm]
%       if (swinco.ne.3) then
%         ldwet = 0.d0
%         spev = 0.d0
%       endif
runon           = 0.0d0;
%       saev = 0.d0
P.qtop          = 0.d0; % water flux through soil surface [cm d-1]
P.q             = zeros(P.nz+1,1);% it was P.flux for simplicity!!
%       do i = 1,macp
%         evp(i) = 0.0d0
%         rfcp(i) = 1.0d0
%       enddo
%       vtair = 0.0d0
%       cQMpLatSs = 0.0d0
%% ! --- Soil physics: tabulated or MualemVanGenuchten functions
%       cofgen = 0.0d0
%       if(swsophy.eq.1) then
% !       tabulated functions (h,theta,k,dthetadh,dkdtheta) tabulated
%         do node = 1,numnod
%           numtab(node) = numtablay(layer(node))
%         end do
%         do node = 1,numnod
%           do i = 1,5
%             do j = 1, numtab(node)
%               sptab(i,node,j) = sptablay(i,layer(node),j)
%             end do
%           end do
% !         assign values to cofgen
%           cofgen(1,node) = 0.0d0                          ! thetar
%           cofgen(2,node) = sptab(2,node,numtab(node))     ! thetas
%           cofgen(3,node) = sptab(3,node,numtab(node))     ! ksat
%         end do
%         do lay = 1,numlay
%           ksatfit(lay) = cofgen(3,nod1lay(lay))
%           thetsl(lay) = cofgen(2,nod1lay(lay))
%         end do
%       else
% !       MvanG functions 
%         do node = 1,numnod
%           lay = layer(node)
%           do i = 1,10
%             cofgen(i,node) = paramvg(i,lay)
%           end do
% !         assign dummy value to alphaw
%           cofgen(8,node) = -9999.9d0
%           cofgen(11,node) = relsatthr(lay)
%           cofgen(12,node) = ksatthr(lay)
%         end do
%         thetsl = 0.0d0
%         do lay = 1, numlay
%           thetsl(lay) = paramvg(2,lay)
%         end do
%       endif
% 
%% ! --- saturated and residual watercontent of each node; hysteresis parameters
%       do node = 1,numnod
%         lay = layer(node)
%         thetar(node) = cofgen(1,node)
%         thetas(node) = cofgen(2,node)
% !!! Kroes: disable combi of swsophy=1 and swhyst=1
%         if (swhyst.eq.1) then
% ! ---     wetting curve  
%           indeks(node) = 1
%           cofgen(4,node) = paramvg(8,lay)
%         elseif (swhyst.eq.0.or.swhyst.eq.2) then
% ! ---     drying branch or simulation without hysteresis
%           indeks(node) = -1
%           cofgen(4,node) = paramvg(4,lay)
%         endif
%       end do
% 
%% ! --- additional Input checks
% !   - ThetCrMP should be less then CofGen(2,Lay): requires additional inputcheck
%       if (swmacro.eq.1) then
%          do lay = 1, NumLay
%             if (SwSoilShr(lay).gt.0 .and. 
%      &          ThetCrMp(lay).gt.CofGen(2,nod1lay(lay))) then
%                messag = ' ThetCrMP.gt.ThetSat'
%                call fatalerr('MacroRead',messag)
%             endif
%          enddo
%       endif
%% h(i)
P.h             = P.hin;
%       if (swinco.eq.1) then
% ! --- pressure head profile is input 
%         do i = 1, nhead
%           tab(i*2) = h(i)
%           tab(i*2-1) = abs(zi(i))
%         end do
%         do i = 1, numnod
%           h(i) = afgen(tab,macp*2,abs(z(i)))
%         end do
%       endif
%       if (swinco.eq.2 .and. swbotb.ne.8) then
%         if (abs(gwli-(z(numnod)-0.5d0*dz(numnod))) .lt.1.0d-4) then
%           messag = 'Initial groundwaterlevel (SWINCO=2) is '//
%      &    'too close to bottom of soil profile'//
%      &    ' must be corrected!'
%           call fatalerr ('soilwater',messag)
%         endif
%       endif
%       if (swinco.eq.3) then
%         if (nhead.ne.numnod) then
%           messag = 'Initial data are read from file (SWINCO=3) and '//
%      &    'number of nodes/compartments is not consistent with NUMNOD'//
%      &    'must be corrected!'
%           call fatalerr ('soilwater',messag)
%         endif
%       endif
%       if (swinco.eq.1.or.swinco.eq.3) then
% ! ---   determine groundwater level
%         if (h(numnod) .gt. -1.d-5) then
%           i = numnod
%           do while ((h(i) .gt. -1.d-5) .and. (i .gt. 1))
%               i = i - 1
%           end do
%           if (h(i) .lt. -1.d-5) then 
%             gwl = z(i+1) + h(i+1) / (h(i+1) - h(i)) * (z(i) - z(i+1))
% ! ---       assume hydrostatic equilibrium in saturated part
%             do j = i+1, numnod
%               h(j) = gwl - z(j)
%             end do
%           endif
%         endif
%       else
%% ! --- pressure head profile is calculated from groundwater level
%         if (swbotb.eq.1) then  
%           gwl = afgen (gwltab,mabbc*2,t1900+dt-1.d0)
% 
%           if(abs(gwl-(z(numnod)-0.5d0*dz(numnod))) .lt.1.0d-4) then
%           messag = 'Groundwaterlevel as bottom boundary (SWBOTB=1) is'//
%      &    'below or to close to bottom of soil profile'//
%      &    ' must be corrected!'
%             call fatalerr ('soilwater',messag)
%           endif
%         else
%           gwl = gwli
%         endif
%         if (gwl.gt.0.0d0) then 
%           pond = gwl
%         else
          pond = 0.0d0;
%         endif
%         do i = 1,numnod
%           h(i) = gwl - z(i)
%         end do
%       endif
% 
%% ! --- in case of preferential flow, adjust Van Genuchten parameters
%       do i = 1, numnod
%         theta(i) = watcon(i,h(i),cofgen(1,i),swsophy,numtab,sptab)
%       end do
% P.teta          = fnteta( P.hin, P.sh, 1:P.nz );
P.teta          = multilayer_fnteta_vgm( P.hin, P.sh, 1:P.nz );
%% ! --- initial soil water storage
if ~W.SwMacro
    %  Fraction of horizontal area of soil matrix per compartment (-):
    P.FrArMtrx(1:P.nz) = 1.d0;
% --- TO BE ACTIVATED
%     volact = 0.0d0;
%     watstor(volm1,volact,numnod,theta,dz,FrArMtrx)
%     volini = volact
end
% pondini = pond;
%% ! --- hydr. conductivities, differential moisture capacities
% !     and mean hydraulic conductivities for each node
%       do node = 1,numnod
%         dimoca(node) = moiscap(node,h(node),cofgen(1,node),dt,
%      &                                             swsophy,numtab,sptab)
% 
%         FrArMtrx(node) = 1.d0
%         k(node) = hconduc (node,theta(node),cofgen(1,node), 
%      &                     swfrost,rfcp(node),swsophy,numtab,sptab)
%         if(swmacro.eq.1)  k(node) = FrArMtrx(node) * k(node)
% 
%         if(node.gt.1) kmean(node) =  hcomean(swkmean,k(node-1),k(node),
%      &                                    dz(node-1),dz(node))
%       end do
%       kmean(numnod+1) = k(numnod)
P.cap           = multilayer_capacity( P.hin,  P.sh, 1:P.nz, P.dt );
% Nodal/Internodal Hydraulic Conductivity [P.K,P.Kim2,P.Kip2] <-- h,teta
multilayer_Kis2
%% ! --- initial groundwater level
%       call calcgwl(logf,swscre,swbotb,flmacropore,numnod,CritUndSatVol,
%      &             gwlinp,h,z,dz,pond,t1900,Theta,ThetaS,gwl,
%      &             nodgwl,bpegwl,npegwl,pegwl,nodgwlflcpzo,gwlflcpzo)
% 
%       return