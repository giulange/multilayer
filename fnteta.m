function teta = fnteta( x, P, ii )
% teta = fnteta( x, P, ii )

df          = P.tetas(ii)-P.tetar(ii);

if x>0 
    teta    = P.tetas(ii);
elseif P.ifr(ii)==1
    em      = 1-1/P.en(ii);
    teta    = P.tetar(ii) + df/(1.+(abs(x)*P.alfvg(ii))^P.en(ii))^em;

elseif P.ifr(ii)==2
    em2     = 1-1/P.en2(ii);
    teta    = P.tetar(ii) + df * ( P.fi(ii) * ...
              ( (1+P.alfrs(ii)*abs(x))*exp(-P.alfrs(ii)*abs(x)) ) + ...
              (1-P.fi(ii))*(1+(P.alfvg2(ii)*abs(x))^P.en2(ii))^(-em2)    );    
%     lx=-1000000:1000:0;
%     ly=tetar+df.*(fi.*((1+alfrs.*abs(lx)).*exp(-alfrs.*abs(lx)))+(1-fi).*(1+(alfvg.*abs(lx)).^en).^(-em));

elseif P.ifr(ii)==3
    em      = 1-1/P.en(ii);
    em2     = 1-1/en2;
    sef1    = (1+(P.alfvg (ii)*abs(x))^P.en(ii))^-em;
    sef2    = (1+(P.alfvg2(ii)*abs(x))^P.en2(ii))^-em2;
    sef     = fi*sef1 + (1-fi)*sef2;
    teta    = P.tetar(ii) + df*sef;
end

return