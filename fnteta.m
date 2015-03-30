function teta = fnteta( x, Psh, ii )
% teta = fnteta( x, P.sh, ii )

df          = Psh.tetas(ii)-Psh.tetar(ii);

if x>0 
    teta    = Psh.tetas(ii);
elseif Psh.ifr(ii)==1
    em      = 1-1/Psh.en(ii);
    teta    = Psh.tetar(ii) + df/(1.+(abs(x)*Psh.alfvg(ii))^Psh.en(ii))^em;

elseif Psh.ifr(ii)==2
    em2     = 1-1/Psh.en2(ii);
    teta    = Psh.tetar(ii) + df * ( Psh.fi(ii) * ...
              ( (1+Psh.alfrs(ii)*abs(x))*exp(-Psh.alfrs(ii)*abs(x)) ) + ...
              (1-Psh.fi(ii))*(1+(Psh.alfvg2(ii)*abs(x))^Psh.en2(ii))^(-em2)    );    
%     lx=-1000000:1000:0;
%     ly=tetar+df.*(fi.*((1+alfrs.*abs(lx)).*exp(-alfrs.*abs(lx)))+(1-fi).*(1+(alfvg.*abs(lx)).^en).^(-em));

elseif Psh.ifr(ii)==3
    em      = 1-1/Psh.en(ii);
    em2     = 1-1/en2;
    sef1    = (1+(Psh.alfvg (ii)*abs(x))^Psh.en(ii))^-em;
    sef2    = (1+(Psh.alfvg2(ii)*abs(x))^Psh.en2(ii))^-em2;
    sef     = fi*sef1 + (1-fi)*sef2;
    teta    = Psh.tetar(ii) + df*sef;
end

return