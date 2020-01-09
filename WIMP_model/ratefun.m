function y=ratefun(x,mx)
%  mx is dark matter mass
mt=63.03655288; % Mass of target =0.932A in GeV/c^2
  Er=2*mt*5885.6226*power(10,-10)*mx*mx/(mx+mt)^2;   %v0=230 km/s
  R0=power(10,-8)*503/(mx*mt);
  y=(R0)*(1/Er)*exp(-x*power(10,-6)/Er);
  %x is in keV,but we need to put it in GeV scale, so 10^-6 factor in exponential
  
endfunction
