function R=erfratefun(Er,mx)  %Er is in keV ; 1 keV=1.609 X 10^(-16) Joules
%  mx is dark matter mass
% A=72.64 for Ge
mt=63.062597376; % Mass of target =0.932A in GeV/c^2
  %Er=2*mt*5885.6226*power(10,-10)*mx*mx/(mx+mt)^2;   %v0=230 km/s
  
  v0=220; %most probable WIMP velocity in km/s
  vE=232; % mean circular velocity of Earth through the halo
  vesc=544; % Galactic escape velocity
 
 uT=mt*mx/(mt+mx);  % units of mass in GeV/c2   % 1 GeV/c2=1.79 X 10^(-27) kg ; 1 keV/c2=1.79 X 10^(-33) kg 
 scalefactor=sqrt(mt/2.0)*(2.998137224*10^2)*(1/uT);
 vmin=Er.^(0.5);
 vmin=scalefactor.*vmin;
 [m,n]=size(vmin)
 k0byk1=1/0.993361484; 
  R=zeros(m,n);
  for i=1:n
  x=vmin(1,i)/v0;
  y=vE/v0;
  z=vesc/v0;
  if x>=0 && x<z-y
  R(1,i)=erf(x+y)-erf(x-y)-(4/sqrt(pi))*y*exp(-z^2); endif
  if x>z-y && x<z+y
  R(1,i)=erf(z)-erf(x-y)-(2/sqrt(pi))*(y+z-x)*exp(-z^2); endif
  if x>y+z
  R(1,i)=0; endif
  i
  fprintf('Value of vmin/v0 \n')
  vmin(1,i)/v0
endfor




%R0=power(10,-8)*503/(mx*mt);
%y=(R0)*(1/Er)*exp(-x*power(10,-6)/Er);
%x is in keV,but we need to put it in GeV scale, so 10^-6 factor in exponential
  
endfunction
