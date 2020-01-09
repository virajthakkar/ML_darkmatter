x_lo=0;
x_hi=160;
y_lo=0;
y_hi=50;
x=[0:0.01:100];



mass_dm=100;
modi_sigmoid=0.9*sigmoid(x-2);
[hax, h1, h2] = plotyy (x, ratefun(x,mass_dm), x, modi_sigmoid,'plot');

set(h1(1),'LineStyle','-','Marker','+','MarkerSize',4,'MarkerEdgeColor','b','Linewidth',2);
%set(h1(2),'LineStyle','none','Marker','o','MarkerSize',6,'MarkerEdgeColor','b');  
set(h2(1),'LineStyle','-','Marker','+','MarkerSize',3,'MarkerEdgeColor','r','Linewidth',2);

mul= modi_sigmoid.*ratefun(x,mass_dm);
hold on
plot(x,mul,'k','Linewidth',2,'Markersize',7)
title("WIMP Model",'FontSize',18,'FontWeight','bold');
xlabel(hax(1),"E_R ( keV)",'FontSize',18,'FontWeight','bold');
 ylabel (hax(1), "WIMP spectrum",'FontSize',18,'FontWeight','bold');
 ylabel (hax(2), "Detector eff.");

 hold on;
%figure 2;
load Cf_data.txt;
[n,p]=size(Cf_data)
A=Cf_data;

xCf=A(:,1);
yCf=0.75*power(10,-7)*A(:,2)/1.68355;
yerr=A(:,3);

plot(xCf,yCf,'g','Linewidth',2,'Markersize',7);
 
%l=legend(sprintf('m=%d GeV/c^2 ',mass_dm) );
l=legend(sprintf('WIMP spectrum m=%d GeV/c^2 ',mass_dm),"Corrected spectrum","Cf data","Det eff(Sigmoid)")

set (l, "fontsize", 14) 


figure 2;

xCf2=A(10:425,1);
yCf2=0.75*power(10,-7)*A(10:425,2)/1.68355;
yw=0.9*sigmoid(xCf2-2).*ratefun(xCf2,mass_dm)./yCf2;  %weight vector

plot(xCf2,yw,'o','Linewidth',2,'Markersize',6)

title("Weight vector distribution",'FontSize',18,'FontWeight','bold');

xlabel({'E_R ','( keV)'},'FontSize',18,'FontWeight','bold');
ylabel({'Weight values',''},'FontSize',18,'FontWeight','bold');


yCf_weighted=yCf2.*yw;  % reweighting Cf data

figure 3;

plot(xCf2,yCf_weighted,'k','Linewidth',2,'LineStyle','-','Marker','+','Markersize',4);
hold on;
plot(xCf2,yCf2,'g','Linewidth',2,'Markersize',7);

hold on;
plot(xCf2,yw.*yCf2,'r','Linewidth',1);
title("Corrected RRQs obtained from weight vector",'FontSize',18,'FontWeight','bold');

xlabel({'E_R ','( keV)'},'FontSize',18,'FontWeight','bold');
ylabel({'Rate','( units)'},'FontSize',18,'FontWeight','bold');

l2=legend("weight vector RRQ","Cf data",sprintf('theory corrected m=%d GeV/c^2 ',mass_dm))

set (l, "fontsize", 14) 


