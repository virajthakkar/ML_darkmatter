
clear all;
load efficiency.txt;

[n,p]=size(efficiency);
eff=efficiency;

x_eff=eff(:,1);
y_eff=eff(:,2);

lerr_eff=eff(:,3); #lower error 
uerr_eff=eff(:,4);

mass_dm=100;

load Cf_data.txt;
[n,p]=size(Cf_data)
A=Cf_data;


figure 1;
[hax, h1, h2] = plotyy (x_eff, ratefun(x_eff,mass_dm), x_eff, y_eff,'plot');

set(h1(1),'LineStyle','-','Marker','+','MarkerSize',4,'MarkerEdgeColor','b','Linewidth',2);
%set(h1(2),'LineStyle','none','Marker','o','MarkerSize',6,'MarkerEdgeColor','b');  
set(h2(1),'LineStyle','-','Marker','+','MarkerSize',3,'MarkerEdgeColor','r','Linewidth',2);
title("WIMP Model",'FontSize',18,'FontWeight','bold');
xlabel(hax(1),"E_R ( keV)",'FontSize',18,'FontWeight','bold');
 ylabel (hax(1), "WIMP spectrum",'FontSize',18,'FontWeight','bold');
 ylabel (hax(2), "Detector eff.");

# multiply theory rate X eff


 
mul= y_eff.*ratefun(x_eff,mass_dm);
hold on
plot(x_eff,mul,'k','Linewidth',2,'Markersize',7)

 hold on;
%figure 2;

xCf=A(:,1);
yCf=0.75*power(10,-7)*A(:,2)/1.68355;
yerr=A(:,3);  # counts of Cf data

plot(xCf,yCf,'g','Linewidth',2,'Markersize',7);

l=legend(sprintf('WIMP spectrum m=%d GeV/c^2 ',mass_dm),"Corrected spectrum","Cf data","Det eff")

set (l, "fontsize", 14) 



#figure 2;
#plot(x_eff,y_eff,'k','Linewidth',2,'Marker','+','Markersize',7)
#hold on;
#figure 3;
#plot(xCf2,det_eff_inter,'r','Linewidth',2,'Marker','*','Markersize',2)

xCf2=A(10:425,1);
det_eff_inter=interp1(x_eff,y_eff,xCf2,'spline');  #interpolated detector efficiency on points of xCf
yCf2=yCf(10:425,1);
mul_inter=interp1(x_eff,mul,xCf2,'spline');

figure 2; plot(x_eff,mul,'r','Linewidth',2,'Marker','*','Markersize',2)
hold on; 
plot(xCf2,mul_inter,'k','LineStyle','-','Linewidth',1,'Marker','+','Markersize',4)
l=legend(sprintf('Corrected spectrum for m=%d GeV/c^2',mass_dm),"Interpolated corrected spectrum at Cf data points")
grid("on");
xlabel({'E_R ','( keV)'},'FontSize',18,'FontWeight','bold');
ylabel({'Rate','( units)'},'FontSize',18,'FontWeight','bold');
title("Corrected spectrum",'FontSize',18,'FontWeight','bold');
#limits= axis ([0 150])
xlim([0 150])
ylim([0 inf])

figure 3;

yw=mul_inter./yCf2;

plot(xCf2,yw,'o','Linewidth',2,'Markersize',5)

title("Weight vector distribution",'FontSize',18,'FontWeight','bold');

xlabel({'E_R ','( keV)'},'FontSize',18,'FontWeight','bold');
ylabel({'Weight values',''},'FontSize',18,'FontWeight','bold');
xlim([0 150])
ylim([0 inf])


