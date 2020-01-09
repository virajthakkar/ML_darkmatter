

x_lo=0;
x_hi=160;
y_lo=0;
y_hi=20;
x=[0:0.01:160];

mass_dm=100;
figure 1;
plot(x,ratefun(x,mass_dm),'r','Linewidth',2,'Markersize',7)

grid("on");
xlabel({'E_R ','( keV)'},'FontSize',18,'FontWeight','bold');
ylabel({'Rate','( units)'},'FontSize',18,'FontWeight','bold');
title("Plot of Event rate of DM-nucleon interaction vs recoil energy",'FontSize',18,'FontWeight','bold');
limits= axis ([x_lo x_hi])
%hold on;
modi_sigmoid=0.9*sigmoid(x-2);
%plot(x,modi_sigmoid,'b','Linewidth',2,'Markersize',7)

mul= modi_sigmoid.*ratefun(x,mass_dm);
hold on
plot(x,mul,'k','Linewidth',2,'Markersize',7)


%l=legend(sprintf('m=%d GeV/c^2 ',mass_dm),sprintf('m=%d GeV/c^2 ',5*mass_dm),sprintf('m=%d GeV/c^2 ',10*mass_dm),
%sprintf('m=%d GeV/c^2 ',15*mass_dm)  );





hold on;
%figure 2;
load Cf_data.txt;
[n,p]=size(Cf_data)
A=Cf_data;

xCf=A(:,1);
yCf=0.75*power(10,-7)*A(:,2)/1.68355;
yerr=A(:,3);

plot(xCf,yCf,'g','Linewidth',2,'Markersize',7);
%hold on;
%errorbar(xCf,yCf,yerr)



l=legend(sprintf('theory dark matter m=%d GeV/c^2 ',mass_dm),"detector eff corrected(used sigmoid)","Cf data")

set (l, "fontsize", 14) 




figure 2;

xCf2=A(10:425,1);
yCf2=0.75*power(10,-7)*A(10:425,2)/1.68355;
yw=0.9*sigmoid(xCf2-2).*ratefun(xCf2,mass_dm)./yCf2;

plot(xCf2,yw,'o','Linewidth',2,'Markersize',6)

title("Weights calculated",'FontSize',18,'FontWeight','bold');

xlabel({'E_R ','( keV)'},'FontSize',18,'FontWeight','bold');
ylabel({'Weights','( units)'},'FontSize',18,'FontWeight','bold');
