workDirectory='/home/eidelyur/mcool/eic';
filename='/home/eidelyur/mcool/eic/frictionForceLong_Nsplit_5-6-7-9.dat';
% filename='/home/eidelyur/mcool/eic/frictionForceLong_Nsplit_2-3-4-5.dat';
cd (workDirectory);
workDirectory='./';
fileID=fopen(filename,'r')


for i=1:13
    tline = fgetl(fileID)
end

k=0;
while ~feof(fileID)
    k=k+1
    tline = fgetl(fileID)
    if tline ~= -1
        [B,count] = fscanf(fileID, '%e %e %e %e %e',inf);
    end
     display(['count=',num2str(count)]);
end
fclose(fileID)
count

lines=count/5-1

Vion=zeros(lines+1,1);
FrctnForce=zeros(lines+1,4);

n=1;
for k=1:lines+1
    Vion(k,1)=B(n);                %
    FrctnForce(k,1)=B(n+1);               %
    FrctnForce(k,2)=B(n+2);               %
    FrctnForce(k,3)=B(n+3);               %
    FrctnForce(k,4)=B(n+4);               %
    n=n+5;
end

Vion
FrctnForce

Nsplit= [9,7,6,5]';
Delta_e_par = 1.e5;
powDelta_e_par = round(log10(Delta_e_par))
mantDelta_e_par= Delta_e_par/(10^powDelta_e_par)

figure(110)
nBeg=1;
nFin=20;
plot(Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,1),'-b', ...
     Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,2),'-m', ...
     Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,3),'-g', ...
     Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,4),'-k', ...
     Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,1),'xr', ...
     Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,2),'xr', ...
  	 Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,3),'xr', ...
     Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,4),'xr','LineWidth',2) 
xlabel('Relative Ion Velocity, V_{i||}/v_{e||}','Color','m','FontSize',14)
ylabel('-F_{||}, eV/m','Color','m','FontSize',14)
title(['Longitudinal Friction Force: V_{e||}=',num2str(mantDelta_e_par,'%3.1f'), ...
       '\cdot10^',int2str(powDelta_e_par),'m/s'],'Color','m','Fontsize',14)
legend({['N_{split}=',int2str(Nsplit(1))],['N_{split}=',int2str(Nsplit(2))], ...
        ['N_{split}=',int2str(Nsplit(3))],['N_{split}=',int2str(Nsplit(4))]}, ...
       'Location','northeast','Fontsize',14)
grid on

figure(120)
nBeg=7;
nFin=20;
loglog(Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,1),'-b', ...
       Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,2),'-m', ...
       Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,3),'-g', ...
       Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,4),'-k', ...
       Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,1),'xr', ...
       Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,2),'xr', ...
  	   Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,3),'xr', ...
       Vion(nBeg:nFin), -FrctnForce(nBeg:nFin,4),'xr','LineWidth',2) 
xlabel('Relative Ion Velocity, V_{i||}/v_{e||}','Color','m','FontSize',14)
ylabel('-F_{||}, eV/m','Color','m','FontSize',14)
xlim([.55,2.])
ylim([150.,5000.])
title(['Longitudinal Friction Force: V_{e||}=',num2str(mantDelta_e_par,'%3.1f'), ...
       '\cdot10^',int2str(powDelta_e_par),'m/s'],'Color','m','Fontsize',14)
legend({['N_{split}=',int2str(Nsplit(1))],['N_{split}=',int2str(Nsplit(2))], ...
        ['N_{split}=',int2str(Nsplit(3))],['N_{split}=',int2str(Nsplit(4))]}, ...
       'Location','southwest','Fontsize',14)
grid on

FrctnForce_1=zeros(14,1);
VionShort=zeros(14,1);
sumV=0.;
sumF=0.;
sumV2=0.;
sumFV=0.;
for k=1:14
   FrctnForce_1(k)=-1./FrctnForce(k+6,1);
   VionShort(k)=Vion(k+6);
   sumV=sumV+log10(VionShort(k));
   sumV2=sumV2+log10(VionShort(k))^2;
   sumF=sumF+log10(FrctnForce_1(k));
   sumFV=sumFV+log10(FrctnForce_1(k))*log10(VionShort(k));
end

delta=14*sumV2-sumV^2;
A=(sumV2*sumF-sumV*sumFV)/delta
B=(14*sumFV-sumV*sumF)/delta

FrctnForce_1_fit=zeros(14,1);
for k=1:14
   FrctnForce_1_fit(k)=A+B*log10(VionShort(k));  
end

figure(130)
nBeg=7;
nFin=20;
% plot(log10(VionShort),log10(FrctnForce_1),'-r',log10(VionShort),FrctnForce_1_fit,'xb','LineWidth',2) 
plot(log10(VionShort),log10(FrctnForce_1),'-r','LineWidth',2) 
xlabel('Relative Ion Velocity, log_{10}(V_{i||}/v_{e||})','Color','m','FontSize',14)
ylabel('log_{10}(-F_{||}^{-1}), eV/m','Color','m','FontSize',14)
% xlim([.55,2.])
% ylim([2.e-4,6.e-3])
% titleHeader = ('Longitudinal Friction Force: V_{e||}=',num2str(mantDelta_e_par,'%3.1f'),%3.1f\cdot10^{%2d} m/s',  ...
%               (mantDelta_e_par,powDelta_e_par))
title(['Longitudinal Friction Force: V_{e||}=',num2str(mantDelta_e_par,'%3.1f'), ...
       '\cdot10^',int2str(powDelta_e_par),'m/s'],'Color','m','Fontsize',14)
grid on

