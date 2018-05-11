%========================================================
%
% This code compairs two approaches: "classical" (from [1]) and
% "magnus" (from [2]).
%
% For "classical" approach the magnetized interaction between ion
% and electron is considered for ion velocities V_i > rmsTrnsvVe.
%
% References:
%
% [1] I. Meshkov, A. Sidorin, A. Smirnov, G. Trubnikov.
%    "Physics guide of BETACOOL code. Version 1.1". C-A/AP/#262, November
%    2006, Brookhaven National Laboratory, Upton, NY 11973.
% [2] 
%
%========================================================

pi=3.14159265358;

%
% Physical constants:
%
m_e=9.10938356e-28;         % electron mass, g
m_p=1.672621898e-24;        % electron mass, g
q_e=4.803204673e-10;        % electron charge, sqrt(g*cm^3/sec^2)
cLight=2.99792458e10;       % speed of light, cm/sec
eVtoErg=1.6021766208e-12;   % 1 eV = 1.6...e-12 erg
CtoPart=2.99792458e9;       % 1 C = 1 A*sec = 2.9...e9 particles

% 
% Electron beam parameters:
%
Ekin=3.0e4;                   % kinetic energy, eV
curBeam=0.5;                  % current density, A/cm^2
dBeam=3.0;                    % beam diameter, cm
angSpread=3.0;                % angular spread, mrad
alpha=0.;                     % angle between V_i and magnetic field
trnsvT=0.5;                   % transversal temperature, eV
longT=0.0;                    % longitudinal temperature, eV (2.0e-4)
fieldB=3.e3;                  % magnetic field, Gs
nField=length(fieldB)
omega_p=1.0e9;                         % plasma frequency, 1/sec
n_e=omega_p^2*m_e/(4.*pi*q_e^2)        % plasma density, 3.1421e+08 cm-3

n_e1=8.e7;                             % plasma density, cm-3
omega_p1=sqrt(4.*pi*n_e1*q_e^2/m_e)    % plasma frequency, 5.0459e+08 1/s  
%
% Cooling system parameter:
%
coolLength=150.0;       % typical length of the coolong section, cm

% 
% Calculated parameters of the electron beam:
%
V0=sqrt(2.*Ekin*eVtoErg/m_e)                % average velocity, cm/s
rmsTrnsvVe=sqrt(2.*trnsvT*eVtoErg/m_e)      % RMS transversal velocity, cm/s
rmsLongVe=sqrt(2.*longT*eVtoErg/m_e)        % RMS longitudinal velocity, cm/s
% dens=curBeam*CtoPart/V0                   % density, 1/cm^3
% omega=sqrt(4.*pi*dens*q_e^2/m_e)          % plasma frequency, 1/s
cyclFreq=q_e*fieldB/(m_e*cLight)            % cyclotron frequency, 1/s
rmsRoLarm=rmsTrnsvVe*cyclFreq.^(-1)         % RMS Larmor radius, cm
dens=omega_p^2*m_e/(4.*pi*q_e^2)            % density, 1/cm^3
likeDebyeR=(3./dens)^(1./3.)                % "Debye" sphere with 3 electrons, cm

coolPassTime=coolLength/V0                  % time pass through cooling section, cm

powV0=round(log10(V0));
mantV0=V0/(10^powV0);
powV0=int2str(powV0);
%
% Formfactor ffForm for friction force:
%
% ffForm = 2*pi*dens*q_e^4/(m_e*V0^2)=
%        = 0.5*omega_p^2*q_e^2/V0^2;
%
% Dimension of ffForm is force:  g*cm/sec^2=erg/cm;
%
% 1 MeV/m = 1.e6*eVtoErg/100. g*cm/sec^2 = 1.e4*eVtoErg erg/cm
MeV_mToErg_cm=1.e4*eVtoErg
ffForm=-.5*omega_p^2*q_e^2/V0^2/MeV_mToErg_cm      % MeV/m
eV_mToErg_m=100.*eVtoErg
ffForm=-.5*omega_p^2*q_e^2/V0^2/eV_mToErg_m        % =-6.8226e-12 eV/m
% eV_mInErg_cm=100.*eVtoErg
ffForm=-.5*omega_p^2*q_e^2/V0^2/eVtoErg            % =-6.8226e-10 eV/cm
ffForm=100.*ffForm                                 % =-6.8226e-08 eV/m

%
% Relative velocities of electrons:
%
relVeTrnsv=rmsTrnsvVe/V0 
relVeLong=rmsLongVe/V0

%
% Scanning over relative ion velocity 
%
%=====================================================================
%
% "Low" area (rmsLongVe < Vion < rmsTrnsvVe):
%
%=====================================================================
minRelVion_l=.1*relVeTrnsv;
maxRelVion_l=relVeTrnsv;
nVion_l=50;

minRelVion_l=.4e-2*relVeTrnsv;
maxRelVion_l=relVeTrnsv;
nVion_l=100;

stepRelVion_l=(maxRelVion_l-minRelVion_l)/(nVion_l-1)

relVion_l=zeros(nVion_l,1);       % relative ion velocity (related V0) 
debyeR_l=zeros(nVion_l,1);        % Debye radius (Betacool), cm
rhoMax_l=zeros(nVion_l,1);        % maximal impact parameter, cm
rhoMin_l=zeros(nVion_l,1);        % minimal impact parameter, cm
rhoPass_l=zeros(nVion_l,1);
rhoFast_l=zeros(nVion_l,nField);

for n=1:nVion_l
    relVion_l(n)=minRelVion_l+stepRelVion_l*(n-1);
    absVion=V0*relVion_l(n);
    rhoPass_l(n)=sqrt(absVion^2+rmsLongVe^2)*coolPassTime;
    debyeR_l(n)=sqrt(absVion^2+rmsLongVe^2+rmsTrnsvVe^2)/omega_p;     
    for k=1:nField
        rhoFast_l(n,k)=sqrt(absVion^2+rmsLongVe^2)/cyclFreq(k);
    end
end

%=====================================================================
%
% "High" area (rmsLongVe < rmsTrnsvVe < Vion):
%
%=====================================================================
minRelVion_h=relVeTrnsv;
maxRelVion_h=10.*relVeTrnsv;
nVion_h=50;

stepRelVion_h=(maxRelVion_h-minRelVion_h)/(nVion_h-1)

relVion_h=zeros(nVion_h,1);       % relative ion velocity (related V0) 
debyeR_h=zeros(nVion_h,1);        % Debye radius (Betacool), cm
rhoMax_h=zeros(nVion_h,1);        % maximal impact parameter, cm
rhoMin_h=zeros(nVion_h,1);        % minimal impact parameter, cm
rhoPass_h=zeros(nVion_h,1);
rhoFast_h=zeros(nVion_h,nField);

for n=1:nVion_h
    relVion_h(n)=minRelVion_h+stepRelVion_h*(n-1);
    absVion=V0*relVion_h(n);
    rhoPass_h(n)=sqrt(absVion^2+rmsLongVe^2)*coolPassTime;
    debyeR_h(n)=sqrt(absVion^2+rmsLongVe^2+rmsTrnsvVe^2)/omega_p;     
    for k=1:nField
        rhoFast_h(n,k)=sqrt(absVion^2+rmsLongVe^2)/cyclFreq(k);
    end
end

xLimit=[.9*minRelVion_l,1.1*maxRelVion_h]

figure(205)
loglog(relVion_l,debyeR_l,'-b',relVion_h,debyeR_h,'-m','LineWidth',2)
grid on
hold on
plot ([relVion_l(1),relVion_h(nVion_h)],[likeDebyeR,likeDebyeR],'Color','k','LineWidth',2)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('R_{Debye} & R_e, cm','Color','m','FontSize',16)
title(['R_{Debye} and R_e: V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s'],'Color','m','FontSize',14)
xlim(xLimit)
ylim([1.e-3,.5])
plot([relVeTrnsv,relVeTrnsv],1.e-3*[1.,1.3],'k','LineWidth',1)
text(2.6e-3,7.e-4,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(2.e-3,3.e-3,'R_e=(3/n_e)^{1/3}','Color','k','FontSize',13)
text(2.e-3,3.5e-2,'R_{Debye}=','Color','k','FontSize',14)
plot([5.e-3,1.6e-2],[3.75e-2,3.75e-2],'Color','k')
text(5.e-3,5.e-2,'<|V_i-\Delta_{e||}|>','Color','k','FontSize',14)
text(8.e-3,3.e-2,'\omega_{ep}','Color','k','FontSize',14)

figure(207)
loglog(relVion_l,rhoPass_l,'-b',relVion_h,rhoPass_h,'-m','LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('R_{Pass}, cm','Color','m','FontSize',16)
title(['V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s, L_{Cool}=',num2str(coolLength,'%4.1f'), ...
       ' cm'],'Color','m','FontSize',14)
xlim(xLimit)
ylim([4.e-2,7.0])
plot([relVeTrnsv,relVeTrnsv],4.e-2*[1.,1.3],'k','LineWidth',1)
text(2.6e-3,3.e-2,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(2.e-3,.2,'R_{Pass}=<|V_i-\Delta_{e||}|>\cdot','Color','k','FontSize',14)
plot([1.4e-2,2.5e-2],[.205,.205],'Color','k')
text(1.4e-2,.26,'L_{Cool}','Color','k','FontSize',14)
text(1.5e-2,.16,'V_{e0}','Color','k','FontSize',14)

figure(209)
loglog(relVion_l,debyeR_l,'-b',relVion_h,debyeR_h,'-m', ...
       relVion_l,rhoPass_l,'-b',relVion_h,rhoPass_h,'-m','LineWidth',2)
grid on
hold on
plot ([relVion_l(1),relVion_h(nVion_h)],[likeDebyeR,likeDebyeR],'Color','k','LineWidth',2)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('R_{Debye} & R_{Pass} & R_e, cm','Color','m','FontSize',16)
title(['R_{Debye} & R_{Pass} & R_e: V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s'],'Color','m','FontSize',14)
xlim(xLimit)
ylim([1.e-3,10.])
plot([relVeTrnsv,relVeTrnsv],1.e-3*[1.,1.4],'k','LineWidth',1)
text(2.6e-3,6.e-4,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(3.e-3,3.25e-3,'R_e','Color','k','FontSize',14)
text(3.e-3,.11,'R_{Debye}','Color','k','FontSize',14)
text(3.e-3,1.1,'R_{Pass}','Color','k','FontSize',14)

rhoMax_l=zeros(nVion_l,1);   
for n=1:nVion_l
    help=max(debyeR_l(n),likeDebyeR);
    rhoMax_l(n)=min(rhoPass_l(n),help);
end
   
rhoMax_h=zeros(nVion_h,1);   
for n=1:nVion_h
    help=max(debyeR_h(n),likeDebyeR);
    rhoMax_h(n)=min(rhoPass_h(n),help);
end
   
figure(215)
loglog(relVion_l,rhoMax_l,'-b',relVion_h,rhoMax_h,'-m','LineWidth',2)
hold on
grid on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('R_{max}, cm','Color','m','FontSize',16)
title(['V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s'],'Color','m','FontSize',14)
xlim(xLimit)
ylim([4.e-2,.5])
plot([relVeTrnsv,relVeTrnsv],4.e-2*[1.,1.1],'k','LineWidth',1)
text(2.6e-3,3.5e-2,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(2.e-3,0.3,'R_{max}=','Color','k','FontSize',14)
text(6.e-4,0.235,'=min\{max\{R_{Debye},R_e\},R_{Pass}\}','Color','k','FontSize',14)

rhoCrit_l=zeros(nVion_l,nField);
for n=1:nVion_l
    for k=1:nField
        rhoCrit_l(n,k)=(q_e^2/m_e/cyclFreq(k)^2)^(1./3);
    end
end

rhoCrit_h=zeros(nVion_h,nField);
for n=1:nVion_h
    for k=1:nField
        rhoCrit_h(n,k)=(q_e^2/m_e/cyclFreq(k)^2)^(1./3);
    end
end

figure(305)
loglog(relVion_l,rhoFast_l(:,1),'-b',relVion_h,rhoFast_h(:,1),'-m', ...
       relVion_l,rhoCrit_l(:,1),'-b',relVion_h,rhoCrit_h(:,1),'-b','LineWidth',2)
hold on
grid on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('R_{Fast} & R_{Crit}, cm','Color','m','FontSize',16)
title(['V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s, B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','m','FontSize',14)
xlim(xLimit)
ylim([4.e-5,.01])
plot([relVeTrnsv,relVeTrnsv],4.e-5*[1.,1.3],'k','LineWidth',1)
text(2.6e-3,3.e-5,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(2.e-3,3.e-4,'R_{Fast}=<|V_i-\DeltaV_{e||}|>/\omega_{Larm}','Color','k','FontSize',13)
text(2.e-3,6.e-5,'R_{Crit}=(q_e^2/m_e/\omega_{Larm}^2)^{1/3}','Color','k','FontSize',13)

rhoMin_l=zeros(nVion_l,1);   
for n=1:nVion_l
    rhoMin_l(n)=q_e^2/m_e/((relVion_l(n)*V0)^2+rmsTrnsvVe^2+rmsLongVe^2);    
end
   
rhoMin_h=zeros(nVion_h,1);   
for n=1:nVion_h
    rhoMin_h(n)=q_e^2/m_e/((relVion_h(n)*V0)^2+rmsTrnsvVe^2+rmsLongVe^2);    
end

figure(307)
loglog(relVion_l,rhoMin_l,'-b',relVion_h,rhoMin_h,'-m','LineWidth',2)
hold on
grid on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('R_{min}, cm','Color','m','FontSize',16)
title(['V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s'],'Color','m','FontSize',14)
xlim(xLimit)
ylim([1.e-9,2.e-7])
plot([relVeTrnsv,relVeTrnsv],1.e-9*[1.,1.3],'k','LineWidth',1)
text(2.6e-3,7.5e-10,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(8.e-4,3.e-8,'R_{min}=','Color','k','FontSize',13)
plot([1.5e-3,6.e-3],[3.1e-8,3.1e-8],'Color','k')
text(2.2e-3,4.e-8,'q_e^2/m_e','Color','k','FontSize',13)
text(1.5e-3,2.4e-8,'(<|V_i-V_{e\perp}|>)^2','Color','k','FontSize',13)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Main calculations
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%=====================================================================
%
% "Low" area (rmsLongVe < Vion < rmsTrnsvVe):
%
%=====================================================================

rhoLarm_l=zeros(nVion_l,nField);

CL_l=zeros(nVion_l,nField);        % Coulomb logarithm

trnsvFF_l=zeros(nVion_l,nField);   % transverse friction force  
longFF_l=zeros(nVion_l,nField);    % longitufinal friction force          

for n=1:nVion_l
    for k=1:nField
        rhoLarm_l(n,k)=rmsTrnsvVe/cyclFreq(k);
        fctr=.5*rhoMax_l(n)/rhoLarm_l(n,k);
        if fctr > 1.
            CL_l(n,k)=log(fctr);
        end
    end
    relVion_l(n)=minRelVion_l+stepRelVion_l*(n-1);
    absVion=V0*relVion_l(n);
    trnsvVion=absVion*sin(alpha);  % transverse ion velocity
    longVion=absVion*cos(alpha);   % longitudinal ion velocity
    
    for k=1:nField
% Transverse force:        
        trnsvFF_l(n,k)=ffForm*CL_l(n,k)/relVion_l(n)^2;
% For future: 
%         trnsvFF_l(n,k)=ffForm*CL_l(n,k)*trnsvVion/absVion^3* ...
%                      (relVeTrnsv^2-2.*relVeLong^2);
% Longitudinal force:        
        longFF_l(n,k)=ffForm*2./relVion_l(n)^2;
% For future: 
%         longFF_l(n,k)=ffForm*longVion/absVion^3* ...
%                      (3*CL_l(n,k)*relVeTrnsv^2+2.);
    end
end

%=====================================================================
%
% "High" area (rmsLongVe < rmsTrnsvVe < Vion):
%
%=====================================================================

rhoLarm_h=zeros(nVion_h,nField);

CL_h=zeros(nVion_h,nField);        % Coulomb logarithm

trnsvFF_h=zeros(nVion_h,nField);   % transverse friction force  
longFF_h=zeros(nVion_h,nField);    % longitufinal friction force          

for n=1:nVion_h
    for k=1:nField
        rhoLarm_h(n,k)=rmsTrnsvVe/cyclFreq(k);
        fctr=.5*rhoMax_h(n)/rhoLarm_h(n,k);
        if fctr > 1.
            CL_h(n,k)=log(fctr);
        end
    end
    relVion_h(n)=minRelVion_h+stepRelVion_h*(n-1);
    absVion=V0*relVion_h(n);
    trnsvVion=absVion*sin(alpha);  % transverse ion velocity
    longVion=absVion*cos(alpha);   % longitudinal ion velocity
    
    for k=1:nField
% Transverse force:        
        trnsvFF_h(n,k)=ffForm*CL_h(n,k)/relVion_h(n)^2;
% For future: 
%         trnsvFF_h(n,k)=ffForm*CL_h(n,k)*trnsvVion/absVion^3* ...
%                      (relVeTrnsv^2-2.*relVeLong^2);
% Longitudinal force:        
        longFF_h(n,k)=ffForm*2./relVion_h(n)^2;
% For future: 
%         longFF_h(n,k)=ffForm*longVion/absVion^3* ...
%                      (3*CL_h(n,k)*relVeTrnsv^2+2.);
    end
end

figure(3151)
loglog(relVion_l,rhoMax_l,'-b',relVion_h,rhoMax_h,'-m', ...
       relVion_l,2.*rhoLarm_l(:,1),'-b',relVion_h,2.*rhoLarm_h(:,1),'-m','LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('R_{max} & 2\cdot<rho_\perp>, cm','Color','m','FontSize',16)
title(['V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s, B=', ...
       num2str(fieldB(1),'%6.1f'),' Gs'],'Color','m','FontSize',14)
xlim(xLimit)
ylim([1.e-3,.6])
plot([relVeTrnsv,relVeTrnsv],1.e-3*[1.,1.3],'k','LineWidth',1)
text(2.6e-3,7.e-4,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(3.e-3,1.95e-3,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(3.e-3,.1,'R_{max}','Color','k','FontSize',16)
text(5.e-4,1.e-2,'<rho_\perp>=','Color','k','FontSize',13)
plot([1.35e-3,2.35e-3],[1.025e-2,1.025e-2],'Color','k')   
text(1.35e-3,1.35e-2,'\DeltaV_{e\perp}','Color','k','FontSize',13)
text(1.35e-3,8.e-3,'\omega_{Larm}','Color','k','FontSize',13)

%----------------------------------------------
%
% Figures of Coulomb logarithms:
%
%----------------------------------------------
figure(3201)
loglog(relVion_l,CL_l(:,1),'-xb',relVion_h,CL_h(:,1),'-m', 'LineWidth',2)
hold on
grid on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{magnetized}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{magnetized}=R_{max}/(2\cdot<rho_\perp>)', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([2.e-2,5.])
plot([relVeTrnsv,relVeTrnsv],1.e-3*[1.,1.3],'k','LineWidth',1)
text(2.6e-3,7.e-4,'\DeltaV_{e\perp}/ V_{e0}','Color','k','FontSize',10)
text(1.e-5,3.5,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0,'} cm/s'], ...
       'Color','m','FontSize',14)
text(1.2e-5,4.e-1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(7.e-5,1.2e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(3.75e-4,3.e-2,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k','FontSize',14)

figure(3201)
loglog(relVion_s,CL_msb(:,1),'-r',relVion_l,CL_mlb(:,1),'-b', ...
       relVion_h,CL_mhb(:,1),'-m',relVion_s,CL_msb(:,2),'-r', ...
       relVion_l,CL_mlb(:,2),'-b',relVion_h,CL_mhb(:,2),'-m', ...
       relVion_s,CL_msb(:,3),'-r',relVion_l,CL_mlb(:,3),'-b', ...
       relVion_h,CL_mhb(:,3),'-m', 'LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{magnetized}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{magnetized}=R_{max}/(2\cdot<rho_\perp>)', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([7.e-3,5.])
text(1.05e-6,3.5,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0,'} cm/s'], ...
       'Color','m','FontSize',14)
text(2.e-6,1.5,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(2.e-6,0.3,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(1.8e-4,0.02,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k','FontSize',14)

figure(330)
loglog(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
hold on
plot([9.1e-6,2.5e-5],[.805,.805],'Color','k')   
plot([1.05e-4,9.e-4],[.805,.805],'Color','k')   
grid on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{adiabatic}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{adiabatic}=2\cdot<rho_\perp>/R_{Fast}', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([.6,6.])
text(1.e-5,5.25,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0,'} cm/s'], ...
       'Color','m','FontSize',14)
text(1.5e-6,.8,'<rho_\perp>=','Color','k','FontSize',13)
text(9.1e-6,.9,'\DeltaV_{e\perp}','Color','k','FontSize',13)
text(9.1e-6,.75,'\omega_{Larm}','Color','k','FontSize',13)
text(2.7e-5,.8,', R_{Fast}=','Color','k','FontSize',13)
text(1.75e-4,.75,'\omega_{Larm}','Color','k','FontSize',13)
text(1.05e-4,.9,'<|V_i-\DeltaV_{e||}|>','Color','k','FontSize',13)
text(1.25e-6,2.,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,1.65,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)

figure(331)
semilogx(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
hold on
grid on
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{adiabatic}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{adiabatic}=2\cdot<rho_\perp>/R_{Fast}', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([.5,5.25])
text(1.e-5,4.9,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0,'} cm/s'], ...
       'Color','m','FontSize',14)
text(1.5e-6,.9,'<rho_\perp>=','Color','k','FontSize',13)
plot([9.1e-6,2.5e-5],[.905,.905],'Color','k')   
text(9.1e-6,1.105,'\DeltaV_{e\perp}','Color','k','FontSize',13)
text(9.1e-6,.75,'\omega_{Larm}','Color','k','FontSize',13)
text(2.7e-5,.9,', R_{Fast}=','Color','k','FontSize',13)
plot([1.05e-4,9.e-4],[.905,.905],'Color','k')   
text(1.75e-4,.75,'\omega_{Larm}','Color','k','FontSize',13)
text(1.05e-4,1.105,'<|V_i-\DeltaV_{e||}|>','Color','k','FontSize',13)
text(1.25e-6,2.5,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,2.1,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)

figure(340)
semilogx(relVion_s,Numb_s,'-r',relVion_l,Numb_l,'-b','LineWidth',2)
hold on
plot([8.e-6,2.1e-5],[3.05,3.05],'Color','r')
plot([1.75e-4,3.5e-4],[3.05,3.05],'Color','b')
grid on
% xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('N_{Larm}','Color','m','FontSize',16)
title('Number of Electron Larmor Turns N_{Larm}', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([0,19])
text(1.e-5,18.,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0,'} cm/s'], ...
       'Color','m','FontSize',14)
text(1.25e-6,3.,'N_{Larm}=1+','Color','r','FontSize',13)
plot([7.5e-6,7.5e-6],[4.5,1.25],'Color','r')
text(8.e-6,3.85,'\Delta_{e\perp}','Color','r','FontSize',13)
text(8.e-6,2.15,'\pi\cdot\Delta_{e||}','Color','r','FontSize',13)
plot([2.2e-5,2.2e-5],[4.5,1.25],'Color','r')
text(2.25e-5,3.,';','Color','k','FontSize',13)
text(2.75e-5,3.,'N_{Larm}=1+','Color','b','FontSize',13)
plot([1.65e-4,1.65e-4],[4.5,1.25],'Color','b')
text(1.75e-4,3.85,'\Delta_{e\perp}','Color','b','FontSize',13)
text(1.75e-4,2.15,'\pi\cdotV_i','Color','b','FontSize',13)
plot([3.65e-4,3.65e-4],[4.5,1.25],'Color','b')

%----------------------------------------------
%
% Figures of Coulomb logarithms (continuation):
%
%----------------------------------------------
figure(350)
loglog(relVion_s,CL_fs(:,1),'-r',relVion_l,CL_fl(:,1),'-b', ...
       relVion_h,CL_fh(:,1),'-m',relVion_s,CL_fs(:,2),'-r', ...
       relVion_l,CL_fl(:,2),'-b',relVion_h,CL_fh(:,2),'-m', ...
       relVion_s,CL_fs(:,3),'-r',relVion_l,CL_fl(:,3),'-b', ...
       relVion_h,CL_fh(:,3),'-m', 'LineWidth',2)
hold on
grid on
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{Fast}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{Fast}=R_{Fast}/R_{min}', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([4.,17.])
text(1.e-5,14.5,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0,'} cm/s'], ...
       'Color','m','FontSize',14)
text(3.e-6,4.95,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(3.e-6,6.65,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(3.e-6,8.5,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(3.e-4,5.75,'R_{Fast}=','Color','k','FontSize',13)
plot([1.05e-3,1.2e-2],[5.8,5.8],'Color','k')   
text(1.05e-3,6.2,'<|V_i-\DeltaV_{e||}|>','Color','k','FontSize',13)
text(2.5e-3,5.6,'\omega_{Larm}','Color','k','FontSize',13)
text(3.e-4,4.5,'R_{min}=','Color','k','FontSize',13)
plot([1.05e-3,1.5e-2],[4.55,4.55],'Color','k')   
text(2.e-3,4.9,'q_e^2/m_e','Color','k','FontSize',13)
text(1.05e-3,4.25,'<|V_i-\DeltaV_{e\perp}|>^2','Color','k','FontSize',13)

figure(351)
semilogx(relVion_s,CL_fs(:,1),'-r',relVion_l,CL_fl(:,1),'-b', ...
         relVion_h,CL_fh(:,1),'-m',relVion_s,CL_fs(:,2),'-r', ...
         relVion_l,CL_fl(:,2),'-b',relVion_h,CL_fh(:,2),'-m', ...
         relVion_s,CL_fs(:,3),'-r',relVion_l,CL_fl(:,3),'-b', ...
         relVion_h,CL_fh(:,3),'-m', 'LineWidth',2)
hold on
plot([1.05e-3,1.2e-2],[6.75,6.75],'Color','k')   
plot([1.05e-3,1.5e-2],[5.05,5.05],'Color','k')   
grid on
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{Fast}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{Fast}=R_{Fast}/R_{min}', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([4.,17.])
text(1.e-5,15.25,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0,'} cm/s'], ...
       'Color','m','FontSize',14)
text(3.e-6,5.1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(3.e-6,6.75,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(3.e-6,8.5,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k','FontSize',14)
text(3.e-4,6.7,'R_{Fast}=','Color','k','FontSize',13)
text(1.05e-3,7.25,'<|V_i-\DeltaV_{e||}|>','Color','k','FontSize',13)
text(2.5e-3,6.5,'\omega_{Larm}','Color','k','FontSize',13)
text(3.e-4,5.0,'R_{min}=','Color','k','FontSize',13)
text(2.e-3,5.6,'q_e^2/m_e','Color','k','FontSize',13)
text(1.05e-3,4.5,'<|V_i-\DeltaV_{e\perp}|>^2','Color','k','FontSize',13)
text(1.25e-6,13.5,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,12.5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-6,11.5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
 
%
% B=100 Gs ("Meshkov" approach):
%
figure(360)
loglog(relVion_s,CL_msm(:,3),'-r',relVion_l,CL_mlm(:,3),'-b', ...
       relVion_h,CL_mhm(:,3),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
loglog(relVion_s,CL_fs(:,3),'-r',relVion_l,CL_fl(:,3),'-b',  ...
       relVion_h,CL_fh(:,3),'-m','LineWidth',2)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Coulomb Logarithm','Color','m','FontSize',16)
title('Coulomb Logarithm: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([2.e-2,20.])
text(1.2e-6,3.e-2,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(3),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(6.e-5,12.,'Fast','Color','k','FontSize',16)
text(4.e-4,4.,'Adiabatic','Color','k','FontSize',16)
text(3.25e-4,.15,'Magnetized','Color','k','FontSize',16)
text(1.25e-6,1.35,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,.8,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-6,.5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=3000 Gs ("Meshkov" approach):
%
figure(361)
loglog(relVion_s,CL_msm(:,1),'-r',relVion_l,CL_mlm(:,1),'-b', ...
       relVion_h,CL_mhm(:,1),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
loglog(relVion_s,CL_fs(:,1),'-r',relVion_l,CL_fl(:,1),'-b',  ...
       relVion_h,CL_fh(:,1),'-m','LineWidth',2)
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Coulomb Logarithm','Color','m','FontSize',16)
title('Coulomb Logarithm: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([2.e-1,13.])
text(1.05e-6,10.,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.5e-4,7.,'Fast','Color','k','FontSize',16)
text(2.25e-4,4.,'Adiabatic','Color','k','FontSize',16)
text(3.e-4,.5,'Magnetized','Color','k','FontSize',16)
text(1.25e-6,2.,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,1.45,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-6,1.05,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% For comparison with figure 4 from [3] only: B=1000 Gs ("Meshkov" approach):
%
if (reference3flag ==1)
    figure(8361)
    semilogx(relVion_s*V0*.01,CL_msm(:,4),'-r',relVion_l*V0*.01,CL_mlm(:,4),'-b', ...
           relVion_h*V0*.01,CL_mhm(:,4),'-m','LineWidth',2)
    grid on
    hold on
    semilogx(relVion_s*V0*.01,CL_as,'-r',relVion_l*V0*.01,CL_al,'-b','LineWidth',2)
    semilogx(relVion_s*V0*.01,CL_fs(:,4),'-r',relVion_l*V0*.01,CL_fl(:,4),'-b',  ...
           relVion_h*V0*.01,CL_fh(:,4),'-m','LineWidth',2)
%     xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
    xlabel('Ion Velocity, V_i, m/s','Color','m','FontSize',16)
    ylabel('Coulomb Logarithm','Color','m','FontSize',16)
    title('Coulomb Logarithm: 3 types of Interaction', ...
          'Color','m','FontSize',13)
%     xlim([2.e-7*V0*.01,3.e-3*V0*.01])
    xlim(V0*.01*xLimit)
    ylim([0.,13.])
    text(1.05e-6,10.,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
         '} cm/s,  B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
           'Color','k','FontSize',16)
    text(1.5e-4,7.,'Fast','Color','k','FontSize',16)
    text(2.25e-4,4.,'Adiabatic','Color','k','FontSize',16)
    text(3.e-4,.5,'Magnetized','Color','k','FontSize',16)
    text(1.25e-6,2.,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
         'Color','r','FontSize',13)
    text(1.25e-6,1.45,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
         'Color','b','FontSize',13)
    text(1.25e-6,1.05,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
         'Color','m','FontSize',13)
end

%
% B=600 Gs ("Meshkov" approach):
%
figure(362)
loglog(relVion_s,CL_msm(:,2),'-r',relVion_l,CL_mlm(:,2),'-b', ...
       relVion_h,CL_mhm(:,2),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
loglog(relVion_s,CL_fs(:,2),'-r',relVion_l,CL_fl(:,2),'-b',  ...
       relVion_h,CL_fh(:,2),'-m','LineWidth',2)
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Coulomb Logarithm','Color','m','FontSize',16)
title('Coulomb Logarithm: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([5.e-2,17.])
text(1.05e-6,12.,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(2),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(4.e-4,6.,'Fast','Color','k','FontSize',16)
text(9.e-5,2.,'Adiabatic','Color','k','FontSize',16)
text(5.5e-5,.07,'Magnetized','Color','k','FontSize',16)
text(1.25e-6,.7,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,.45,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-6,.3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=100 Gs ("Betacool" approach):
%
figure(3601)
loglog(relVion_s,CL_msb(:,3),'-r',relVion_l,CL_mlb(:,3),'-b', ...
       relVion_h,CL_mhb(:,3),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
loglog(relVion_s,CL_fs(:,3),'-r',relVion_l,CL_fl(:,3),'-b',  ...
       relVion_h,CL_fh(:,3),'-m','LineWidth',2)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Coulomb Logarithm','Color','m','FontSize',16)
title('Coulomb Logarithm: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([6.e-3,20.])
text(1.2e-6,1.e-2,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(3),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(6.e-5,12.,'Fast','Color','k','FontSize',16)
text(4.e-4,4.,'Adiabatic','Color','k','FontSize',16)
text(1.5e-4,.03,'Magnetized','Color','k','FontSize',16)
text(1.25e-6,1.35,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,.8,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-6,.45,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=3000 Gs ("Betacool" approach):
%
figure(3611)
loglog(relVion_s,CL_msb(:,1),'-r',relVion_l,CL_mlb(:,1),'-b', ...
       relVion_h,CL_mhb(:,1),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
loglog(relVion_s,CL_fs(:,1),'-r',relVion_l,CL_fl(:,1),'-b',  ...
       relVion_h,CL_fh(:,1),'-m','LineWidth',2)
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Coulomb Logarithm','Color','m','FontSize',16)
title('Coulomb Logarithm: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([5.5e-1,13.])
text(1.05e-6,10.,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.5e-4,7.,'Fast','Color','k','FontSize',16)
text(1.5e-3,2.,'Adiabatic','Color','k','FontSize',16)
text(2.e-6,2.3,'Magnetized','Color','k','FontSize',16)
text(1.25e-6,1.5,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,1.2,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-6,.95,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=600 Gs ("Betacool" approach):
%
figure(3621)
loglog(relVion_s,CL_msb(:,2),'-r',relVion_l,CL_mlb(:,2),'-b', ...
       relVion_h,CL_mhb(:,2),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,CL_as,'-r',relVion_l,CL_al,'-b','LineWidth',2)
loglog(relVion_s,CL_fs(:,2),'-r',relVion_l,CL_fl(:,2),'-b',  ...
       relVion_h,CL_fh(:,2),'-m','LineWidth',2)
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Coulomb Logarithm','Color','m','FontSize',16)
title('Coulomb Logarithm: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([.4,17.])
text(1.05e-6,13.,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(2),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(4.e-4,7.,'Fast','Color','k','FontSize',16)
text(2.e-4,4.,'Adiabatic','Color','k','FontSize',16)
text(6.e-5,.6,'Magnetized','Color','k','FontSize',16)
text(1.25e-6,2.5,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-6,1.9,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-6,1.4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%----------------------------------------------
%
% Figures of impact parameters:
%
%----------------------------------------------
%
% B=3000 Gs ("Meshkov" approach):
%
figure(370)
loglog(relVion_s,rhoMax_sm,'-r',relVion_l,rhoMax_lm,'-b', ...
       relVion_h,rhoMax_hm,'-m','LineWidth',2)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,1),'-r',relVion_l,2.*rhoLarm_l(:,1),'-b', ...
       relVion_h,2.*rhoLarm_h(:,1),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,1),'-r',relVion_l,rhoFast_l(:,1),'-b', ...
       relVion_h,rhoFast_h(:,1),'-m','LineWidth',2)
loglog(relVion_s,rhoMin_sm,'-r',relVion_l,rhoMin_lm,'-b', ...
       relVion_h,rhoMin_hm,'-m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,1),'--r',relVion_l,rhoCrit_l(:,1),'--b', ...
       relVion_h,rhoCrit_h(:,1),'--m','LineWidth',2)
grid on
xlim(xLimit)
ylim([7.e-9,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,8.e-2,['B=',num2str(fieldB(1),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,6.e-3,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,8.e-4,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,7.e-6,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,1.1e-4,'R_{Crit}','Color','k','FontSize',14)

%
% Figure (370) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%
nVion_h_max=nVion_h;
for i=1:nVion_h
    if (2.*rhoLarm_h(i,1) < rhoFast_h(i,1)) && (nVion_h_max == nVion_h)
        nVion_h_max=i;
    end
end

% Figure (370) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%
figure(371)
loglog(relVion_s,rhoMax_sm,'-k',relVion_l,rhoMax_lm,'-k', ...
       relVion_h,rhoMax_hm,'-k','LineWidth',1)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,1),'-k',relVion_l,2.*rhoLarm_l(:,1),'-k', ...
       relVion_h,2.*rhoLarm_h(:,1),'-k','LineWidth',1)
loglog(relVion_s,rhoFast_s(:,1),'-k',relVion_l,rhoFast_l(:,1),'-k', ...
       relVion_h(1:nVion_h_max),rhoFast_h(1:nVion_h_max,1),'-k','LineWidth',1)
loglog(relVion_s,rhoMin_sm,'-k',relVion_l,rhoMin_lm,'-k', ...
       relVion_h,rhoMin_hm,'-k','LineWidth',1)
loglog(relVion_s,rhoCrit_s(:,1),'--k',relVion_l,rhoCrit_l(:,1),'--k', ...
       relVion_h,rhoCrit_h(:,1),'--k','LineWidth',1)
grid on
xlim(xLimit)
ylim([7.e-9,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,8.e-2,['B=',num2str(fieldB(1),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,6.e-3,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,8.e-4,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,7.e-6,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,1.1e-4,'R_{Crit}','Color','k','FontSize',14)
text(1.5e-5,1.e-2,'Screened','Color','r','FontSize',18)
text(6.e-4,3.5e-3,'Magnetized','Color','r','FontSize',18)
text(3.e-5,4.e-4,'Adiabatic','Color','r','FontSize',18)
text(1.e-4,2.e-6,'Fast','Color','r','FontSize',18)
text(4.e-5,4.e-8,'Weak','Color','r','FontSize',18)

%----------------------------------------------
%
% B=3000 Gs ("Betacool" approach):
%
figure(3701)
loglog(relVion_s,rhoMax_sb,'-r',relVion_l,rhoMax_lb,'-b', ...
       relVion_h,rhoMax_hb,'-m','LineWidth',2)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,1),'-r',relVion_l,2.*rhoLarm_l(:,1),'-b', ...
       relVion_h,2.*rhoLarm_h(:,1),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,1),'-r',relVion_l,rhoFast_l(:,1),'-b', ...
       relVion_h,rhoFast_h(:,1),'-m','LineWidth',2)
loglog(relVion_s,rhoMin_sb,'-r',relVion_l,rhoMin_lb,'-b', ...
       relVion_h,rhoMin_hb,'-m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,1),'--r',relVion_l,rhoCrit_l(:,1),'--b', ...
       relVion_h,rhoCrit_h(:,1),'--m','LineWidth',2)
grid on
xlim(xLimit)
ylim([7.e-9,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,10.e-2,['B=',num2str(fieldB(1),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,3.e-2,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,8.e-4,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,7.e-6,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,1.1e-4,'R_{Crit}','Color','k','FontSize',14)

%
% Figure (3701) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%
nVion_h_max=nVion_h;
for i=1:nVion_h
    if (2.*rhoLarm_h(i,1) < rhoFast_h(i,1)) && (nVion_h_max == nVion_h)
        nVion_h_max=i;
    end
end

% Figure (3701) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%
figure(3711)
loglog(relVion_s,rhoMax_sb,'-k',relVion_l,rhoMax_lb,'-k', ...
       relVion_h,rhoMax_hb,'-k','LineWidth',1)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,1),'-k',relVion_l,2.*rhoLarm_l(:,1),'-k', ...
       relVion_h,2.*rhoLarm_h(:,1),'-k','LineWidth',1)
loglog(relVion_s,rhoFast_s(:,1),'-k',relVion_l,rhoFast_l(:,1),'-k', ...
       relVion_h(1:nVion_h_max),rhoFast_h(1:nVion_h_max,1),'-k','LineWidth',1)
loglog(relVion_s,rhoMin_sb,'-k',relVion_l,rhoMin_lb,'-k', ...
       relVion_h,rhoMin_hb,'-k','LineWidth',1)
loglog(relVion_s,rhoCrit_s(:,1),'--k',relVion_l,rhoCrit_l(:,1),'--k', ...
       relVion_h,rhoCrit_h(:,1),'--k','LineWidth',1)
grid on
xlim(xLimit)
ylim([7.e-9,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,10.e-2,['B=',num2str(fieldB(1),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,3.e-2,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,8.e-4,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,7.e-6,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,1.1e-4,'R_{Crit}','Color','k','FontSize',14)
text(7.e-6,3.e-2,'Screened','Color','r','FontSize',18)
text(3.e-5,5.e-3,'Magnetized','Color','r','FontSize',18)
text(3.e-5,4.e-4,'Adiabatic','Color','r','FontSize',18)
text(3.e-5,2.e-6,'Fast','Color','r','FontSize',18)
text(3.e-5,4.e-8,'Weak','Color','r','FontSize',18)

%----------------------------------------------
%
% B=600 Gs ("Meshkov" approach):
%
figure(373)
loglog(relVion_s,rhoMax_sm,'-r',relVion_l,rhoMax_lm,'-b', ...
       relVion_h,rhoMax_hm,'-m','LineWidth',2)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,2),'-r',relVion_l,2.*rhoLarm_l(:,2),'-b', ...
       relVion_h,2.*rhoLarm_h(:,2),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,2),'-r',relVion_l,rhoFast_l(:,2),'-b', ...
       relVion_h,rhoFast_h(:,2),'-m','LineWidth',2)
loglog(relVion_s,rhoMin_sm,'-r',relVion_l,rhoMin_lm,'-b', ...
       relVion_h,rhoMin_hm,'-m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,2),'--r',relVion_l,rhoCrit_l(:,2),'--b', ...
       relVion_h,rhoCrit_h(:,2),'--m','LineWidth',2)
grid on
xlim(xLimit)
ylim([6.e-9,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,.1,['B=',num2str(fieldB(2),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,1.e-3,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,1.5e-2,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,4.e-5,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,3.5e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,3.e-4,'R_{Crit}','Color','k','FontSize',14)
%
% Figure (374) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%
nVion_l_min=1;
for i=1:nVion_l
    if (2.*rhoLarm_l(i,2) < rhoMax_lm(i))
        nVion_l_min=i;
        break;
    end
end

nVion_h_max=1;
for i=1:nVion_h
    if (2.*rhoLarm_h(i,2) < rhoFast_h(i,2))
        nVion_h_max=i;
        break;
    end
end

figure(374)
loglog(relVion_s,rhoMax_sm,'-k',relVion_l,rhoMax_lm,'-k', ...
       relVion_h,rhoMax_hm,'-k','LineWidth',1)
hold on
loglog(relVion_l(nVion_l_min:nVion_l),2.*rhoLarm_l(nVion_l_min:nVion_l,2),'-k', ...
       relVion_h(1:nVion_h_max),2.*rhoLarm_h(1:nVion_h_max,2),'-k', 'LineWidth',1)
loglog(relVion_s,rhoFast_s(:,2),'-k',relVion_l,rhoFast_l(:,2),'-k', ...
       relVion_h,rhoFast_h(:,2),'-k','LineWidth',1)
loglog(relVion_s,rhoMin_sm,'-k',relVion_l,rhoMin_lm,'-k', ...
       relVion_h,rhoMin_hm,'-k','LineWidth',1)
loglog(relVion_s,rhoCrit_s(:,2),'--k',relVion_l,rhoCrit_l(:,2),'--k', ...
       relVion_h,rhoCrit_h(:,2),'--k','LineWidth',1)
grid on
xlim(xLimit)
ylim([1.e-8,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(2.e-6,.1,['B=',num2str(fieldB(2),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,5.e-3,'R_{max}','Color','k','FontSize',14)
text(1.5e-4,2.e-2,'2\cdot<rho_\perp>','Color','k','FontSize',14)
arrw=annotation('arrow',[.7 .73],[.8 .76]);
text(1.5e-6,4.e-5,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,3.5e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,3.e-4,'R_{Crit}','Color','k','FontSize',14)
text(1.e-5,6.e-3,'Screened','Color','r','FontSize',18)
text(1.e-5,5.e-4,'Adiabatic','Color','r','FontSize',18)
text(1.e-5,4.e-6,'Fast','Color','r','FontSize',18)
text(1.e-5,5.e-8,'Weak','Color','r','FontSize',18)
text(2.5e-4,.1,'Magnetized','Color','r','FontSize',18)
arrw=annotation('arrow',[.73 .85],[.845 .8]);
arrw.Color=[1 0 0];

%----------------------------------------------
%
% B=600 Gs ("Betacool" approach):
%
figure(3731)
loglog(relVion_s,rhoMax_sb,'-r',relVion_l,rhoMax_lb,'-b', ...
       relVion_h,rhoMax_hb,'-m','LineWidth',2)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,2),'-r',relVion_l,2.*rhoLarm_l(:,2),'-b', ...
       relVion_h,2.*rhoLarm_h(:,2),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,2),'-r',relVion_l,rhoFast_l(:,2),'-b', ...
       relVion_h,rhoFast_h(:,2),'-m','LineWidth',2)
loglog(relVion_s,rhoMin_sb,'-r',relVion_l,rhoMin_lb,'-b', ...
       relVion_h,rhoMin_hb,'-m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,2),'--r',relVion_l,rhoCrit_l(:,2),'--b', ...
       relVion_h,rhoCrit_h(:,2),'--m','LineWidth',2)
grid on
xlim(xLimit)
ylim([1.e-8,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,.1,['B=',num2str(fieldB(2),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,3.e-2,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,4.e-3,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,4.e-5,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,3.5e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,3.e-4,'R_{Crit}','Color','k','FontSize',14)
%
% Figure (3741) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%

nVion_h_max=1;
for i=1:nVion_h
    if (2.*rhoLarm_h(i,2) < rhoFast_h(i,2))
        nVion_h_max=i;
        break;
    end
end

figure(3741)
loglog(relVion_s,rhoMax_sb,'-k',relVion_l,rhoMax_lb,'-k', ...
       relVion_h,rhoMax_hb,'-k','LineWidth',1)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,2),'-k',relVion_l,2.*rhoLarm_l(:,2),'-k', ...
       relVion_h(1:nVion_h_max),2.*rhoLarm_h(1:nVion_h_max,2),'-k', 'LineWidth',1)
loglog(relVion_s,rhoFast_s(:,2),'-k',relVion_l,rhoFast_l(:,2),'-k', ...
       relVion_h,rhoFast_h(:,2),'-k','LineWidth',1)
loglog(relVion_s,rhoMin_sb,'-k',relVion_l,rhoMin_lb,'-k', ...
       relVion_h,rhoMin_hb,'-k','LineWidth',1)
loglog(relVion_s,rhoCrit_s(:,2),'--k',relVion_l,rhoCrit_l(:,2),'--k', ...
       relVion_h,rhoCrit_h(:,2),'--k','LineWidth',1)
grid on
xlim(xLimit)
ylim([1.e-8,2.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(2.e-6,.1,['B=',num2str(fieldB(2),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,3.e-2,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,4.e-3,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,4.e-5,'R_{Fast}','Color','k','FontSize',14)
text(1.5e-6,3.5e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,3.e-4,'R_{Crit}','Color','k','FontSize',14)
text(2.5e-4,.1,'Screened','Color','r','FontSize',18)
text(1.e-5,5.e-4,'Adiabatic','Color','r','FontSize',18)
text(1.e-5,4.e-6,'Fast','Color','r','FontSize',18)
text(1.e-5,5.e-8,'Weak','Color','r','FontSize',18)
text(2.5e-4,2.e-2,'Magnetized','Color','r','FontSize',18)

%----------------------------------------------
%
% B=100 Gs ("Meshkov" approach):
%
figure(377)
loglog(relVion_s,rhoMax_sm,'-r',relVion_l,rhoMax_lm,'-b', ...
       relVion_h,rhoMax_hm,'-m','LineWidth',2)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,3),'-r',relVion_l,2.*rhoLarm_l(:,3),'-b', ...
       relVion_h,2.*rhoLarm_h(:,3),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,3),'-r',relVion_l,rhoFast_l(:,3),'-b', ...
       relVion_h,rhoFast_h(:,3),'-m','LineWidth',2)
loglog(relVion_s,rhoMin_sm,'-r',relVion_l,rhoMin_lm,'-b', ...
       relVion_h,rhoMin_hm,'-m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,3),'--r',relVion_l,rhoCrit_l(:,3),'--b', ...
       relVion_h,rhoCrit_h(:,3),'--m','LineWidth',2)
grid on
xlim(xLimit)
ylim([7.e-9,3.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,.15,['B=',num2str(fieldB(3),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,6.e-3,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,1.2e-1,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,2.e-4,'R_{Fast}','Color','k','FontSize',13)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,1.e-3,'R_{Crit}','Color','k','FontSize',14)

%
% Figure (378) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%
nVion_h_min=1;
for i=1:nVion_h
    if (2.*rhoLarm_h(i,3) < rhoMax_hm(i))
        nVion_h_min=i;
        break;
    end
end

nVion_h_max=nVion_h;
for i=1:nVion_h
    if (2.*rhoLarm_h(i,3) < rhoFast_h(i,3))
        nVion_h_max=i;
        break;
    end
end

figure(378)
loglog(relVion_s,rhoMax_sm,'-k',relVion_l,rhoMax_lm,'-k', ...
       relVion_h,rhoMax_hm,'-k','LineWidth',1)
hold on
loglog(relVion_h(nVion_h_min:nVion_h), ...
       2.*rhoLarm_h(nVion_h_min:nVion_h,3),'-k','LineWidth',1)
loglog(relVion_s,rhoFast_s(:,3),'-k',relVion_l,rhoFast_l(:,3),'-k', ...
       relVion_h(1:nVion_h_max),rhoFast_h(1:nVion_h_max,3),'-k','LineWidth',1)
loglog(relVion_s,rhoMin_sm,'-k',relVion_l,rhoMin_lm,'-k', ...
       relVion_h,rhoMin_hm,'-k','LineWidth',1)
loglog(relVion_s,rhoCrit_s(:,3),'--k',relVion_l,rhoCrit_l(:,3),'--k', ...
       relVion_h,rhoCrit_h(:,3),'--k','LineWidth',1)
grid on
xlim(xLimit)
ylim([7.e-9,3.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(2.e-6,.15,['B=',num2str(fieldB(3),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,6.e-3,'R_{max}','Color','k','FontSize',14)
text(2.e-3,6.e-3,'2\cdot<rho_\perp>','Color','k','FontSize',14)
arrw=annotation('arrow',[.84 .83],[.75 .82]);
text(1.5e-6,2.e-4,'R_{Fast}','Color','k','FontSize',13)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,2.e-4,'R_{Crit}','Color','k','FontSize',14)
text(1.e-5,2.e-2,'Screened','Color','r','FontSize',18)
text(3.e-6,1.e-3,'Adiabatic','Color','r','FontSize',18)
text(4.e-5,2.e-5,'Fast','Color','r','FontSize',18)
text(4.e-5,4.e-8,'Weak','Color','r','FontSize',18)
text(3.e-4,.15,'Magnetized','Color','r','FontSize',18)
arrw=annotation('arrow',[.84 .875],[.86 .835]);
arrw.Color=[1 0 0];

%----------------------------------------------
%
% B=100 Gs and 3000 Gs ("Meshkov" approach):
%
figure(380)
loglog(relVion_s,rhoMax_sm,'-.r',relVion_l,rhoMax_lm,'-.b', ...
       relVion_h,rhoMax_hm,'-.m','LineWidth',3)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,3),'-r',relVion_l,2.*rhoLarm_l(:,3),'-b', ...
       relVion_h,2.*rhoLarm_h(:,3),'-m','LineWidth',2)
loglog(relVion_s,2.*rhoLarm_s(:,1),'-r',relVion_l,2.*rhoLarm_l(:,1),'-b', ...
       relVion_h,2.*rhoLarm_h(:,1),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,3),'--r',relVion_l,rhoFast_l(:,3),'--b', ...
       relVion_h,rhoFast_h(:,3),'--m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,1),'--r',relVion_l,rhoFast_l(:,1),'--b', ...
       relVion_h,rhoFast_h(:,1),'--m','LineWidth',2)
loglog(relVion_s,rhoMin_sm,'-.r',relVion_l,rhoMin_lm,'-.b', ...
       relVion_h,rhoMin_hm,'-.m','LineWidth',3)
loglog(relVion_s,rhoCrit_s(:,3),'-.r',relVion_l,rhoCrit_l(:,3),'-.b', ...
       relVion_h,rhoCrit_h(:,3),'-.m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,1),'-.r',relVion_l,rhoCrit_l(:,1),'-.b', ...
       relVion_h,rhoCrit_h(:,1),'-.m','LineWidth',2)
grid on
xlim(xLimit)
ylim([7.e-9,3.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Impact Parameter (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-4,2.e-2,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,.09,'2\cdot<rho_\perp>_1','Color','k','FontSize',12)
text(1.5e-6,4.e-3,'2\cdot<rho_\perp>_2','Color','k','FontSize',12)
annotation('arrow',[.25 .3],[.705 .675]);
text(3.e-4,3.5e-7,'R_{min}','Color','k','FontSize',14)
text(2.e-3,5.e-3,'R_{Fast_1}','Color','k','FontSize',14)
% text(3.e-4,1.5e-4,'R_{Fast_2}','Color','k','FontSize',14)
text(1.1e-4,1.5e-4,'R_{Fast_2}','Color','k','FontSize',14)
text(5.e-3,2.e-4,'R_{Crit_1}','Color','k','FontSize',14)
text(5.e-3,1.7e-5,'R_{Crit_2}','Color','k','FontSize',14)
text(3.e-6,6.e-6,['Index "1": B=',num2str(fieldB(3),'%6.1f'),'  Gs'], ...
     'Color','k','FontSize',14)
text(3.e-6,2.e-6,['Index "2": B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
     'Color','k','FontSize',14)

%----------------------------------------------
%
% B=100 Gs ("Betacool" approach):
%
figure(3771)
loglog(relVion_s,rhoMax_sb,'-r',relVion_l,rhoMax_lb,'-b', ...
       relVion_h,rhoMax_hb,'-m','LineWidth',2)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,3),'-r',relVion_l,2.*rhoLarm_l(:,3),'-b', ...
       relVion_h,2.*rhoLarm_h(:,3),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,3),'-r',relVion_l,rhoFast_l(:,3),'-b', ...
       relVion_h,rhoFast_h(:,3),'-m','LineWidth',2)
loglog(relVion_s,rhoMin_sb,'-r',relVion_l,rhoMin_lb,'-b', ...
       relVion_h,rhoMin_hb,'-m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,3),'--r',relVion_l,rhoCrit_l(:,3),'--b', ...
       relVion_h,rhoCrit_h(:,3),'--m','LineWidth',2)
grid on
xlim(xLimit)
ylim([1.e-8,3.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-5,.15,['B=',num2str(fieldB(3),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,6.e-3,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,1.2e-1,'2\cdot<rho_\perp>','Color','k','FontSize',14)
text(1.5e-6,2.e-4,'R_{Fast}','Color','k','FontSize',13)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,1.e-3,'R_{Crit}','Color','k','FontSize',14)

%
% Figure (3781) for drawing of 3 areas (fast, adiabatic and magnetized
% collisions):
%
nVion_l_min=1;
for i=1:nVion_l
    if (2.*rhoLarm_l(i,3) < rhoMax_lb(i))
        nVion_l_min=i;
        break;
    end
end

figure(3781)
loglog(relVion_s,rhoMax_sb,'-k',relVion_l,rhoMax_lb,'-k', ...
       relVion_h,rhoMax_hb,'-k','LineWidth',1)
hold on
loglog(relVion_l(nVion_l_min:nVion_l),2.*rhoLarm_l(nVion_l_min:nVion_l,3), ...
       relVion_h,2.*rhoLarm_h(:,3),'-k','LineWidth',1)
loglog(relVion_s,rhoFast_s(:,3),'-k',relVion_l,rhoFast_l(:,3),'-k', ...
       relVion_h(1:nVion_h_max),rhoFast_h(1:nVion_h_max,3),'-k','LineWidth',1)
loglog(relVion_s,rhoMin_sb,'-k',relVion_l,rhoMin_lb,'-k', ...
       relVion_h,rhoMin_hb,'-k','LineWidth',1)
loglog(relVion_s,rhoCrit_s(:,3),'--k',relVion_l,rhoCrit_l(:,3),'--k', ...
       relVion_h,rhoCrit_h(:,3),'--k','LineWidth',1)
grid on
xlim(xLimit)
ylim([7.e-9,3.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Cooling Types Interaction (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(2.e-6,.15,['B=',num2str(fieldB(3),'%6.2f'),' Gs'], ...
                  'Color','k','FontSize',14)
text(1.5e-6,6.e-3,'R_{max}','Color','k','FontSize',14)
text(2.e-3,6.e-3,'2\cdot<rho_\perp>','Color','k','FontSize',14)
arrw=annotation('arrow',[.84 .83],[.75 .82]);
text(1.5e-6,2.e-4,'R_{Fast}','Color','k','FontSize',13)
text(1.5e-6,4.e-7,'R_{min}','Color','k','FontSize',14)
text(5.e-3,2.e-4,'R_{Crit}','Color','k','FontSize',14)
text(3.e-6,4.e-2,'Screened','Color','r','FontSize',18)
text(1.e-5,3.e-3,'Adiabatic','Color','r','FontSize',18)
text(4.e-5,2.e-5,'Fast','Color','r','FontSize',18)
text(4.e-5,4.e-8,'Weak','Color','r','FontSize',18)
text(3.e-4,.15,'Magnetized','Color','r','FontSize',18)
arrw=annotation('arrow',[.84 .875],[.86 .835]);
arrw.Color=[1 0 0];

%----------------------------------------------
%
% B=100 Gs and 3000 Gs ("Betacool" approach):
%
figure(3801)
loglog(relVion_s,rhoMax_sb,'-.r',relVion_l,rhoMax_lb,'-.b', ...
       relVion_h,rhoMax_hb,'-.m','LineWidth',3)
hold on
loglog(relVion_s,2.*rhoLarm_s(:,3),'-r',relVion_l,2.*rhoLarm_l(:,3),'-b', ...
       relVion_h,2.*rhoLarm_h(:,3),'-m','LineWidth',2)
loglog(relVion_s,2.*rhoLarm_s(:,1),'-r',relVion_l,2.*rhoLarm_l(:,1),'-b', ...
       relVion_h,2.*rhoLarm_h(:,1),'-m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,3),'--r',relVion_l,rhoFast_l(:,3),'--b', ...
       relVion_h,rhoFast_h(:,3),'--m','LineWidth',2)
loglog(relVion_s,rhoFast_s(:,1),'--r',relVion_l,rhoFast_l(:,1),'--b', ...
       relVion_h,rhoFast_h(:,1),'--m','LineWidth',2)
loglog(relVion_s,rhoMin_sb,'-.r',relVion_l,rhoMin_lb,'-.b', ...
       relVion_h,rhoMin_hb,'-.m','LineWidth',3)
loglog(relVion_s,rhoCrit_s(:,3),'-.r',relVion_l,rhoCrit_l(:,3),'-.b', ...
       relVion_h,rhoCrit_h(:,3),'-.m','LineWidth',2)
loglog(relVion_s,rhoCrit_s(:,1),'-.r',relVion_l,rhoCrit_l(:,1),'-.b', ...
       relVion_h,rhoCrit_h(:,1),'-.m','LineWidth',2)
grid on
xlim(xLimit)
ylim([1.e-8,3.e-1])
xlabel('Relative Ion Velocity, V_i/V_{e0}','Color','m','FontSize',16)
ylabel('Impact Parameter, cm','Color','m','FontSize',16)
title(['Impact Parameter (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',14)
text(3.e-4,2.e-2,'R_{max}','Color','k','FontSize',14)
text(1.5e-6,.09,'2\cdot<rho_\perp>_1','Color','k','FontSize',12)
text(1.5e-6,4.e-3,'2\cdot<rho_\perp>_2','Color','k','FontSize',12)
annotation('arrow',[.25 .3],[.705 .675]);
text(3.e-4,3.5e-7,'R_{min}','Color','k','FontSize',14)
text(2.e-3,5.e-3,'R_{Fast_1}','Color','k','FontSize',14)
% text(3.e-4,1.5e-4,'R_{Fast_2}','Color','k','FontSize',14)
text(1.1e-4,1.5e-4,'R_{Fast_2}','Color','k','FontSize',14)
text(5.e-3,2.e-4,'R_{Crit_1}','Color','k','FontSize',14)
text(5.e-3,1.7e-5,'R_{Crit_2}','Color','k','FontSize',14)
text(3.e-6,6.e-6,['Index "1": B=',num2str(fieldB(3),'%6.1f'),'  Gs'], ...
     'Color','k','FontSize',14)
text(3.e-6,2.e-6,['Index "2": B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
     'Color','k','FontSize',14)

%--------------------------------------------------------------
%
%    F  i g u r e s   f o r   T r a n s v e r s e   F r i c t i o n    F o r c e
%
%--------------------------------------------------------------

%
% B=3000 Gs (Three types of interaction; "Meshkov" approach):
%       
figure(390)
loglog(relVion_s,abs(trnsvFF_fs(:,1)),'-r',relVion_l,abs(trnsvFF_fl(:,1)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,1)),'--r',relVion_l,abs(trnsvFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msm(:,1)),'-.r',relVion_l,abs(trnsvFF_mlm(:,1)),'-.b', ...
       relVion_h,abs(trnsvFF_mhm(:,1)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([1.e-5,8.e-0])
text(2.e-6,4.e-0,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(2.e-4,1.5e-3,'Fast','Color','k','FontSize',16)
text(3.e-5,2.e-2,'Adiabatic','Color','k','FontSize',16)
text(1.65e-5,4.e-1,'Magnetized','Color','k','FontSize',16)
text(1.25e-5,1.5e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-5,5.e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-5,2.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% For comparison with figure 4 from [3] only: B=1000 Gs ("Meshkov" approach):
%
if (reference3flag ==1)
figure(8390)
loglog(relVion_s,abs(trnsvFF_fs(:,4)),'-r',relVion_l,abs(trnsvFF_fl(:,4)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,4)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,4)),'--r',relVion_l,abs(trnsvFF_al(:,4)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msm(:,4)),'-.r',relVion_l,abs(trnsvFF_mlm(:,4)),'-.b', ...
       relVion_h,abs(trnsvFF_mhm(:,4)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(V0*.01*xLimit)
ylim([2.e-5,3.e+2])
text(4.e-7,1.e+2,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(4),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-5,8.e-3,'Fast','Color','k','FontSize',16)
text(1.e-5,1.e-1,'Adiabatic','Color','k','FontSize',16)
text(2.e-6,1.e+1,'Magnetized','Color','k','FontSize',16)
text(2.e-6,1.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(2.e-6,2.25e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(2.e-6,7.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
end

%
% Magnetized interaction (three values of magnetic field; "Meshkov" approach):
%       
figure(395)
loglog(relVion_s,abs(trnsvFF_msm(:,1)),'-r',relVion_l,abs(trnsvFF_mlm(:,1)),'-b', ...
       relVion_h,abs(trnsvFF_mhm(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_msm(:,2)),'--r',relVion_l,abs(trnsvFF_mlm(:,2)),'--b', ...
       relVion_h,abs(trnsvFF_mhm(:,2)),'--m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msm(:,3)),'-.r',relVion_l,abs(trnsvFF_mlm(:,3)),'-.b', ...
       relVion_h,abs(trnsvFF_mhm(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: Magnetized Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([5.e-5,1.e+1])
text(2.e-5,5.e-0,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(1.55e-5,3.e-1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(5.e-5,1.e-2,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(7.e-4,1.e-3,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.05e-6,5.5e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.05e-6,2.25e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.05e-6,9.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
annotation('arrow',[.765 .8175],[.31 .25]);

%
% B=600 Gs (Three types of interaction; "Meshkov" approach):
%       
figure(400)
loglog(relVion_s,abs(trnsvFF_fs(:,2)),'-r',relVion_l,abs(trnsvFF_fl(:,2)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,2)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,2)),'--r',relVion_l,abs(trnsvFF_al(:,2)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msm(:,2)),'-.r',relVion_l,abs(trnsvFF_mlm(:,2)),'-.b', ...
       relVion_h,abs(trnsvFF_mhm(:,2)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([1.5e-5,2.e-1])
text(2.e-6,1.2e-1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(2),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-4,1.e-3,'Fast','Color','k','FontSize',16)
text(3.e-5,2.e-2,'Adiabatic','Color','k','FontSize',16)
text(5.e-4,2.e-3,'Magnetized','Color','k','FontSize',16)
text(1.25e-5,1.05e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-5,5.e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-5,2.5e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% Transverse friction force (adiabatic):
%       
figure(405)
loglog(relVion_s,abs(trnsvFF_as(:,1)),'-r',relVion_l,abs(trnsvFF_al(:,1)),'-b', ...
       'LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: Adiabatic Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([2.e-4,2.e-2])
text(2.e-5,1.5e-2,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(1.e-5,4.e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,2.75e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)

%
% B=100 Gs (Three types of interaction; "Meshkov" approach):
%       
figure(410)
loglog(relVion_s,abs(trnsvFF_fs(:,3)),'-r',relVion_l,abs(trnsvFF_fl(:,3)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,3)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,3)),'--r',relVion_l,abs(trnsvFF_al(:,3)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msm(:,3)),'-.r',relVion_l,abs(trnsvFF_mlm(:,3)),'-.b', ...
       relVion_h,abs(trnsvFF_mhm(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([2.e-5,2.e-1])
text(2.e-6,1.2e-1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(3),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.2e-4,1.5e-3,'Fast','Color','k','FontSize',16)
text(3.e-5,2.e-6,'Adiabatic','Color','k','FontSize',16)
text(3.5e-4,4.e-4,'Magnetized','Color','k','FontSize',16)
text(1.e-5,1.e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,5.5e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,3.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% Transverse friction force (fast; three values of magnetic field; "Meshkov" approach):
%       
figure(415)
loglog(relVion_s,abs(trnsvFF_fs(:,1)),'-r',relVion_l,abs(trnsvFF_fl(:,1)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_fs(:,2)),'--r',relVion_l,abs(trnsvFF_fl(:,2)),'--b', ...
       relVion_h,abs(trnsvFF_fh(:,2)),'--m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_fs(:,3)),'-.r',relVion_l,abs(trnsvFF_fl(:,3)),'-.b', ...
       relVion_h,abs(trnsvFF_fh(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: Fast Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([1.e-5,2.e-5])
text(2.e-5,1.e-5,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(1.8e-5,1.e-2,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-4,3.e-4,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.5e-4,2.e-3,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,6.e-5,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,3.e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,1.5e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
arrw=annotation('arrow',[.5 .425],[.42 .46]);

%
% Transverse friction force (all; "Meshkov" approach):
%       
figure(420)
loglog(relVion_s,abs(trnsvFF_msm(:,1)),'-r',relVion_l,abs(trnsvFF_mlm(:,1)),'-b', ...
       relVion_h,abs(trnsvFF_mhm(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_msm(:,2)),'-r',relVion_l,abs(trnsvFF_mlm(:,2)),'-b', ...
       relVion_h,abs(trnsvFF_mhm(:,2)),'-m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msm(:,3)),'-r',relVion_l,abs(trnsvFF_mlm(:,3)),'-b', ...
       relVion_h,abs(trnsvFF_mhm(:,3)),'-m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_as(:,1)),'--r',relVion_l,abs(trnsvFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_fs(:,1)),'-.r',relVion_l,abs(trnsvFF_fl(:,1)),'-.b', ...
       relVion_h,abs(trnsvFF_fh(:,1)),'-.m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_fs(:,2)),'-.r',relVion_l,abs(trnsvFF_fl(:,2)),'-.b', ...
       relVion_h,abs(trnsvFF_fh(:,2)),'-.m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_fs(:,3)),'-.r',relVion_l,abs(trnsvFF_fl(:,3)),'-.b', ...
       relVion_h,abs(trnsvFF_fh(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp','Color','m','FontSize',13)
xlim(xLimit)
ylim([1.e-5,1.e+1])
text(2.e-5,5.e-0,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(3.e-5,2.e-2,'Adiabatic','Color','k','FontSize',16)
text(1.5e-5,7.e-1,'Magnetized:','Color','k','FontSize',16)
text(1.5e-5,3.e-1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(8.e-4,7.e-1,'Magnetized:','Color','k','FontSize',16)
text(8.e-4,3.e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(3.e-4,9.e-4,'Magnetized:','Color','k','FontSize',16)
text(3.e-4,4.e-4,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,1.25e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,5.e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,2.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
annotation('arrow',[.8 .825],[.7 .49]);
dim=[.39 .39 .15 .09];
annotation('ellipse',dim,'Color','k');
text(3.e-5,1.75e-3,'Fast','Color','k','FontSize',16)

%
% Transverse friction force (total; "Meshkov" approach):
%       
trnsvFF_tsm=zeros(nVion_s,nField);
trnsvFF_tlm=zeros(nVion_l,nField);
trnsvFF_thm=zeros(nVion_h,nField);

for k=1:nField
    for i=1:nVion_s
        trnsvFF_tsm(i,k)=trnsvFF_msm(i,k)+trnsvFF_fs(i,k)+trnsvFF_as(i);  
    end
    for i=1:nVion_l
        trnsvFF_tlm(i,k)=trnsvFF_mlm(i,k)+trnsvFF_fl(i,k)+trnsvFF_al(i);  
    end
    for i=1:nVion_h
        trnsvFF_thm(i,k)=trnsvFF_mhm(i,k)+trnsvFF_fh(i,k);  
    end
end

figure(430)
loglog(relVion_s,abs(trnsvFF_tsm(:,1)),'-.r',relVion_l,abs(trnsvFF_tlm(:,1)),'-.b', ...
       relVion_h,abs(trnsvFF_thm(:,1)),'-.m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_tsm(:,2)),'-r',relVion_l,abs(trnsvFF_tlm(:,2)),'-b', ...
       relVion_h,abs(trnsvFF_thm(:,2)),'-m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_tsm(:,3)),'--r',relVion_l,abs(trnsvFF_tlm(:,3)),'--b', ...
       relVion_h,abs(trnsvFF_thm(:,3)),'--m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title(['Total Transverse Friction Force F_\perp (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
xlim(xLimit)
ylim([2.e-4,4.e-0])
text(1.5e-5,3.5e-1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(7.e-4,1.e-2,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
annotation('arrow',[.75 .7],[.48 .55]);
text(1.5e-5,5.e-2,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
annotation('arrow',[.58 .7],[.575 .575]);
text(1.5e-5,2.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.5e-5,9.5e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.5e-5,4.5e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
%
% For comparison with figure 4 from [3] only: B=1000 Gs ("Meshkov" approach):
%
if (reference3flag ==1)
figure(8430)
loglog(relVion_s,abs(trnsvFF_tsm(:,4)),'-.r',relVion_l,abs(trnsvFF_tlm(:,4)),'-.b', ...
       relVion_h,abs(trnsvFF_thm(:,4)),'-.m','LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title(['Total Transverse Friction Force F_\perp (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
xlim(V0*.01*xLimit)
ylim([4.e-2,1.5e+2])
text(7.e-6,1.e+2,['B=',num2str(fieldB(4),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.5e-5,2.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.5e-5,9.5e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.5e-5,4.5e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
end
%
% B=3000 Gs (Three types of interaction; "Betacool" approach):
%       
figure(3901)
loglog(relVion_s,abs(trnsvFF_fs(:,1)),'-r',relVion_l,abs(trnsvFF_fl(:,1)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,1)),'--r',relVion_l,abs(trnsvFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msb(:,1)),'-.r',relVion_l,abs(trnsvFF_mlb(:,1)),'-.b', ...
       relVion_h,abs(trnsvFF_mhb(:,1)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([1.e-5,1.e+2])
text(2.e-6,4.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-4,6.e-4,'Fast','Color','k','FontSize',16)
text(3.e-5,2.e-2,'Adiabatic','Color','k','FontSize',16)
text(1.65e-5,2.e-0,'Magnetized','Color','k','FontSize',16)
text(1.25e-5,1.5e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-5,5.e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-5,2.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% For comparison with figure 4 from [3] only: B=1000 Gs ("BETACOOL" approach):
%
if (reference3flag == 1)
figure(83901)
loglog(relVion_s,abs(trnsvFF_fs(:,4)),'-r',relVion_l,abs(trnsvFF_fl(:,4)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,4)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,4)),'--r',relVion_l,abs(trnsvFF_al(:,4)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msb(:,4)),'-.r',relVion_l,abs(trnsvFF_mlb(:,4)),'-.b', ...
       relVion_h,abs(trnsvFF_mhb(:,4)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(V0*.01*xLimit)
ylim([1.e-5,1.25e+2])
text(4.e-7,7.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(4),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-5,8.e-3,'Fast','Color','k','FontSize',16)
text(1.e-5,1.e-1,'Adiabatic','Color','k','FontSize',16)
text(3.e-6,4.e+0,'Magnetized','Color','k','FontSize',16)
text(2.e-6,1.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(2.e-6,2.25e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(2.e-6,7.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
 end
 
%
% Magnetized area (three values of magnetic field; "Betacool" approach):
%       
figure(3951)
loglog(relVion_s,abs(trnsvFF_msb(:,1)),'-r',relVion_l,abs(trnsvFF_mlb(:,1)),'-b', ...
       relVion_h,abs(trnsvFF_mhb(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_msb(:,2)),'--r',relVion_l,abs(trnsvFF_mlb(:,2)),'--b', ...
       relVion_h,abs(trnsvFF_mhb(:,2)),'--m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msb(:,3)),'-.r',relVion_l,abs(trnsvFF_mlb(:,3)),'-.b', ...
       relVion_h,abs(trnsvFF_mhb(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: Magnetized Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([7.e-5,1.e+2])
text(2.e-5,4.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(4.e-4,2.e-0,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,4.e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.6e-4,4.e-4,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.1e-6,2.e-2,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.1e-6,8.e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.1e-6,3.e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=600 Gs (Three types of interaction; "Betacool" approach):
%       
figure(4001)
loglog(relVion_s,abs(trnsvFF_fs(:,2)),'-r',relVion_l,abs(trnsvFF_fl(:,2)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,2)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,2)),'--r',relVion_l,abs(trnsvFF_al(:,2)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msb(:,2)),'-.r',relVion_l,abs(trnsvFF_mlb(:,2)),'-.b', ...
       relVion_h,abs(trnsvFF_mhb(:,2)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([1.5e-5,3.e+1])
text(2.e-6,1.25e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(2),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-4,8.e-4,'Fast','Color','k','FontSize',16)
text(3.e-5,2.e-2,'Adiabatic','Color','k','FontSize',16)
text(2.e-5,6.e-1,'Magnetized','Color','k','FontSize',16)
text(1.25e-5,1.5e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-5,6.e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-5,2.5e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=100 Gs (Three types of interaction; "Betacool" approach):
%       
figure(4101)
loglog(relVion_s,abs(trnsvFF_fs(:,3)),'-r',relVion_l,abs(trnsvFF_fl(:,3)),'-b', ...
       relVion_h,abs(trnsvFF_fh(:,3)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_as(:,3)),'--r',relVion_l,abs(trnsvFF_al(:,3)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msb(:,3)),'-.r',relVion_l,abs(trnsvFF_mlb(:,3)),'-.b', ...
       relVion_h,abs(trnsvFF_mhb(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([2.e-5,2.e-1])
text(2.e-6,1.2e-1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(3),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.2e-4,1.5e-3,'Fast','Color','k','FontSize',16)
text(3.e-5,2.e-2,'Adiabatic','Color','k','FontSize',16)
text(1.5e-4,3.e-4,'Magnetized','Color','k','FontSize',16)
text(1.e-5,1.e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,5.5e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,3.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% Transverse friction force (fast; three values of magnetic field):
% This type on interaction does not depend on "Meshkov" or "Betacool"
% approach. So, future figure 4151 is the same as 415!
%       

%
% Transverse friction force (all; "Betacool" approach):
%       
figure(4201)
loglog(relVion_s,abs(trnsvFF_msb(:,1)),'-r',relVion_l,abs(trnsvFF_mlb(:,1)),'-b', ...
       relVion_h,abs(trnsvFF_mhb(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_msb(:,2)),'-r',relVion_l,abs(trnsvFF_mlb(:,2)),'-b', ...
       relVion_h,abs(trnsvFF_mhb(:,2)),'-m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_msb(:,3)),'-r',relVion_l,abs(trnsvFF_mlb(:,3)),'-b', ...
       relVion_h,abs(trnsvFF_mhb(:,3)),'-m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_as(:,1)),'--r',relVion_l,abs(trnsvFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(trnsvFF_fs(:,1)),'-.r',relVion_l,abs(trnsvFF_fl(:,1)),'-.b', ...
       relVion_h,abs(trnsvFF_fh(:,1)),'-.m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_fs(:,2)),'-.r',relVion_l,abs(trnsvFF_fl(:,2)),'-.b', ...
       relVion_h,abs(trnsvFF_fh(:,2)),'-.m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_fs(:,3)),'-.r',relVion_l,abs(trnsvFF_fl(:,3)),'-.b', ...
       relVion_h,abs(trnsvFF_fh(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp','Color','m','FontSize',13)
xlim(xLimit)
ylim([1.e-5,1.e+2])
text(2.e-5,4.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(3.e-5,2.e-2,'Adiabatic','Color','k','FontSize',16)
text(6.e-4,3.e-0,'Magnetized:','Color','k','FontSize',16)
text(6.e-4,1.e-0,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.e-5,6.e-1,'Magnetized:','Color','k','FontSize',16)
text(1.e-5,2.5e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.e-4,9.e-4,'Magnetized:','Color','k','FontSize',16)
text(1.38e-4,3.5e-8,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,1.25e-4,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,5.e-5,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,2.e-5,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
annotation('arrow',[.8 .825],[.675 .47]);
dim=[.325 .3 .15 .1];
annotation('ellipse',dim,'Color','k');
text(1.25e-5,7.e-4,'Fast','Color','k','FontSize',16)

%
% Transverse friction force (total; "Betacool" approach):
%       
trnsvFF_tsb=zeros(nVion_s,nField);
trnsvFF_tlb=zeros(nVion_l,nField);
trnsvFF_thb=zeros(nVion_h,nField);

for k=1:nField
    for i=1:nVion_s
        trnsvFF_tsb(i,k)=trnsvFF_msb(i,k)+trnsvFF_fs(i,k)+trnsvFF_as(i);  
    end
    for i=1:nVion_l
        trnsvFF_tlb(i,k)=trnsvFF_mlb(i,k)+trnsvFF_fl(i,k)+trnsvFF_al(i);  
    end
    for i=1:nVion_h
        trnsvFF_thb(i,k)=trnsvFF_mhb(i,k)+trnsvFF_fh(i,k);  
    end
end

figure(4301)
loglog(relVion_s,abs(trnsvFF_tsb(:,1)),'-.r',relVion_l,abs(trnsvFF_tlb(:,1)),'-.b', ...
       relVion_h,abs(trnsvFF_thb(:,1)),'-.m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(trnsvFF_tsb(:,2)),'-r',relVion_l,abs(trnsvFF_tlb(:,2)),'-b', ...
       relVion_h,abs(trnsvFF_thb(:,2)),'-m','LineWidth',2)
loglog(relVion_s,abs(trnsvFF_tsb(:,3)),'--r',relVion_l,abs(trnsvFF_tlb(:,3)),'--b', ...
       relVion_h,abs(trnsvFF_thb(:,3)),'--m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title(['Total Transverse Friction Force F_\perp (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
xlim(xLimit)
ylim([2.e-4,3.e+1])
text(2.e-4,7.e-0,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.e-5,8.e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.5e-5,2.5e-2,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.5e-5,2.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.5e-5,9.5e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.5e-5,4.5e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% For comparison with figure 4 from [3] only: B=1000 Gs ("BETACOOL" approach):
%
if (reference3flag == 1)
figure(84301)
loglog(relVion_s,abs(trnsvFF_tsb(:,4)),'-.r',relVion_l,abs(trnsvFF_tlb(:,4)),'-.b', ...
       relVion_h,abs(trnsvFF_thb(:,4)),'-.m','LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('F_\perp, eV/m','Color','m','FontSize',16)
title(['Total Transverse Friction Force F_\perp (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
xlim(V0*.01*xLimit)
ylim([2.e-4,3.e+1])
text(2.e-4,7.e-0,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.e-5,8.e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.5e-5,2.5e-2,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.5e-5,2.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.5e-5,9.5e-4,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.5e-5,4.5e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
end
%--------------------------------------------------------------
%
%  F  i g u r e s   f o r   L o n g i t u d i n a l    F r i c t i o n    F o r c e
%
%--------------------------------------------------------------
%
% Longitudinal friction force for B=3000 Gs (Three types of interaction):
%       
%
% B=3000 Gs (Three types of interaction; "Meshkov" approach):
%       
figure(490)
loglog(relVion_s,abs(longFF_fs(:,1)),'-r',relVion_l,abs(longFF_fl(:,1)),'-b', ...
       relVion_h,abs(longFF_fh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_as(:,1)),'--r',relVion_l,abs(longFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_msm(:,1)),'-.r',relVion_l,abs(longFF_ml(:,1)),'-.b', ...
       relVion_h,abs(longFF_mh(:,1)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msb(nVion_s,1)),abs(longFF_ml(1,1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([6.e-4,1.e+2])
text(2.e-6,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(3.e-5,1.e-2,'Fast','Color','k','FontSize',16)
text(2.75e-5,1.25e-1,'Adiabatic','Color','k','FontSize',16)
text(1.75e-5,3.e-0,'Magnetized','Color','k','FontSize',16)
text(1.e-5,5.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,2.4e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,1.1e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% For [3] only: B=1000 Gs (Three types of interaction; "Meshkov" approach):
%
if (reference3flag == 1)
figure(8490)
loglog(relVion_s*V0*.01,abs(longFF_fs(:,4)),'-r',relVion_l*V0*.01,abs(longFF_fl(:,4)),'-b', ...
       relVion_h*V0*.01,abs(longFF_fh(:,4)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s*V0*.01,abs(longFF_as(:,4)),'--r',relVion_l*V0*.01,abs(longFF_al(:,4)),'--b', ...
       'LineWidth',2)
loglog(relVion_s*V0*.01,abs(longFF_msm(:,4)),'-.r',relVion_l*V0*.01,abs(longFF_ml(:,4)),'-.b', ...
       relVion_h*V0*.01,abs(longFF_mh(:,4)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s)*V0*.01,relVion_s(nVion_s)*V0*.01], ...
     [abs(longFF_msb(nVion_s,4)),abs(longFF_ml(1,4))],'-.k','LineWidth',1)   
xlabel('Ion Velocity, V_i. m/s','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(V0*.01*xLimit)
ylim([2.e-3,3.e+2])
text(1.e+2,1.5e+2,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(4),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(3.e-5,1.e-6,'Fast','Color','k','FontSize',16)
text(2.75e-5,1.25e-5,'Adiabatic','Color','k','FontSize',16)
text(1.75e-5,3.e-4,'Magnetized','Color','k','FontSize',16)
text(1.e-5,5.e-7,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,2.4e-7,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,1.1e-7,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
end
%
% Magnetized interaction (three values of magnetic field; "Meshkov" approach):
%       
figure(495)
loglog(relVion_s,abs(longFF_msb(:,1)),'-r',relVion_l,abs(longFF_ml(:,1)),'-b', ...
       relVion_h,abs(longFF_mh(:,1)),'-m','LineWidth',2)
grid on
hold on
%
% Values longFF_ml and longFF_mh udon't depend on magnetic field; so, values  
% longFF_ml(:,2) and longFF_mh(:,2)are multiplyed on 1.1 and values
% longFF_ml(:,3 ) and longFF_mh(:,3)are multiplyed on .9 that would be 
% distinguishable from the longFF_ml(:,1) and longFF_mh(:,1) of the figure:
%
loglog(relVion_s,abs(longFF_msb(:,2)),'--r',relVion_l,abs(1.1*longFF_ml(:,2)),'--b', ...
       relVion_h,abs(1.1*longFF_mh(:,2)),'--m','LineWidth',2)
loglog(relVion_l,abs(.9*longFF_ml(:,3)),'-.b', ...
      relVion_h,abs(.9*longFF_mh(:,3)),'-.m',relVion_s,abs(longFF_msb(:,3)),'-.r','LineWidth',2)
% loglog(relVion_l,abs(.9*longFF_ml(:,3)),'-.b', ...
%       relVion_h,abs(.9*longFF_mh(:,3)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msb(nVion_s,2)),abs(longFF_ml(1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: Magnetized Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([6.e-4,1.e+2])
text(2.e-5,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(1.5e-6,1.e+1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(3.e-6,1.e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
% text(7.e-4,1.e-3,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
%      'FontSize',16)
dim=[.575 .55 .2 .1];
annotation('ellipse',dim,'Color','k');
text(3.e-4,7.e-1,'All Fields','Color','k','FontSize',16)
text(1.1e-6,5.5e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.1e-6,2.3e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.1e-6,1.e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=600 Gs (Three types of interaction; "Meshkov" approach):
%       
figure(500)
loglog(relVion_s,abs(longFF_fs(:,2)),'-r',relVion_l,abs(longFF_fl(:,2)),'-b', ...
       relVion_h,abs(longFF_fh(:,2)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_as(:,2)),'--r',relVion_l,abs(longFF_al(:,2)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_msb(:,2)),'-.r',relVion_l,abs(longFF_ml(:,2)),'-.b', ...
       relVion_h,abs(longFF_mh(:,2)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msb(nVion_s,2)),abs(longFF_ml(1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([8.e-4,1.e+2])
text(2.e-6,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(2),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-5,2.e-2,'Fast','Color','k','FontSize',16)
text(2.5e-5,1.5e-1,'Adiabatic','Color','k','FontSize',16)
text(3.e-5,1.2e-0,'Magnetized','Color','k','FontSize',16)
text(1.25e-5,7.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-5,3.e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-5,1.25e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% Longitudinal friction force (adiabatic):
%       
figure(505)
loglog(relVion_s,abs(longFF_as(:,1)),'-r',relVion_l,abs(longFF_al(:,1)),'-b', ...
       'LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: Adiabatic Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([5.e-3,1.e-0])
text(2.e-5,7.e-1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(2.e-6,1.e-2,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(2.e-6,7.e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)

%
% B=100 Gs (Three types of interaction; "Meshkov" approach):
%       
figure(510)
loglog(relVion_s,abs(longFF_fs(:,3)),'-r',relVion_l,abs(longFF_fl(:,3)),'-b', ...
       relVion_h,abs(longFF_fh(:,3)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_as(:,3)),'--r',relVion_l,abs(longFF_al(:,3)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_msm(:,3)),'-.r',relVion_l,abs(longFF_ml(:,3)),'-.b', ...
       relVion_h,abs(longFF_mh(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([8.e-4,1.e+2])
text(2.e-6,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(3),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(7.e-5,1.1e-1,'Fast','Color','k','FontSize',16)
text(2.e-5,1.e-0,'Adiabatic','Color','k','FontSize',16)
text(1.e-5,1.e+1,'Magnetized','Color','k','FontSize',16)
text(1.e-5,7.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,3.e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,1.3e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% Longitudinal friction force (fast; three values of magnetic field; "Meshkov" approach):
%       
figure(515)
loglog(relVion_s,abs(longFF_fs(:,1)),'-r',relVion_l,abs(longFF_fl(:,1)),'-b', ...
       relVion_h,abs(longFF_fh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_fs(:,2)),'--r',relVion_l,abs(longFF_fl(:,2)),'--b', ...
       relVion_h,abs(longFF_fh(:,2)),'--m','LineWidth',2)
loglog(relVion_s,abs(longFF_fs(:,3)),'-.r',relVion_l,abs(longFF_fl(:,3)),'-.b', ...
       relVion_h,abs(longFF_fh(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: Fast Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([6.e-4,2.e-1])
text(2.e-5,1.2e-1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(2.5e-5,1.e-2,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(3.e-4,3.e-2,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.15e-6,2.e-2,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,2.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,1.3e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,8.5e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
arrw=annotation('arrow',[.59 .5],[.67 .75]);

%
% Longitudinal friction force (all; "Meshkov" approach):
%       
figure(520)
loglog(relVion_s,abs(longFF_msm(:,1)),'-r',relVion_l,abs(longFF_ml(:,1)),'-b', ...
       relVion_h,abs(longFF_mh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_msm(:,2)),'-.r',relVion_l,abs(1.1*longFF_ml(:,2)),'-.b', ...
       relVion_h,abs(1.1*longFF_mh(:,2)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_msm(:,3)),'-.r',relVion_l,abs(.9*longFF_ml(:,3)),'-.b', ...
       relVion_h,abs(.9*longFF_mh(:,3)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_as(:,1)),'--r',relVion_l,abs(longFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_fs(:,1)),'-.r',relVion_l,abs(longFF_fl(:,1)),'-.b', ...
       relVion_h,abs(longFF_fh(:,1)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_fs(:,2)),'-.r',relVion_l,abs(longFF_fl(:,2)),'-.b', ...
       relVion_h,abs(longFF_fh(:,2)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_fs(:,3)),'-.r',relVion_l,abs(longFF_fl(:,3)),'-.b', ...
       relVion_h,abs(longFF_fh(:,3)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msm(nVion_s,1)),abs(longFF_ml(1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}','Color','m','FontSize',13)
xlim(xLimit)
ylim([5.e-4,1.e+2])
text(2.e-5,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(2.5e-5,2.e-1,'Adiabatic','Color','k','FontSize',16)
text(2.e-6,2.e-0,'Magnetized:','Color','k','FontSize',16)
text(2.e-6,.9e-4,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,4.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,1.75e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,7.e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
% dim=[.35 .39 .15 .09];
% annotation('ellipse',dim,'Color','k');
% text(1.75e-5,5.e-6,'Fast','Color','k','FontSize',16)
dim=[.35 .29 .15 .19];
annotation('ellipse',dim,'Color','k');
text(2.e-5,5.e-2,'Fast','Color','k','FontSize',16)
text(2.5e-5,2.e-2,'All','Color','k','FontSize',16)
text(2.e-5,1.e-2,'Fields','Color','k','FontSize',16)
dim=[.5 .65 .3 .19];
annotation('ellipse',dim,'Color','k');
text(2.e-4,1.e+1,'Magnetized:','Color','k','FontSize',16)
text(4.e-4,4.e-0,'All Fields','Color','k','FontSize',16)

%
% Longitudinal friction force (total; "Meshkov" approach):
%       
longFF_tsm=zeros(nVion_s,nField);
longFF_tlm=zeros(nVion_l,nField);
longFF_thm=zeros(nVion_h,nField);

for k=1:nField
    for i=1:nVion_s
        longFF_tsm(i,k)=longFF_msm(i,k)+longFF_fs(i,k)+longFF_as(i,k);  
    end
    for i=1:nVion_l
        longFF_tlm(i,k)=longFF_ml(i,k)+longFF_fl(i,k)+longFF_al(i,k);  
    end
    for i=1:nVion_h
        longFF_thm(i,k)=longFF_mh(i,k)+longFF_fh(i,k);  
    end
end

figure(530)
loglog(relVion_s,abs(longFF_tsm(:,1)),'-.r',relVion_l,abs(longFF_tlm(:,1)),'-.b', ...
       relVion_h,abs(longFF_thm(:,1)),'-.m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_tsm(:,2)),'-r',relVion_l,abs(longFF_tlm(:,2)),'-b', ...
       relVion_h,abs(longFF_thm(:,2)),'-m','LineWidth',2)
loglog(relVion_s,abs(.9*longFF_tsm(:,3)),'--r',relVion_l,abs(longFF_tlm(:,3)),'--b', ...
       relVion_h,abs(longFF_thm(:,3)),'--m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title(['Total Longitudinal Friction Force F_{||} (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_tsm(nVion_s,2)),abs(longFF_tlm(1))],'-.k','LineWidth',1)   
xlim(xLimit)
ylim([1.e-2,2.5e+1])
text(1.5e-6,2e-0,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
% annotation('arrow',[.2 .3],[.58 .46]);
% annotation('arrow',[.2 .52],[.63 .8]);
text(7.e-5,1.5e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
annotation('arrow',[.47 .36],[.415 .415]);
% annotation('arrow',[.7 .85],[.41 .265]);
text(7.e-5,.85e-1,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
annotation('arrow',[.47 .325],[.355 .355]);
% annotation('arrow',[.7 .85],[.35 .25]);
text(6.e-6,4.e-2,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(6.e-6,2.25e-2,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(6.e-6,1.35e-2,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
dim=[.5 .6 .2 .124];
annotation('ellipse',dim,'Color','k');
text(1.1e-4,2.e-0,'All  Fields','Color','k','FontSize',16)

%
% For [3]: longitudinal friction force (total; "Meshkov" approach):
%       
figure(1530)
loglog(relVion_s*V0*.01,abs(longFF_tsm(:,4)),'-.r',relVion_l*V0*.01,abs(longFF_tlm(:,4)),'-.b', ...
       relVion_h*V0*.01,abs(longFF_thm(:,4)),'-.m','LineWidth',2)
grid on
hold on
xlabel('Ion Velocity, V_i, m/s}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title(['Total Longitudinal Friction Force F_{||} (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
plot([relVion_s(nVion_s)*V0*.01,relVion_l(1)*V0*.01], ...
     [abs(longFF_tsm(nVion_s,2)),abs(longFF_tlm(1))],'-.k','LineWidth',1)   
xlim([2.e-7*V0*.01,3.e-3*V0*.01])
ylim([3.e-6,2.e-2])
text(1.5e-6,2e-4,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
% annotation('arrow',[.2 .3],[.58 .46]);
% annotation('arrow',[.2 .52],[.63 .8]);
text(7.e-5,1.5e-5,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
annotation('arrow',[.47 .36],[.415 .415]);
% annotation('arrow',[.7 .85],[.41 .265]);
text(7.e-5,.85e-5,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
annotation('arrow',[.47 .325],[.355 .355]);
% annotation('arrow',[.7 .85],[.35 .25]);
text(6.e-6,4.e-6,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(6.e-6,2.25e-6,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(6.e-6,1.35e-6,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
dim=[.5 .6 .2 .124];
annotation('ellipse',dim,'Color','k');
text(1.1e-4,2.e-4,'All  Fields','Color','k','FontSize',16)

%
% B=3000 Gs (Three types of interaction; "Betacool" approach):
%       
figure(4901)
loglog(relVion_s,abs(longFF_fs(:,1)),'-r',relVion_l,abs(longFF_fl(:,1)),'-b', ...
       relVion_h,abs(longFF_fh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_as(:,1)),'--r',relVion_l,abs(longFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_msb(:,1)),'-.r',relVion_l,abs(longFF_ml(:,1)),'-.b', ...
       relVion_h,abs(longFF_mh(:,1)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: 3 types of Interaction', ...
      'Color','m','FontSize',13)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msb(nVion_s,1)),abs(longFF_ml(1))],'-.k','LineWidth',1)   
xlim(xLimit)
ylim([6.e-4,1.e+2])
text(2.e-6,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(1),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-4,3.e-2,'Fast','Color','k','FontSize',16)
text(2.7e-5,1.5e-1,'Adiabatic','Color','k','FontSize',16)
text(1.65e-5,3.e-0,'Magnetized','Color','k','FontSize',16)
text(1.25e-5,4.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-5,1.7e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-5,8.e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% Magnetized area (three values of magnetic field; "Betacool" approach):
%       
figure(4951)
loglog(relVion_s,abs(longFF_msb(:,1)),'-r',relVion_l,abs(longFF_ml(:,1)),'-b', ...
       relVion_h,abs(longFF_mh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_msb(:,2)),'--r',relVion_l,abs(1.1*longFF_ml(:,2)),'--b', ...
       relVion_h,abs(1.1*longFF_mh(:,2)),'--m','LineWidth',2)
% loglog(relVion_s,abs(.9*longFF_msb(:,3)),'-.r',relVion_l,abs(.9*longFF_ml(:,3)),'-.b', ...
%       relVion_h,abs(.9*longFF_mh(:,3)),'-.m','LineWidth',2)
loglog(relVion_l,abs(.9*longFF_ml(:,3)),'-.b', ...
       relVion_h,abs(.9*longFF_mh(:,3)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msb(nVion_s,2)),abs(longFF_ml(1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: Magnetized Case', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([7.e-4,1.e+2])
text(2.e-5,6.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(1.1e-6,1.e+1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,3.e-1,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
dim=[.5 .6 .2 .124];
annotation('ellipse',dim,'Color','k');
text(1.1e-4,2.e-0,'All  Fields','Color','k','FontSize',16)
text(2.e-6,5.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(2.e-6,2.25e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(2.e-6,1.e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=600 Gs (Three types of interaction; "Betacool" approach):
%       
figure(5001)
loglog(relVion_s,abs(longFF_fs(:,2)),'-r',relVion_l,abs(longFF_fl(:,2)),'-b', ...
       relVion_h,abs(longFF_fh(:,2)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_as(:,2)),'--r',relVion_l,abs(longFF_al(:,2)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_msb(:,2)),'-.r',relVion_l,abs(longFF_ml(:,2)),'-.b', ...
       relVion_h,abs(longFF_mh(:,2)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msb(nVion_s,2)),abs(longFF_ml(1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([8.e-4,1.e+2])
text(2.e-6,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(2),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-4,3.5e-2,'Fast','Color','k','FontSize',16)
text(2.5e-5,1.5e-1,'Adiabatic','Color','k','FontSize',16)
text(2.5e-5,1.5e-0,'Magnetized','Color','k','FontSize',16)
text(1.25e-5,6.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.25e-5,2.5e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.25e-5,1.15e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% B=100 Gs (Three types of interaction; "Betacool" approach):
%       
figure(5101)
loglog(relVion_s,abs(longFF_fs(:,3)),'-r',relVion_l,abs(longFF_fl(:,3)),'-b', ...
       relVion_h,abs(longFF_fh(:,3)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_as(:,3)),'--r',relVion_l,abs(longFF_al(:,3)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_msb(:,3)),'-.r',relVion_l,abs(longFF_ml(:,3)),'-.b', ...
       relVion_h,abs(longFF_mh(:,3)),'-.m','LineWidth',2)
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}: 3 types of Interaction', ...
      'Color','m','FontSize',13)
xlim(xLimit)
ylim([8.e-4,1.e+2])
text(2.e-6,5.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s,  B=',num2str(fieldB(3),'%6.1f'),' Gs'], ...
       'Color','k','FontSize',16)
text(1.e-4,5.e-2,'Fast','Color','k','FontSize',16)
text(2.e-5,1.e-0,'Adiabatic','Color','k','FontSize',16)
text(1.e-5,1.e+1,'Magnetized','Color','k','FontSize',16)
text(1.e-5,6.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,2.5e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,1.2e-3,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)

%
% Longitudinal friction force (fast; three values of magnetic field):
% This type on interaction does not depend on "Meshkov" or "Betacool"
% approach. So, future figure 4151 is the same as 415!
%       

%
% Longitudinal friction force (all; "Betacool" approach):
%       
figure(5201)
loglog(relVion_s,abs(longFF_msb(:,1)),'-r',relVion_l,abs(longFF_ml(:,1)),'-b', ...
       relVion_h,abs(longFF_mh(:,1)),'-m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_msb(:,2)),'-.r',relVion_l,abs(1.1*longFF_ml(:,2)),'-.b', ...
       relVion_h,abs(1.1*longFF_mh(:,2)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_msb(:,3)),'--r',relVion_l,abs(.9*longFF_ml(:,3)),'--b', ...
       relVion_h,abs(.9*longFF_mh(:,3)),'--m','LineWidth',2)
loglog(relVion_s,abs(longFF_as(:,1)),'--r',relVion_l,abs(longFF_al(:,1)),'--b', ...
       'LineWidth',2)
loglog(relVion_s,abs(longFF_fs(:,1)),'-.r',relVion_l,abs(longFF_fl(:,1)),'-.b', ...
       relVion_h,abs(longFF_fh(:,1)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_fs(:,2)),'-.r',relVion_l,abs(longFF_fl(:,2)),'-.b', ...
       relVion_h,abs(longFF_fh(:,2)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_fs(:,3)),'-.r',relVion_l,abs(longFF_fl(:,3)),'-.b', ...
       relVion_h,abs(longFF_fh(:,3)),'-.m','LineWidth',2)
loglog(relVion_s,abs(longFF_msb(:,2)),'-.r',relVion_l,abs(longFF_ml(:,2)),'-.b', ...
       relVion_h,abs(longFF_mh(:,2)),'-.m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_l(1)], ...
     [abs(longFF_msb(nVion_s,2)),abs(longFF_ml(1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title('Longitudinal Friction Force F_{||}','Color','m','FontSize',13)
xlim(xLimit)
ylim([6.e-4,1.e+2])
text(2.e-5,6.e+1,['V_{e0}=',num2str(mantV0,'%4.2f'),'\cdot10^{',powV0, ...
     '} cm/s'],'Color','k','FontSize',16)
text(2.5e-5,1.5e-1,'Adiabatic','Color','k','FontSize',16)
text(3.e-6,2.e+1,'Magnetized:','Color','k','FontSize',16)
text(3.e-6,1.e+1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,2.e-0,'Magnetized:','Color','k','FontSize',16)
text(1.e-5,1.e-0,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.e-5,4.e-3,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.e-5,1.75e-3,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.e-5,8e-4,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
% dim=[.35 .3 .15 .1];
% annotation('ellipse',dim,'Color','k');
% text(3.e-5,1.2e-2,'Fast','Color','k','FontSize',16)
dim=[.35 .29 .15 .175];
annotation('ellipse',dim,'Color','k');
text(2.e-5,4.5e-2,'Fast','Color','k','FontSize',16)
text(2.5e-5,2.2e-2,'All','Color','k','FontSize',16)
text(2.15e-5,1.e-2,'Fields','Color','k','FontSize',16)
dim=[.5 .65 .3 .19];
annotation('ellipse',dim,'Color','k');
text(2.e-4,1.e+1,'Magnetized:','Color','k','FontSize',16)
text(4.e-4,4.e-0,'All Fields','Color','k','FontSize',16)

%
% Longitudinal friction force (total; "Betacool" approach):
%       
longFF_tsb=zeros(nVion_s,3);
longFF_tlb=zeros(nVion_l,3);
longFF_thb=zeros(nVion_h,3);

for k=1:3
    for i=1:nVion_s
        longFF_tsb(i,k)=longFF_msb(i,k)+longFF_fs(i,k)+longFF_as(i);  
    end
    for i=1:nVion_l
        longFF_tlb(i,k)=longFF_ml(i,k)+longFF_fl(i,k)+longFF_al(i);  
    end
    for i=1:nVion_h
        longFF_thb(i,k)=longFF_mh(i,k)+longFF_fh(i,k);  
    end
end

figure(5301)
loglog(relVion_s,abs(longFF_tsb(:,1)),'-.r',relVion_l,abs(longFF_tlb(:,1)),'-.b', ...
       relVion_h,abs(longFF_thb(:,1)),'-.m','LineWidth',2)
grid on
hold on
loglog(relVion_s,abs(longFF_tsb(:,2)),'-r',relVion_l,abs(longFF_tlb(:,2)),'-b', ...
       relVion_h,abs(longFF_thb(:,2)),'-m','LineWidth',2)
loglog(relVion_s,abs(longFF_tsb(:,3)),'--r',relVion_l,abs(longFF_tlb(:,3)),'--b', ...
       relVion_h,abs(longFF_thb(:,3)),'--m','LineWidth',2)
plot([relVion_s(nVion_s),relVion_s(nVion_s)], ...
     [abs(longFF_tsb(nVion_s,3)),abs(longFF_tsb(nVion_s,1))],'-.k','LineWidth',1)   
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title(['Total Longitudinal Friction Force F_{||} (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
xlim(xLimit)
ylim([1.e-2,3.e+1])
text(1.2e-6,1.e+1,['B=',num2str(fieldB(1),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(2.e-5,1.e-0,['B=',num2str(fieldB(2),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(3.e-5,2.e-1,['B=',num2str(fieldB(3),'%6.1f'),' Gs'],'Color','k', ...
     'FontSize',16)
text(1.5e-5,5.e-2,'Area "Superlow": V_i < \DeltaV_{e||} < \DeltaV_{e\perp}', ...
     'Color','r','FontSize',13)
text(1.5e-5,2.75e-2,' Area "Low":          \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','b','FontSize',13)
text(1.5e-5,1.6e-2,'Area "High":         \DeltaV_{e||} < V_i < \DeltaV_{e\perp}', ...
     'Color','m','FontSize',13)
dim=[.5 .65 .3 .09];
annotation('ellipse',dim,'Color','k');
text(3.e-4,3.e-0,'All Fields','Color','k','FontSize',16)

 %==================================================
 %
 % For fitting attempts
 %
 %==================================================

 longFF_3000tb=zeros(nVion_s+nVion_l+nVion_h,3);
relVion_t=zeros(nVion_s+nVion_l+nVion_h,1);
for n=1:nVion_s
   relVion_t(n)=relVion_s(n);
   longFF_3000tb(n,1)=longFF_tsb(n,1); 
end

for n=1:nVion_l
   relVion_t(nVion_s+n)=relVion_l(n);
   longFF_3000tb(nVion_s+n,1)=longFF_tlb(n,1); 
end
for n=1:nVion_h
   relVion_t(nVion_s+nVion_l+n)=relVion_h(n);
   longFF_3000tb(nVion_s+nVion_l+n,1)=longFF_thb(n,1); 
end

figure(53011)
loglog(relVion_t,abs(longFF_3000tb(:,1)),'-r','LineWidth',2)
grid on
hold on
xlabel('Relative Ion Velocity, V_{i||}/V_{e0}','Color','m','FontSize',16)
ylabel('F_{||}, eV/m','Color','m','FontSize',16)
title(['Total Longitudinal Friction Force F_{||} (V_{e0}=',num2str(mantV0,'%4.2f'), ...
       '\cdot10^{',powV0,'} cm/s)'],'Color','m','FontSize',12)
xlim(xLimit)
ylim([1.e-6,5.e-3])
 
 
