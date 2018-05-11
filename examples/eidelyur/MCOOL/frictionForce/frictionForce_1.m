%
% Physical constants:
%
m_e=9.10938356e-28;         % electron mass, g
q_e=4.803204673e-10;        % electron charge, sqrt(g*cm^3/sec^2)
cLight=2.99792458e10;       % speed of light, cm/sec
eVtoErg=1.6021766208e-12;   % 1 eV = 1.6...e-12 erg
CtoPart=2.99792458e9;       % 1 C = 1 A*sec = 2.9...e9 particles

pi=3.14159265358


% 
% Electron beam parameters (Table III):
%
Ekin=3.0e4;                   % kinetic energy, eV
curBeam=0.5;                  % current density, A/cm^2
dBeam=3.0;                    % beam diameter, cm
angSpread=3.0;                % angular spread, mrad
trnsvT=0.5;                   % transversal temperature, eV
longT=2.0e-4;                 % longitudinal temperature, eV
fieldB=1.e3*[3. 0.6 0.1]';    % magnetic field, Gs
omega_p=1.0e9;                % plasma frequency, 1/sec

% 
% Calculated parameters of the electron beam:
%
V0=sqrt(2.*Ekin*eVtoErg/m_e)                % average velocity, cm/sec
rmsTrnsvVe=sqrt(2.*trnsvT*eVtoErg/m_e)      % RMS transversal velocity, cm/sec
rmsLongVe=sqrt(2.*longT*eVtoErg/m_e)        % RMS longitudinal velocity, cm/sec
% dens=curBeam*CtoPart/V0                   % density, 1/cm^3
% omega=sqrt(4.*pi*dens*q_e^2/m_e)          % plasma frequency, 1/sec
cyclFreq=q_e*fieldB/(m_e*cLight)            % cyclotron frequency, 1/sec
rmsRoLarm=rmsTrnsvVe*cyclFreq.^(-1)         % RMS Larmor radius, cm
dens=omega_p^2*m_e/(4.*pi*q_e^2)            % density, 1/cm^3
%
% Formfactor ffForm for friction force:
%
% ffForm = 2*pi*dens*q_e^4/(m_e*V0^2)=
%        = 0.5*omega_p^2*q_e^2/V0^2;
%
% Dimension of ffForm is  g*cm/sec^2;
% 1 MeV/m = 1.e6*eVtoErg*.01 g*cm/sec^2 ==>
MeV_mToErg_cm=1.e4*eVtoErg
ffForm=.5*omega_p^2*q_e^2/V0^2/MeV_mToErg_cm      % MeV/m

%
% Relative velocities of electrons:
%
relVeTrnsv=V0/rmsTrnsvVe 
relVeLong=V0/rmsLongVe
%
% Scanning over relative ion velocity from 1e-6 till 1e-2: 
%
logVionMin=-6.;
logVionMax=-2.;
nVion=100;
stepLogVion=(logVionMax-logVionMin)/(nVion-1);
       
relVion=zeros(nVion);       % relative ion velocity (related V0)
%    
% Number interactions between ion electron during its passing near ion:
%
N_l=zeros(nVion);             % low 
N_s=zeros(nVion);             % superlow 
trnsvFF_sl=zeros(nVion,3);    % transversal friction force (superlow)
trnsvFF_l=zeros(nVion,3);     % transversal friction force (low)
trnsvFF_h=zeros(nVion,3);     % transversal friction force (high)
longFF_sl=zeros(nVion,3);     % longitufinal friction force (superlow)
lonhFF_l=zeros(nVion,3);      % longitudinal friction force (low)
longFF_h=zeros(nVion,3);      % lomgitudinal friction force (high)
%
% Indices: first describes the type of collision 
% and second characrerises ion velovity relative
% the electron.
%
% Values of the indices:
% First = f means fast), or 
%       = a (adiabatical), or 
%       = m (magnetized);
% Second = s (means superlow), or
%        = l (low), or
%        = h (high)
%
% Coulomb logarithms:
%
CL_as=log(2.*rmsTrnsvVe/rmsLongVe);                           % value approx 5.
CL_al=log(2.);                                                % value from log2 till CL_as
CL_fs=log(m_e*rmsLongVe*rmsTrnsvVe^2/q_e^2*cyclFreq.^(-1));   % 3 values approx 5.
CL_ms=log(.5*(3./dens)^(1./3.)*rmsRoLarm.^(-1));              % 3 values approx 2.
CL_fl=zeros(nVion,3);                                         % nVion x 3 values approx 5.
CL_fh=zeros(nVion,3);                                         % nVion x 3 values approx 10.
CL_ml=zeros(nVion,3);                                         % nVion x 3 values approx from 3.0 till 1.0
CL_mh=zeros(nVion,3);                                         % nVion x 3 values  approx 2.

for n=1:nVion
    curLogVion=logVionMin+stepLogVion*(n-1);
    relVion(n)=10^curLogVion;
    absVion=V0*relVion(n);
    if (absVion > rmsLongVe)
        N_l(n)=1+round(rmsTrnsvVe/(pi*absVion));
    else
        N_s(n)=1+round(rmsTrnsvVe/(pi*rmsLongVe));
    end
%
% Coulomb logarithms:
%
    for k=1:3
        CL_mh(n,k)=log(absVion*cyclFreq(k)/(2.*rmsLongVe*omega_p));
        CL_fh(n,k)=log(m_e*absVion^3/(q_e^2*cyclFreq(k)));
        num=min(absVion/omega_p,(3./dens)^(1./3.));
        CL_ml(n,k)=log(.5*num/rmsRoLarm(k));
        CL_fl(n,k)=CL_fh(n,k);
% 
% Friction forces:
%
        trnsvFF_sl(n,k)=abs(ffForm*relVion(n)* ...
                     (CL_ms(k)/relVion(n)^3+2*relVeTrnsv^3*(CL_fs(k)+N_s(n)*CL_as)));
        trnsvFF_l(n,k)=abs(ffForm*relVion(n)* ...
                     (CL_ml(n,k)/relVion(n)^3+2*relVeTrnsv^3*(CL_fl(n,k)+N_l(n)*CL_al)));
        trnsvFF_h(n,k)=abs(ffForm*(2.*CL_fh(n,k)+CL_mh(n,k))/relVion(n)^2);
        longFF_sl(n,k)=abs(ffForm*relVion(n)* ...
                     (CL_ms(k)/relVion(n)^3+2*relVeTrnsv^2*relVeLong*(CL_fs(k)+N_s(n)*CL_as)));
        longFF_l(n,k)=abs(ffForm*relVion(n)* ...
                     (2./relVion(n)^3+2*relVeTrnsv^2*(CL_fl(n,k)+N_l(n)*CL_al)));
        longFF_h(n,k)=abs(ffForm*(2.*CL_fh(n,k)+2.)/relVion(n)^2);
    end
end

figure(10)
semilogx(relVion,CL_mh(:,1),'-r',relVion,CL_mh(:,2),'-b',relVion,CL_mh(:,3),'-m','LineWidth',2)
grid on
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{mh}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{mh}','Color','m','FontSize',16)
legend(['B=',num2str(fieldB(1),'%6.1f'),' kG'],['B=',num2str(fieldB(2),'%6.1f'),' kG'], ...
       ['B=',num2str(fieldB(3),'%6.1f'),' kG'],'Location','s')

figure(20)
semilogx(relVion,CL_fh(:,1),'-r',relVion,CL_fh(:,2),'-b',relVion,CL_fh(:,3),'-m','LineWidth',2)
grid on
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{fh} and CL_{fl}','Color','m','FontSize',16)
title('Coulomb Logarithms CL_{fh} and CL_{fl}','Color','m','FontSize',16)
legend(['B=',num2str(fieldB(1),'%6.1f'),' kG'],['B=',num2str(fieldB(2),'%6.1f'),' kG'], ...
       ['B=',num2str(fieldB(3),'%6.1f'),' kG'],'Location','s')

figure(30)
semilogx(relVion,CL_ml(:,1),'-r',relVion,CL_ml(:,2),'-b',relVion,CL_ml(:,3),'-m','LineWidth',2)
grid on
xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','Color','m','FontSize',16)
ylabel('CL_{ml}','Color','m','FontSize',16)
title('Coulomb Logarithm CL_{ml}','Color','m','FontSize',16)
legend(['B=',num2str(fieldB(1),'%6.1f'),' kG'],['B=',num2str(fieldB(2),'%6.1f'),' kG'], ...
       ['B=',num2str(fieldB(3),'%6.1f'),' kG'],'Location','s')


figure(100)
loglog(relVion,trnsvFF_sl(:,1),'-r','LineWidth',2)
grid on
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','FontSize',16,'Color','m')
ylabel('F_\perp, MeV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp','Color','m','FontSize',16)
   
figure(110)
loglog(relVion,trnsvFF_l(:,1),'-b','LineWidth',2)
grid on
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','FontSize',16,'Color','m')
ylabel('F_\perp, MeV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp','Color','m','FontSize',16)

figure(120)
loglog(relVion,trnsvFF_h(:,1),'-m','LineWidth',2)
grid on
% xlabel('Relative Ion Velocity, V_{i\perp}/V_{e0}','FontSize',16,'Color','m')
ylabel('F_\perp, MeV/m','Color','m','FontSize',16)
title('Transverse Friction Force F_\perp','Color','m','FontSize',16)




