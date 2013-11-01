function []=stt_ver4()

close all;

% Velocities of polymerization and depolymerization in mic per min
vg1=7.33;
vg2=1.07;
vs1=18.98;
vs2=1.07;

% Transition frequencies for the plus end
pgp1=0.001*60;
pgs1=0.013*60;
ppg1=0.055*60;
pps1=0.062*60;
psg1=0.027*60;
psp1=0.003*60;

% Transition frequencies for the minus end
pgp2=0.014*60;
pgs2=0.025*60;
ppg2=0.003*60;
pps2=0.032*60;
psg2=0.001*60;
psp2=0.013*60;

MTs=500:1000;
mg1=MTs.*0.15;
ms1=MTs.*0.84;
mp1=MTs-(mg1+ms1);
mg2=MTs.*0.10;
ms2=MTs.*0.95;
mp2=MTs-(mg2+ms2);

% Velocities relationships at steady state
y1=vg1+vs1; % sum of the velocities at the plus end
y2=vg2+vs2; % sum of the velocities at the minus end
y3=vg2-vs1; % difference between growing vel at the minus end and shortening vel at the plus end

y1=1./mg1.*(-y3.*MTs + y2.*ms2 + vg2.*mp2 - vs1.*mp1);
mean(y1-(vg1+vs1))

(vg1.*mg1+vg2.*mg2)-(vs1.*ms1+vs2.*ms2)

% + end
N1=(pgs1+pgp1)*mg1-(ppg1*mp1+psg1*ms1);
E1=((-ppg1-pps1-(pps1*pgp1)/pgs1)*mp1+(psp1+(psg1+psp1)*(pgp1/pgs1))*ms1)*(-pgs1/pgp1);
% - end
E2=((pps2*psg2)/ppg2-psp2-psg2+(2*pgs2*psg2)/(pgp2+pgs2))*ms2;
N2=-(((-ppg2-pps2-(pps2*pgp2)/pgs2)*mp2+(psp2+(psg2+psp2)*(pgp2/pgs2))*ms2)+E2*(pgp2/pgs2));

test=-(pgs2+pgp2)*mg2+ppg2*mp2+psg2*ms2;

figure
plot(MTs,N1,'r');
hold on;
plot(MTs,E1,'g');
hold on;
plot(MTs,E2,'b');
hold on;
plot(MTs,N2,'m');
hold on;
plot(MTs,test,'y');

set(gca,'fontweight','b','fontsize',16);
xlabel('MTs','fontweight','b','fontsize',16);
ylabel('Nucleation / Extinction (min^{-1})','fontweight','b','fontsize',16);
legend('Nucleation','Extinction');