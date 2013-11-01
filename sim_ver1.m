function []=sim_ver1()

%clc;
%clear all;
%close all;

% ---------------------------- Parameters ----------------------------

% Parameters for nucleation
k1=0.01;       % nuc rate slope
k2=0.005;       % nuc rate ordinate

% Parameters for assembly
k3=1750e-06;  % pol vel slope (+ end)
k4=1e-06;    % pol vel ordinate (+ end)
k5=18e-06;    % depol vel (+ end)
k6=750e-06;  % pol vel slope (- end)
k7=2.5e-06;    % pol vel ordinate (- end)
k8=1e-06;     % depol vel (- end)

% Parameters for transitions at the MT plus end
k9=0.0002;      % gp transition slope
k10=0.002;     % gp transition ordinate
k11=0.001;     % gs transition slope
k12=0.015;     % gs transition ordinate

k13=0.002;     % pg transition slope
k14=0.040;     % pg transition ordinate
k15=0.0017;     % ps transition slope
k16=0.07;     % ps transition ordinate

k17=0.002;     % sg transition slope
k18=0.015;     % sg transition ordinate
k19=0.0005;     % sp transition slope
k20=0.001;     % sp transition ordinate

% Parameters for transitions at the MT minus end
k21=0.0018;      % gp transition slope
k22=0.022;     % gp transition ordinate
k23=0.001;     % gs transition slope
k24=0.030;     % gs transition ordinate

k25=0.0005;     % pg transition slope
k26=0.001;     % pg transition ordinate
k27=0.0015;     % ps transition slope
k28=0.04;     % ps transition ordinate

k29=0.0002;     % sg transition slope
k30=0.0005;     % sg transition ordinate
k31=0.001;     % sp transition slope
k32=0.005;     % sp transition ordinate

% ------------------------------- Model -------------------------------

ndimer=1634;              % Dimers per micrometer of microtubule
avg=6.022e23;             % Avogadro's number
volume=3000;              % Volume of cell in cubic micrometers
initial_tubulin=10;          % Initial tubulin concentration in micromolar

%initial_free_dimers=init_tubulin*1e-6*avg*volume*1e-18/0.001;
fmt=0:.01:initial_tubulin;  % Concentration range for subunits in micromolar
free_tubulin_dimers=fmt*1e-6*avg*volume*1e-18/0.001; % # of free subunits
ft=free_tubulin_dimers/ndimer*1e-6;

% Sets nucleation
N=k1*fmt-k2;

% Linear models for assembly at both MT ends
vg_p=k3*ft-k4;
vs_p=ones(1,length(ft))*k5;
vg_m=k6*ft-k7;
vs_m=ones(1,length(ft))*k8;

% Linear models for the transition frequencies of both plus and minus end
tfp={'' '-k9*fmt+k10' '-k11*fmt+k12'; 'k13*fmt+k14' '' '-k15*fmt+k16'; 'k17*fmt+k18' 'k19*fmt+k20' ''};
tfm={'' '-k21*fmt+k22' '-k23*fmt+k24'; 'k25*fmt+k26' '' '-k27*fmt+k28'; 'k29*fmt+k30' 'k31*fmt+k32' ''};

% Obtain transition probabilities
tps=zeros(12,length(fmt));
tm_aux={'tfp' 'tfm'}; % Sets the two possible transition frequencies matrices

count=1;
for k=1:2
    for m=1:3
        for n=1:3
            if ~(m==n)
                tps(count,:)=eval(eval(strcat(tm_aux{k},'{m,n}')));
                count=count+1;
            end
        end
    end
end

% --------------- Plot1 (Velocities and Nucleation) -----------------

% Velocities of polymerization and depolymerization on both ends
YMatrix(1,:)=vg_p*1e6;
YMatrix(2,:)=vs_p*1e6;
YMatrix(3,:)=vg_m*1e6;
YMatrix(4,:)=vs_m*1e6;

figure1=figure;

% Plots velocities of polymerization and depolymerization
axes1=axes('Parent',figure1);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(fmt,YMatrix,'Parent',axes1);
legend();

% Plots nucleation rate
axes2 = axes('Parent',figure1,...
    'YAxisLocation','right',...
    'YColor','m',...
    'Color','none');
hold(axes2,'all');
plot(fmt,N,'Color','m','Parent',axes2);

% Create ylabels
ylabel(axes1,'Velocity (micrometer per min)','FontWeight','bold','FontSize',16);
ylabel(axes2,'Nucleation Probability','FontWeight','bold','FontSize',16);
xlabel('Free tubulin (micromolar)','FontWeight','bold','FontSize',16);

% ------------ Plot2 and Plot3 (Transition Probabilities) -------------

% Plus end

figure2=figure;

% Plots transition probabilities
axes2=axes('Parent',figure2);
box(axes2,'on');
grid(axes2,'on');
hold(axes2,'all');
plot(fmt,tps(1:6,:),'Parent',axes2);
legend('Pgp','Pgs','Ppg','Pps','Psg','Psp');

% Create ylabels
ylabel(axes2,'Tansition probabilities (+ end)','FontWeight','bold','FontSize',16);
xlabel('Free tubulin (micromolar)','FontWeight','bold','FontSize',16);

% Minus end

figure3=figure;

% Plots transition probabilities
axes3=axes('Parent',figure3);
box(axes3,'on');
grid(axes3,'on');
hold(axes3,'all');
plot(fmt,tps(7:12,:),'Parent',axes3);
legend('Pgp','Pgs','Ppg','Pps','Psg','Psp');

% Create ylabels
ylabel(axes3,'Tansition probabilities (- end)','FontWeight','bold','FontSize',16);
xlabel('Free tubulin (micromolar)','FontWeight','bold','FontSize',16);

end














