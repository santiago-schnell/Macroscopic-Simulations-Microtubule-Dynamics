function []=stt_ver3()

clc;
clear all;
close all;

% Velocities of polymerization and depolymerization in mic per min
vg_p=7.33;
vg_m=1.07;
vs_p=18.98;
vs_m=1.07;

% Velocities relationships at steady state
y1=vg_p+vs_p; % sum of the velocities at the plus end
y2=vg_m+vs_m; % sum of the velocities at the minus end
y3=vg_m-vs_p; % difference between growing vel at the minus end and shortening vel at the plus end

% Main functions
exp1()

% This function shows the relationships between #MTs growing and shortening
% at steady state as a function of the total # of MTs and # of MTs
% shortening at the minus end
    function exp1()
        mp_p=0:0.01:0.5; % fraction of MTs pausing at the plus end
        mp_m=0:0.01:0.5; % fraction of MTs pausing at the minus end
        ms_m=0:0.1:1; % fraction of MTs shortening at the minus end
        T=2500;
        
        for i=1:length(ms_m)
            [X,Y]=meshgrid(mp_m,mp_p);
            mg_p=1/y1*(y2*ms_m(i)-y3+vg_m.*X-vs_p.*Y); % obtains the fraction of MTs growing at the plus end
            figure(1);
            hold on;
            surf(X,Y,mg_p);
            view(-75,40);
            
            % I used the lines below just to confirm that my calculations were correct and
            % that this is indeed the expression holding the steady state of the system
            mg_m=T-ms_m(i)*T-Y*T;
            ms_p=T-mg_p*T-X*T;
            pol=vg_p.*mg_p*T+vg_m.*mg_m;
            depol=vs_p.*ms_p+vs_m.*ms_m(i)*T;
            res(i)=mean(mean(pol-depol));
        end
        
        %Difference between polymerization and depolymerization at steady
        %state for different fraction of MTs shortening at the minus end
        res
        
        figure(1);
        set(gca,'fontweight','b','fontsize',16);
        xlabel('mp_p','fontweight','b','fontsize',16);
        ylabel('mp_m','fontweight','b','fontsize',16);
        zlabel('mg_p');
        colorbar;
        grid on;
        
    end

% Transition matrix
% tmv=[-(pgs+pgp) ppg psg; pgp -(ppg+pps) psp; pgs pps -(psg+psp)];

%S=solve('(-(pgs+pgp)*mg)+(pgp*mp)+(psg*ms)=-N',mp);
%mp_p=-(N - mg_p*(pgp_p + pgs_p) + ms_p*psg_p)/pgp_p;

%S=solve('(-(pgs+pgp)*mg)+(pgp*mp)+(psg*ms)=0',mp);
%mp_m=(mg_m*(pgp_m + pgs_m) - ms_m*psg_m)/pgp_m;

end














