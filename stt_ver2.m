function []=stt_ver2()

clc;
clear all;
close all;

% Velocities of polymerization and depolymerization in mic per min
vg_p=7.33;
vg_m=1.07;
vs_p=18.98;
vs_m=1.07;

% Catastrophe and rescue frequencies at the minus end of a MT
pgs_m=0.035*60;
psg_m=0.015*60;

% Velocities relationships at steady state
y1=vg_p+vs_p; % sum of the velocities at the plus end 
y2=vg_m+vs_m; % sum of the velocities at the minus end
y3=vg_m-vs_p; % difference between growing vel at the minus end and shortening vel at the plus end

% Main functions
exp1()
%exp2()
%exp3()

% This function shows the relationships between #MTs growing and shortening
% at steady state as a function of the total # of MTs and # of MTs
% shortening at the minus end
    function exp1()
        T=100:50:2500; % total #MTs
        ms_m=10:10:50; % #MTs shortening at the minus end
        
        ct={'b+','bs','bd','b*','bo'};
        for i=1:length(ms_m)
            mg_p=1/y1*(y2*ms_m(i)-y3.*T); % obtains the #MTs growing at the plus end
            ms_p=T-mg_p; % obtains the #MTs shortening at the plus end
            mg_m=T-ms_m(i); % obtains the #MTs growing at the minus end
            
            figure(1);
            hold on;
            plot(T,mg_p./ms_p,ct{i}); % plots the ratio of growing and shortening MTs at the plus end
            figure(2);
            hold on;
            plot(T,mg_m./ms_m(i),ct{i}); % plots the ratio of growing and shortening MTs at the minus end
            
            % I used the lines below just to confirm that my
            % calculations were correct and this is indeed the expression holding
            % the steady state of the system
            
            pol=vg_p.*mg_p+vg_m.*mg_m;
            depol=vs_p.*ms_p+vs_m.*ms_m(i);
            total=pol-depol % This should be very close to zero at the steady state
        end
        
        figure(1);
        set(gca,'fontweight','b','fontsize',16);
        xlabel('T (total #MTs)','fontweight','b','fontsize',16);
        ylabel('MG_+ / MS_+','fontweight','b','fontsize',16);
        legendCell = cellstr(num2str(ms_m', 'MS_- = %-d'));
        legend(legendCell);
        grid on;
        
        figure(2);
        set(gca,'fontweight','b','fontsize',16);
        xlabel('T (total #MTs)','fontweight','b','fontsize',16);
        ylabel('MG_- / MS_-','fontweight','b','fontsize',16);
        legendCell = cellstr(num2str(ms_m', 'MS_- = %-d'));
        grid on;
        legend(legendCell);
    end

% This function shows the relationship between transition frequencies at
% steady state for the plus end of the MT
    function exp2()
        T=2500; % total #MTs
        ms_m=0.1*T; % #MT shortening at the minus end
        pgs_p=(0.001:0.001:0.05)*60; % catastrophe array of values for the plus end
        
        mg_p=1/y1*(y2*ms_m-y3*T); % obtains the #MTs growing at the plus end
        ms_p=T-mg_p; % obtains the #MTs shortening at the plus end
        mg_m=T-ms_m; % obtains the #MTs growing at the minus end
        
        % I used the three lines below to calculate what would be the
        % extinction probability at steady state for the current total # of
        % MTs
        %psg_p=0.028*60;
        %pgs_p=0.013*60;
        %NE=-(T*pgs_p*psg_m*y3 - mg_m*pgs_m*pgs_p*y2 + ms_p*psg_m*psg_p*y1)/(pgs_p*y2 + psg_m*y1)
        %NE/T
        
        Ep=0.06:0.06:0.18; % Extinction probability array of values
        NE=T*Ep; % Obtains the nucleation/extinction value at ss from the extinction probability
        
        figure(3);
        for i=1:length(NE)
            % Obtains the rescue frequency at ss for the plus end of the MT
            psg_p=-(NE(i)*pgs_p*y2 + NE(i)*psg_m*y1 - mg_m*pgs_m*pgs_p*y2 + pgs_p*psg_m*T*y3)/(ms_p*psg_m*y1);
            hold on;
            plot(pgs_p/60,psg_p/60); % plots the rescue frequency (+ end) at ss as a function of the catastrophe frequency (+ end)
        end
        
        hold on;
        plot(0.013,0.028,'r*'); % plots the catrastrophe and rescue frequencies reported in the paper
        set(gca,'fontweight','b','fontsize',16);
        xlabel('Catastrophe (+ MT end)','fontweight','b','fontsize',16);
        ylabel('Rescue (+ MT end)','fontweight','b','fontsize',16);
        grid on;
    end

% This function shows the relationship between transition frequencies at
% steady state for the plus end of the MT and extinction probability
% This is similar to the previous function, but in 3D
    function exp3()
        T=1:50:2500; % total #MTs
        ms_m=0.1*T; % #MT shortening at the minus end
        NE=1:50; % Extinction probability array of values
        pgs_p=(0.001:0.001:0.05)*60; % catastrophe frequencies array of values for the + end of the MT
        
        [X,Y]=meshgrid(NE./T,pgs_p);
        NE_aux=repmat(NE,length(X),1);
        T_aux=repmat(T,length(X),1);
        ms_m_aux=repmat(ms_m,length(X),1);
        
        mg_p=1/y1*(y2.*ms_m_aux-y3.*(NE_aux./X)); % obtains the #MTs growing at the plus end
        ms_p=(NE_aux./X)-mg_p; % obtains the #MTs shortening at the plus end
        mg_m=(NE_aux./X)-ms_m_aux; % obtains the #MTs growing at the minus end
        
        % Obtains the rescue frequency at ss for the plus end of the MT
        psg_p=-(X.*T_aux.*Y*y2 + X.*T_aux.*psg_m*y1 - mg_m.*pgs_m.*Y*y2 + Y.*psg_m.*T_aux*y3)./(ones(length(X),length(X)).*(ms_p.*psg_m*y1));
        
        figure(4);
        surf(X./60,Y./60,psg_p/60); % plots thre relationship between pgs_p, psg_p and NE (extinction probability)
        Z=zeros(length(X),length(X));
        Z(13,:)=0.028;
        hold on;
        surf(X./60,Y./60,Z); % plots the catrastrophe and rescue frequencies reported in the paper
        set(gca,'fontweight','b','fontsize',16);
        xlabel('E/T','fontweight','b','fontsize',16);
        ylabel('Cat','fontweight','b','fontsize',16);
        zlabel('Res');
        view(-70,20);
        colorbar;
    end

% syms pgs_p psg_p pgs_m psg_m NE ms_p mg_m y1 y2 y3 T;
% S=solve('1/pgs_p*(ms_p*psg_p+NE)=1/y1*(y2/psg_m*(mg_m*pgs_m-NE)-y3*T)',psg_p);
% S=-(NE*pgs_p*y2 + NE*psg_m*y1 + T*pgs_p*psg_m*y3 - mg_m*pgs_m*pgs_p*y2)/(ms_p*psg_m*y1)

end














