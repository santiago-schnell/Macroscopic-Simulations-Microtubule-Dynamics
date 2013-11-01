function []=sim_ver3(parts)

% This file contains all the functions needed to interpret the microtubule 
% plant cell cortex simulation results

% Part:

%   1 - Plots MT length, interaction and orientation distribution over time
%   as well as average MT length and number over time

%   2 - Plots free/bound molar tubulin, nucleation/extinction, phase times
%   and run lengths over time

%   3 - Plots velocities over free molar tubulin for all 6 faces and the
%   total average

%   4 - Plots average transitions coming in and going out from all faces of 
%   the cell over time

%   5 - Plots average number of tips growing, shortening and in pause in 
%   all 6 faces of the cell for the plus end over time

%   6 - Displays all MTs as they show up in the last iteration of the
%   simulation

%path='./';
path='Results/08-21-2013/';
it_savs=1; % Results are from iterations at multiples of it_savs
niterations=1000; % Number of iterations
ntimes=10; % Number of time points to plot results
init_times=it_savs:it_savs:niterations; % Vector of result times
fig=4;

if sum(parts==1)>0
    plot_results_part1;
end
if sum(parts==2)>0
    plot_results_part2;
end
if sum(parts==3)>0
    plot_results_part3;
end
if sum(parts==4)>0
    plot_results_part4;
end
if sum(parts==5)>0
    plot_results_part5;
end
if sum(parts==6)>0
    plot_results_part6;
end
if sum(parts==7)>0
    plot_results_part7;
end
if sum(parts==8)>0
    plot_results_part8;
end

% -------------- Code for plotting results-part1 --------------

    function plot_results_part1
        pos=floor(length(init_times)/ntimes); % Find spacing to generate vector of positions
        pos_vec=pos:pos:length(init_times); % Vector of positions in the times vector
        times=init_times(pos_vec); % Final vector of times to generate plots
        
        f1=strcat(path,'Results6-Part1.txt');
        matrix=load(f1);
        orientations=matrix(:,4);
        matrix(:,4)=orientations-180*(floor(orientations/180));
        
        % Defines axis labels and axis scales
        axis_labels={'Polymer length (m)' 'Polymer interactions' 'Polymer orientations'}; % Label for plots
        axis_scale=[min(matrix(:,2)) max(matrix(:,2)) 0 200;0 max(matrix(:,3)) 0 200;0 max(matrix(:,4)) 0 200];
        nbins=[10 10 10]; % Number of bins for the plots
        
        % Generate plots
        for j=1:3
            figure(j);
            for i=1:length(times)
                subplot(4,3,i);
                aux=matrix(find(matrix(:,1)==times(i)),:); % Obtains matrix at relevant time positions
                hist(aux(:,j+1),nbins(j));
                axis(axis_scale(j,:));
                legend(num2str(times(i)));
                xlabel(axis_labels{j},'fontweight','b','fontsize',16);
                ylabel('Polymer number','fontweight','b','fontsize',16);
            end
        end
      
        count=1;
        for i=1:it_savs:niterations
            aux=find(matrix(:,1)==i);
            num(count)=length(aux); % Gets number of MTs at steady state
            len(count)=sum(matrix(aux,2))/length(aux); % Gets average MT length
            count=count+1;
        end
        
        figure(4);
        hold on;
        % Obtain microtubule length over time
        subplot(1,2,1);
        plot(1:it_savs:niterations,len);
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('MT length','fontweight','b','fontsize',16);
        
        % Obtain microtubule number over time
        subplot(1,2,2);
        plot(1:it_savs:niterations,num);
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('MT number','fontweight','b','fontsize',16);
    end

% -------------- End of Code for plotting results-part1 --------------

% -------------- Code for plotting results-part2 --------------

    function plot_results_part2
        f1=strcat(path,'Results6-Part2');
        f2='.txt';
        f3=strcat(f1,f2);
        matrix=load(f3);
        
        % Generate plots
        fig=fig+1;
        figure(fig);
        xaxis=matrix(:,1); % Get iterations
        
        % Plots free and bound molar tubulin over time
        subplot(1,2,1);
        plot(xaxis,matrix(:,2),'color','r');
        hold on;
        plot(xaxis,matrix(:,3),'color','b');
        legend('fmt','bmt');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Molar tubulin ([ ])','fontweight','b','fontsize',16);
        
        % Plots nucleations and extinctions over time
        subplot(1,2,2);
        plot(xaxis,matrix(:,4),'color','r');
        hold on;
        plot(xaxis,matrix(:,5),'color','b');
        legend('Nucleations','Extinctions');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Nuc and Ext','fontweight','b','fontsize',16);
        
        fig=fig+1;
        figure(fig);
        % Plots phase times for the minus end over time
        subplot(2,2,1);
        plot(xaxis,matrix(:,6),'color','r');
        hold on;
        plot(xaxis,matrix(:,7),'color','g');
        hold on;
        plot(xaxis,matrix(:,8),'color','b');
        legend('Ptgm','Ptpm','Ptsm');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Phase time(-end)','fontweight','b','fontsize',16);
        
        % Plots phase times for the plus end over time
        subplot(2,2,2);
        plot(xaxis,matrix(:,9),'color','r');
        hold on;
        plot(xaxis,matrix(:,10),'color','g');
        hold on;
        plot(xaxis,matrix(:,11),'color','b');
        legend('Ptgp','Ptpp','Ptsp');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Phase time(+end)','fontweight','b','fontsize',16);
        
        % Plots run lengths for the minus end over time
        subplot(2,2,3);
        plot(xaxis,matrix(:,12),'color','r');
        hold on;
        plot(xaxis,matrix(:,13),'color','g');
        hold on;
        plot(xaxis,matrix(:,14),'color','b');
        legend('Rlgm','Rlpm','Rlsm');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Run Length(-end)','fontweight','b','fontsize',16);
        
        % Plots run lengths for the plus end over time
        subplot(2,2,4);
        plot(xaxis,matrix(:,15),'color','r');
        hold on;
        plot(xaxis,matrix(:,16),'color','g');
        hold on;
        plot(xaxis,matrix(:,17),'color','b');
        legend('Rlgp','Rlpp','Rlsp');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Run Length(+end)','fontweight','b','fontsize',16);
    end

% -------------- End of Code for plotting results-part2 --------------

% -------------- Code for plotting results-part3 --------------
      
    function plot_results_part3
        f0=strcat(path,'Results6-Part2.txt');
        f1=strcat(path,'Results6-Part3.txt');
        matrix0=load(f0);
        matrix=load(f1);
        [x,~]=size(matrix);
        
        fig=fig+1;
        figure(fig);
        for i=1:1
            plot(matrix0(:,2),matrix(:,4*(i-1)+1).*1e6*60,'color','r');
            hold on;
            plot(matrix0(:,2),matrix(:,4*(i-1)+2).*1e6*60,'color','g');
            hold on;
            plot(matrix0(:,2),matrix(:,4*(i-1)+3).*1e6*60,'color','b');
            hold on;
            plot(matrix0(:,2),matrix(:,4*(i-1)+4).*1e6*60,'color','m');
            legend('POL-','POL+','DEPOL-','DEPOL+');
            xlabel('Free Molar Tubulin [ ]','fontweight','b','fontsize',16);
            aux=strcat('Face ',i,'-Vels(\mum per min)');
            ylabel(aux,'fontweight','b','fontsize',16);
        end
        
%         % Plot total
%         res=zeros(x,4);
%         for j=1:x
%             for k=1:4
%                 res(j,k)=sum(matrix(j,k:4:24))./6;
%             end
%         end
%         fig=fig+1;
%         figure(fig);
%         plot(matrix0(:,2),res(:,1).*1e6*60,'color','r');
%         hold on;
%         plot(matrix0(:,2),res(:,2).*1e6*60,'color','g');
%         hold on;
%         plot(matrix0(:,2),res(:,3).*1e6*60,'color','b');
%         hold on;
%         plot(matrix0(:,2),res(:,4).*1e6*60,'color','m');
%         legend('POL-','POL+','DEPOL-','DEPOL+');
%         xlabel('Free Molar Tubulin [ ]','fontweight','b','fontsize',16);
%         ylabel('Vels(\mum per min)','fontweight','b','fontsize',16);
%         res(x,1)
%         res(x,2)
%         res(x,3)
%         res(x,4)
    end

% -------------- End of Code for plotting results-part3 --------------
        
% -------------- Code for plotting results-part4 --------------

    function plot_results_part4
        f1=strcat(path,'Results6-Part4.txt');
        matrix=load(f1);
        [x,~]=size(matrix);
        in=zeros(x,6);
        out=zeros(x,6);
        for j=1:x
            for i=1:6
                aux1=matrix(j,i:6:36);
                aux2=matrix(j,i*6-5:i*6);
                in(j,i)=sum(aux1);
                out(j,i)=sum(aux2);
            end
        end
        
        fig=fig+1;
        figure(fig);
        subplot(1,2,1);
        plot((1:x)*it_savs,in);
        axis([1 x 0 0.32]);
        grid on;
        legend('Trs. to 1','Trs. to 2','Trs. to 3','Trs. to 4','Trs. to 5','Trs. to 6');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Avg. # of transitions coming in per MT','fontweight','b','fontsize',16);
        
        subplot(1,2,2);
        plot((1:x)*it_savs,out);
        axis([1 x 0 0.32]);
        grid on;
        legend('Trs. from 1','Trs. from 2','Trs. from 3','Trs. from 4','Trs. from 5','Trs. from 6');
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Avg. # of transitions going out per MT','fontweight','b','fontsize',16);
        
%         %aux=mean(out(niterations-100:niterations,:))./mean(in(niterations-100:niterations,:));
%         aux=mean(out(niterations-10:niterations,:))./mean(in(niterations-10:niterations,:));
%         fig=fig+1;
%         figure(fig);
%         hold on;
%         plot(1:6,aux,'b*','MarkerSize',16);
%         hold on;
%         plot(1:6,ones(1,6),'k--');
%         grid on;
%         xlabel('Face','fontweight','b','fontsize',16);
%         ylabel('MTs out / MTs in','fontweight','b','fontsize',16);
%         legend('Homogeneous nucleation','Heterogeneous nucleation');
    end

% -------------- End of Code for plotting results-part4 --------------

% -------------- Code for plotting results-part5 --------------
      
    function plot_results_part5
        f1=strcat(path,'Results6-Part5.txt');
        matrix=load(f1);
        [x,~]=size(matrix);
        avgs=zeros(6,3,2);
        for j=1:2
            fig=fig+1;
            figure(fig);
            for i=1:6
                subplot(2,3,i);
                plot((1:x/2)*it_savs,matrix(j:2:x,3*(i-1)+1),'color','r');
                %avgs(i,1,j)=mean(matrix(x-2000+j:2:x,3*(i-1)+1));
                avgs(i,1,j)=mean(matrix(x-10+j:2:x,3*(i-1)+1));
                hold on;
                plot((1:x/2)*it_savs,matrix(j:2:x,3*(i-1)+2),'color','g');
                avgs(i,2,j)=mean(matrix(x-10+j:2:x,3*(i-1)+2));
                hold on;
                plot((1:x/2)*it_savs,matrix(j:2:x,3*(i-1)+3),'color','b');
                avgs(i,3,j)=mean(matrix(x-10+j:2:x,3*(i-1)+3));
                axis([1 x/2 0 0.25]);
                legend('Growth','Pause','Short');
                xlabel('Time (Sec)','fontweight','b','fontsize',16);
                aux=strcat('Face ',num2str(i),'-Avg.# Tips per area');
                ylabel(aux,'fontweight','b','fontsize',16);
                if j==1
                    title('Data for the minus end','fontweight','b','fontsize',16);
                else
                    title('Data for the plus end','fontweight','b','fontsize',16);
                end
            end
        end
        
        areas=[15^2 15^2 15*40 15*40 15*40 15*40];
        nmts=zeros(6,3,2);
        for j=1:2
           for i=1:6
               nmts(i,:,j)=avgs(i,:,j).*areas(i);
           end
        end
        
        nmts(:,:,1)
        nmts(:,:,2)
        system_pol=(nmts(:,1,1).*1.0109e-8)+(nmts(:,1,2).*1.0414e-7) % polymerization
        %system_pol=(nmts(:,1,1).*1.3948e-8)+(nmts(:,1,2).*1.1310e-7) % polymerization
        system_depol=(nmts(:,3,1).*1.6667e-8)+(nmts(:,3,2).*3e-7) % depolymerization
        system_ss_diff=system_pol-system_depol
        system_ss_ratio=system_pol./system_depol
        system_ss_tt=sum(system_ss_diff)
        
        totals(:,1:3)=nmts(:,:,1);
        totals(:,4:6)=nmts(:,:,2);
        totals(:,7)=system_pol;
        totals(:,8)=system_depol;
        save('homo_data.mat','totals');
    end

% -------------- End of Code for plotting results-part5 --------------

% -------------- Code for plotting results-part6 --------------

    function plot_results_part6
        
        % Creates the faces of the cube... each face has an identifier, a position,
        % a dimension, info about transition faces and transition sides
        face1=Face(1,[-15 0],[15 15],[5 3 6 4],[0 0 1 1 1 1]);
        face2=Face(2,[40 0],[15 15],[3 5 6 4],[0 0 2 2 2 2]);
        face3=Face(3,[0 0],[40 15],[1 2 6 4],[2 1 0 3 0 4]);
        face4=Face(4,[0 15],[40 15],[1 2 3 5],[4 4 4 0 3 0]);
        face5=Face(5,[0 30],[40 15],[1 2 4 6],[1 2 0 4 0 3]);
        face6=Face(6,[0 45],[40 15],[1 2 5 3],[3 3 3 0 4 0]);
        faces={face1 face2 face3 face4 face5 face6}; % Stores all faces
        
        f1=strcat(path,'Results2-Part6');
        f2='.txt';
        f3=strcat(f1,f2);
        fid = fopen(f3);
        
        % Creates the MTs from the file
        MTs={};
        tline = fgetl(fid);
        %aux_tt=zeros(1,100);
        count=1;
        %print=0;
        while ischar(tline)
            line=str2num(tline); % Obtains the array of numbers from the line
            for i=1:(length(line)/4)
                frontier=line(4*i-3); % Obtains the frontier point
                pos=Posit([line(4*i-1) line(4*i)],0,0,frontier,0,0); % Creates the position
                aux=pos.get_position;
                x=aux(1);
                y=aux(2);
                if i==1
                    mt=MT({pos},0,0,0); % Creates the polymer with the initial position
                    %if (sum(aux==[0.00000516 0.00002980])==2)
                    %    print=1;
                    %end
                else
                    %last_pos=mt.get_plus_end_pos;
                    %if ~(pos.get_frontier_point && last_pos.get_frontier_point)
                    %    aux_tt(count)=aux_tt(count)+norm(aux-last_pos.get_position);
                    %end
                    mt.add_position(pos);
                end
            end
            %if (print)
                MTs{end+1}=mt;
                %print=0;
            %end
            tline = fgetl(fid);
            count=count+1;
        end
        
        fclose(fid);
        create_frame(faces); % Create rectangular boxes
        for i=1:length(MTs)
            MTs{i}.display_MT(); % Display the MT
        end
        axis([-15 60 0 60]*1e-6);
    end

% This function creates the rectangular box
    function [h]=create_frame(faces)
        h=figure;
        for i=1:6
            rect_faces=[faces{i}.get_position faces{i}.get_dimensions];
            rectangle('Position',rect_faces);
            hold on;
        end
    end

% -------------- End of Code for plotting results-part6 --------------

        
% -------------- Code for plotting results-part7 --------------

    function plot_results_part7
%         f1=strcat(path,'Results6-Part8.txt');
%         matrix=load(f1);
%         [x,~]=size(matrix);
%         fig=fig+1;
%         figure(fig);
%         hold on;
%         aux=diff(matrix(1:x,11));
%         plot(2:x,aux,'r');
%         matrix(find(aux>40e-6),:)
%         matrix(find(aux>40e-6)+1,:)
%         matrix(find(aux>40e-6)+2,:)
        
        f1=strcat(path,'Results3-Part7.txt');
        matrix=load(f1);
        [x,~]=size(matrix);
        fig=fig+1;
        figure(fig);
        hold on;
        plot(matrix(1:x,1),matrix(1:x,2),'g');
        hold on;
        plot(matrix(1:x,1),matrix(1:x,3),'r');
        set(gca,'FontWeight','bold','FontSize',16);
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Total polymer in the region','fontweight','b','fontsize',16);
        fig=fig+1;
        figure(fig);
        %format long
        %aux=[matrix(2:x,1) diff(matrix(1:x,3))]
        plot(matrix(2:x,1),diff(matrix(1:x,3)));
        axis([1 x -10e-5 10e-5]);
        set(gca,'FontWeight','bold','FontSize',16);
        xlabel('Time (sec)','fontweight','b','fontsize',16);
        ylabel('Derivative','fontweight','b','fontsize',16);
        
        aux=diff(matrix(1:x,3));
        aux(end)
        
    end

% -------------- End of Code for plotting results-part7 --------------

% -------------- Code for plotting results-part8 --------------

    function plot_results_part8
        mat1=load('homo_data.mat');
        mat2=load('hetero_data.mat');
        
%         fig=fig+1;
%         figure(fig);
%         plot(1:6,mat1.totals(:,1)./mat2.totals(:,1),'r.',1:6,mat1.totals(:,2)./mat2.totals(:,2),'g*',1:6,mat1.totals(:,3)./mat2.totals(:,3),'bo');
%         xlabel('Face','fontweight','b','fontsize',16);
%         ylabel('Homo / Hetero number (minus end)','fontweight','b','fontsize',16);
%         hold on;
%         plot(1:6,ones(1,6),'k--');
%         axis([1 6 0.6 1.8]);
%         grid on;
%         legend('Growing','Pause','Shortening');
%         
%         fig=fig+1;
%         figure(fig);
%         plot(1:6,mat1.totals(:,4)./mat2.totals(:,4),'r.',1:6,mat1.totals(:,5)./mat2.totals(:,5),'g*',1:6,mat1.totals(:,6)./mat2.totals(:,6),'bo');
%         hold on;
%         plot(1:6,ones(1,6),'k--');
%         xlabel('Face','fontweight','b','fontsize',16);
%         ylabel('Homo / Hetero number (plus end)','fontweight','b','fontsize',16);
%         axis([1 6 0.6 1.8]);
%         grid on;
%         legend('Growing','Pause','Shortening');
        
        fig=fig+1;
        figure(fig);
        plot(1:6,mat1.totals(:,1)./mat1.totals(:,3),'r.',1:6,mat2.totals(:,1)./mat2.totals(:,3),'b*');
        hold on;
        plot(1:6,ones(1,6),'k--');
        xlabel('Face','fontweight','b','fontsize',16);
        ylabel('Growing / Shortening number (minus end)','fontweight','b','fontsize',16);
        axis([1 6 0.03 0.06]);
        grid on;
        legend('Homogeneous nucleation','Heterogeneous nucleation');

        fig=fig+1;
        figure(fig);
        plot(1:6,mat1.totals(:,4)./mat1.totals(:,6),'r.',1:6,mat2.totals(:,4)./mat2.totals(:,6),'b*');
        hold on;
        plot(1:6,ones(1,6),'k--');
        xlabel('Face','fontweight','b','fontsize',16);
        ylabel('Growing / Shortening number (plus end)','fontweight','b','fontsize',16);
        axis([1 6 3 3.8]);
        grid on;
        legend('Homogeneous nucleation','Heterogeneous nucleation');
        
        fig=fig+1;
        figure(fig);
        plot(1:6,mat1.totals(:,7)./mat1.totals(:,8),'r.',1:6,mat2.totals(:,7)./mat2.totals(:,8),'b*','MarkerSize',16);
        hold on;
        plot(1:6,ones(1,6),'k--');
        grid on;
        xlabel('Face','fontweight','b','fontsize',16);
        ylabel('Net pol / Net depol','fontweight','b','fontsize',16);
        legend('Homogeneous nucleation','Heterogeneous nucleation');
        ss_value_homo=sum(mat1.totals(:,7))-sum(mat1.totals(:,8))
        ss_value_hetero=sum(mat2.totals(:,7))-sum(mat2.totals(:,8))
        
    end

% -------------- End of Code for plotting results-part8 --------------

end





