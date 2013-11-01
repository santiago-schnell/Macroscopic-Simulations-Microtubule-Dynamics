function []=stt_ver5(NIterations,NMTs,NSimulations)

% Plus and minus end tip state
% 0 - Growth
% 1 - Pause
% 2 - Short

% Velocities of polymerization and depolymerization in mic per sec
vg1=7.33/60;
vg2=1.07/60;
vs1=15/60;
vs2=5/60;

% vg1=7.33/60;
% vg2=1.07/60;
% vs1=18.98/60;
% vs2=1.07/60;

% Transition frequencies per sec for the plus end
tfp=[0.986 0.001 0.013;0.055 0.883 0.062;0.027 0.003 0.969];
% Transition frequencies per sec for the minus end
tfm=[0.961 0.014 0.025;0.003 0.965 0.032;0.001 0.013 0.985];

MTs=[]; % Stores all MTs
MTs=create_MTs(str2num(NMTs),MTs); % Creates initial MTs

% Main cycle
for tt=1:str2num(NSimulations)
    file=strcat('lifetime_',num2str(tt),'.txt');
    fid=fopen(file,'w');
    new_MTs=MTs;
    ext_times=[]; % Stores MT extinction times
    for i=1:str2num(NIterations)
        new_MTs=create_MTs(1,new_MTs);
        rem=[];
        phase_time_plus=[];
        phase_time_minus=[];
        run_length_plus=[];
        run_length_minus=[];
        for j=1:length(new_MTs)
            new_MTs(j).state_hist_plus=[new_MTs(j).state_hist_plus new_MTs(j).plus_end_tip_state];
            new_MTs(j).state_hist_minus=[new_MTs(j).state_hist_minus new_MTs(j).minus_end_tip_state];
            new_MTs(j)=update_state_MT(new_MTs(j)); % updates the state of both tips of a MT
            new_MTs(j)=update_length_MT(new_MTs(j)); % update the length of a MT
            for k=1:3
                phase_time_plus(j,k)=sum(new_MTs(j).state_hist_plus==k)/length(new_MTs(j).state_hist_plus);
                phase_time_minus(j,k)=sum(new_MTs(j).state_hist_minus==k)/length(new_MTs(j).state_hist_minus);
                run_length_plus(j,k)=1/(1-phase_time_plus(j,k));
                run_length_minus(j,k)=1/(1-phase_time_minus(j,k));
            end
            if (new_MTs(j).length)<=0
                ext_times=[ext_times i];
                rem=[rem j];
            end
        end
        
        res1(i,:)=mean(phase_time_plus);
        res2(i,:)=mean(phase_time_minus);
        res3(i,:)=mean(run_length_plus);
        res4(i,:)=mean(run_length_minus);
        
        new_MTs(rem)=[]; % Remove MTs with length < 0
        nEXT(i)=length(rem);
        nMTs(i)=length(new_MTs);
    end
    
    size(res1)
    
    figure(1)
    plot(1:str2num(NIterations),res1)
    legend('G','P','S');
    figure(2)
    plot(1:str2num(NIterations),res2)
    legend('G','P','S');
    figure(3)
    plot(1:str2num(NIterations),res3)
    legend('G','P','S');
    figure(4)
    plot(1:str2num(NIterations),res4)
    legend('G','P','S');
    figure(5);
    plot(1:str2num(NIterations),nMTs);
    figure(6);
    plot(1:str2num(NIterations),nEXT);
    
    fprintf(fid,'%g\n',ext_times);
    fclose(fid);
end

    % This function creates and initializes a MT
    function [MTs]=create_MTs(NMTs,MTs)
        for xxx=1:NMTs
            MT.length=vg1; % MTs start growing
            MT.plus_end_tip_state=1; % The MT plus end state is initially on growth
            MT.minus_end_tip_state=2; % The MT minus end state is initially on pause
            MT.state_hist_plus=[];
            MT.state_hist_minus=[];
            MTs=[MTs MT]; % Adds new MT to the array of MTs
        end
    end

    % This function updates the length of the MT
    function [MT]=update_length_MT(MT)
        % Update the MT plus end 
        if MT.plus_end_tip_state==1
            MT.length=MT.length+vg1;
        else if MT.plus_end_tip_state==3
                MT.length=MT.length-vs1;
            end
        end
        % Update the MT minus end 
        if MT.minus_end_tip_state==1
            MT.length=MT.length+vg2;
        else if MT.plus_end_tip_state==3
                MT.length=MT.length-vs2;
            end
        end
    end
    
    % This function updates the state of the MT
    function [MT]=update_state_MT(MT)
        tm_aux1={'tfp' 'tfm'}; % Sets the two possible transition frequencies matrices
        tm_aux2={'plus_end_tip_state' 'minus_end_tip_state'}; % Sets the two possible end state fields
        for k=1:2 % Update both ends of the MT
            val=rand(1);
            aux_comp=1:3;
            ets=MT.(tm_aux2{k}); % Get the end tip state of the MT
            aux_comp(ets)=[]; % Sets possible comparisons for the transition
            tf1=eval(strcat(tm_aux1{k},'(ets,aux_comp(1))')); 
            if val<=tf1 % Check first transition probability
                MT.(tm_aux2{k})=aux_comp(1); % Change state of the tip to first possible transition
            else
                if val<=(tf1+eval(strcat(tm_aux1{k},'(ets,aux_comp(2))'))) % Check first transition probability
                    MT.(tm_aux2{k})=aux_comp(2); % Change state of the tip to second possible transition
                end
            end
        end
    end
 
end