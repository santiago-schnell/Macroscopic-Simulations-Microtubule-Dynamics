function []=sim_ver2(niterations,nsites,prob1,prob2,prob3,prob4,prob5,folder)

probs=[prob1 prob2 prob3 prob4 prob5]; % Orientation probabilities on MT dependent nucleation
%nfprobs=[0.1 0.1 0.1 0.1 0.5 0.1]; % Probabilities of nucleation in cell faces (1-6)
nfprobs=[0.1667 0.1667 0.1667 0.1667 0.1667 0.1667]; % Probabilities of nucleation in cell faces (1-6)
%nfprobs=[0 0 0 0 1 0]; % Probabilities of nucleation in cell faces (1-6)
%sfprobs=[1 1 1 1 1 1]; % Probabilities of severing in cell faces (1-6)
frap_time=1; % Indicator time of when the FRAP simulation starts

%Probs: 1- Probability of nucleation on an existent MT
%       2- Probability of a new MT to orient towards the right
%       3- Probability of a new MT to orient towards the left
%       4- Probability of a new MT to orient towards the MT+ end
%       5- Probability of a new MT to orient towards the MT- end

clc;
close all;

% ---------------------------- Tasks ---------------------------------
graphics=1; % Enables graphics

% --------------------------- Parameters -----------------------------

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

% Constants
ndimer=1634;              % Dimers per micrometer of microtubule
avg=6.022e23;             % Avogadro's number
volume=3000;              % Volume of cell in cubic micrometers
initial_tubulin=10;       % Initial tubulin concentration in micromolar
orientations=0:10:350;    % Possible polymer orientations

% Linear models for the transition frequencies of both minus and plus end
tfm={'' '-k21*fmt+k22' '-k23*fmt+k24'; 'k25*fmt+k26' '' '-k27*fmt+k28'; 'k29*fmt+k30' 'k31*fmt+k32' ''};
tfp={'' '-k9*fmt+k10' '-k11*fmt+k12'; 'k13*fmt+k14' '' '-k15*fmt+k16'; 'k17*fmt+k18' 'k19*fmt+k20' ''};

% Creates the faces of the cube... each face has an identifier, a position,
% a dimension, info about transition faces and transition sides
face1=Face(1,[-15 0],[15 15],[5 3 6 4],[0 0 1 1 1 1]);
face2=Face(2,[40 0],[15 15],[3 5 6 4],[0 0 2 2 2 2]);
face3=Face(3,[0 0],[40 15],[1 2 6 4],[2 1 0 3 0 4]);
face4=Face(4,[0 15],[40 15],[1 2 3 5],[4 4 4 0 3 0]);
face5=Face(5,[0 30],[40 15],[1 2 4 6],[1 2 0 4 0 3]);
face6=Face(6,[0 45],[40 15],[1 2 5 3],[3 3 3 0 4 0]);
%face7=Face(7,[15 35],[6 4],[],[]); % Defines FRAP region
face7=Face(7,[0 30],[40 15],[],[]);
faces={face1 face2 face3 face4 face5 face6 face7}; % Stores all faces

% Initiate graphics if needed
if ~graphics
    gh=0;
else
    gh=create_frame(); % Create rectangular boxes
end

% --------------------------- Main Algorithm --------------------------

it_savs=1; % simulation returns results from iterations at multiples of it_savs
branch_id=0; % Stores the branch ID

fid1=fopen(strcat(folder,'/Results-Part1.txt'),'w+');
fid2=fopen(strcat(folder,'/Results-Part2.txt'),'w+');

free_sites=nsites;
fmt=initial_tubulin;
initial_free_subunits=fmt*1e-6*avg*volume*1e-18/0.001;
free_subunits=initial_free_subunits;
tt_pol_aux=zeros(1,4);
problem=0;

tic
MTs={}; % Stores all the MTs in the system
for iteration=1:niterations % 1 iteration simulates 1 second of polymer dynamics
    display(sprintf('Iteration: %d\n',iteration))
    ft=free_subunits/ndimer*1e-6; % Converts free subunits to meters
    vg_m=(k6*ft-k7)/60; % Determines the polymerization velocity at the minus end (in secs)
    vs_m=k8/60; % Determines the depolymerization velocity at the minus end (in secs)
    vg_p=(k3*ft-k4)/60; % Determines the polymerization velocity at the plus end (in secs)
    vs_p=k5/60; % Determines the depolymerization velocity at the plus end (in secs)    
    prob_nuc=k1*fmt-k2; % Determines probability of nucleation (in sec)
    N=round(prob_nuc*free_sites); % Determines number of new polymers (in sec)
%     if (iteration==1)
%         N=1;
%     else
%         N=0;
%     end
    create_MTs(N,vg_p,vs_p,iteration,probs,nfprobs); % Creates N MTs
    free_sites=free_sites-N; % Updates the number of new free sites
    trans_m=evaluate(tfm); % Determines transition frequencies for the minus end
    trans_p=evaluate(tfp); % Determines transition frequencies for the plus end
    
    remov=[];
    newMTs={};
    free_sites_iter=0;
    parfor k=1:length(MTs) % Updates all MTs
        rel_nuc_site1=0;
        rel_nuc_site2=0;
        MTs{k}=MTs{k}.update_state(trans_m,trans_p); % updates the state of both MT tips
        MTs{k}=MTs{k}.init_mt_length_vars; % Initializes MT length variations
        [MTs{k},extinct,rel_nuc_site1]=MTs{k}.update_position(1,[vg_m vs_m],faces); % updates the MT minus end position
        if ~extinct
            [MTs{k},extinct,rel_nuc_site2]=MTs{k}.update_position(2,[vg_p vs_p],faces); % updates the MT plus end position
        end
        if (rel_nuc_site1 || rel_nuc_site2)
            free_sites_iter=free_sites_iter+1;
        end
        if extinct
            remov=[remov k];
%         else
%             if ((iteration==5) && (length(MTs{k}.get_positions)>4)) % I am requiring that the MT needs to have at least 5 points in order to be severed
%                 [MTs{k},newMT]=MTs{k}.severe(sfprobs,gh); % Severes the MT with a probability
%                 if ~(newMT==-1)
%                     newMTs=[newMTs {newMT}];
%                 end
%             end
        end
    end
    
    % Removes extinct MTs from the display
    if gh
        for k=1:length(remov)
            MTs{remov(k)}.display_MT();
        end
    end
    MTs(remov)=[]; % Removes extinct MT from the list of MTs
    MTs=[MTs newMTs]; % Add newly formed MTs from severing to the list
    free_sites=free_sites+free_sites_iter;
    
    length_MTs=0;
    MTs_vars=zeros(1,4);
    phase_times=zeros(2,3);
    run_lengths=zeros(2,3);
    length_MTs_vars=zeros(1,4);
    face_transitions=zeros(6,6);
    npolends_minus=zeros(6,3);
    npolends_plus=zeros(6,3);
    arrayfun(@(x)faces{x}.set_MTs({}),1:6); % Clears all MTs within the faces of the cube
    for k=1:length(MTs) % Updates the MTs within the faces of the cube
        face_ID=MTs{k}.get_plus_end_pos.get_face_ID;
        state_MT=MTs{k}.get_tip_state;
        MTs_vars=MTs_vars+(MTs{k}.get_mt_length_vars>0);
        length_MTs_vars=length_MTs_vars+MTs{k}.get_mt_length_vars;
        phase_time_MT=MTs{k}.get_state_hist/(1+iteration-MTs{k}.get_init_it);
        run_length_MT=MTs{k}.get_run_length;
        phase_times=phase_times+phase_time_MT;
        run_lengths=run_lengths+run_length_MT;
        face_transitions=face_transitions+(MTs{k}.get_face_transitions);
        arrayfun(@(x)faces{x}.add_MT(MTs{k}),unique(MTs{k}.get_faces_ID));
        length_MTs=length_MTs+MTs{k}.get_mt_length;
        npolends_minus(face_ID,state_MT(1))=npolends_minus(face_ID,state_MT(1))+1;
        npolends_plus(face_ID,state_MT(2))=npolends_plus(face_ID,state_MT(2))+1;
        if gh
            MTs{k}.display_MT(); % Displays the MT if graphics are enabled
        end
    end
    if gh
        drawnow; % Flushes the graphic handles if graphics enabled
    end
    
    % -------------------- Update tubulin concentration ------------------
    bound_subunits=length_MTs*ndimer/1e-6; % Converts length in meters to bound subunits
    free_subunits=initial_free_subunits-bound_subunits; % Obtains free subunits
    fmt=free_subunits/(1e-6*avg*volume*1e-18/0.001); % Determines free molar concentration
    bmt=bound_subunits/(1e-6*avg*volume*1e-18/0.001); % Determines bound molar concentration
    
    % -------------------- Collect results --------------------
    if rem(iteration,it_savs)==0
        write_results(iteration); % writes to file
    end
    
    % -------------------- FRAP simulation --------------------
    
    if iteration>=frap_time
        if iteration==frap_time
            fid7=fopen(strcat(folder,'/Results-Part7.txt'),'w');
        else if iteration==(frap_time+599)
                bleach_region();
             end
        end
        tt_pol=get_polymer_within_region(0,0,0);
        tt_pol_aux(1)=tt_pol_aux(3);
        tt_pol_aux(2)=tt_pol_aux(4);
        tt_pol_aux(3:4)=tt_pol;
        
        fprintf(fid7,'%6.8f %6.8f %6.8f\n',iteration,tt_pol); % iteration,total,bleached
        
        % Writes MT info to file for debug purposes
        %if (iteration==5)
        if (abs(tt_pol_aux(1)-tt_pol_aux(3))>25e-6)
            write_debug_info(iteration,tt_pol_aux);
            problem=1;
            break;
        end
    end
    
end

if ~problem
    write_MT_positions(MTs);
end

fclose(fid1);
fclose(fid2);
if iteration>=frap_time
    fclose(fid7);
end

toc

% ------------- This section sets auxiliary functions --------------

% This function writes MT info to file for debug purposes
    function write_debug_info(i,tt_pol_aux)
        face_MTs=faces{5}.get_MTs();
        write_MT_positions(face_MTs); % Writes all face 5 MT positions to file so we are able to visualize the MTs
        get_polymer_within_region(1,i,tt_pol_aux);
    end
       
% This function tells if a position is within the FRAP region
    function res=in_region(new_pos)
        res=0;
        pos_face=faces{7}.get_position;
        dim_face=faces{7}.get_dimensions;
        if (new_pos(1)>=pos_face(1) && new_pos(1)<=(pos_face(1)+dim_face(1)) && ...
                new_pos(2)>=pos_face(2) && new_pos(2)<=(pos_face(2)+dim_face(2)))
            res=1;
        end
    end

% This function modifies positions of a MT between ind_pos1 and ind_pos2
    function bleaches_MT(ind_pos1,ind_pos2,mt_pos)
        for n=ind_pos1:ind_pos2
            mt_pos{n}.set_bleached(1);
        end
    end

% This function sums segments of a MT between ind_pos1 and ind_pos2
    function res=sum_pol_MT(int_pos1,ind_pos1,ind_pos2,mt_pos)
        res=zeros(1,2);
        for n=ind_pos1:ind_pos2-1
            dist=norm(mt_pos{n+1}.get_position-mt_pos{n}.get_position);
            res(1)=res(1)+dist;
            if (~(mt_pos{n+1}.get_bleached) || ~(mt_pos{n}.get_bleached))
                res(2)=res(2)+dist;
            end
        end
        if (length(int_pos1)>1)
            dist=norm(mt_pos{ind_pos1}.get_position-int_pos1);
            res(1)=res(1)+dist;
            if ~(mt_pos{ind_pos1}.get_bleached)
                res(2)=res(2)+dist;
            end
        end
    end

% This function obtains the intersections of the MT with the FRAP region
    function output=find_intersections(mt_pos)
        mat_pos=ones(length(mt_pos)-1,4).*(-10);
        for m=2:length(mt_pos)
            pos1_obj=mt_pos{m-1};
            pos2_obj=mt_pos{m};
            if ~(pos1_obj.get_frontier_point && pos2_obj.get_frontier_point)
                pos1=pos1_obj.get_position;
                pos2=pos2_obj.get_position;
                mat_pos(m-1,:)=[pos1(1) pos1(2) pos2(1) pos2(2)]; % Generate matrix of positions for intersection
            end
        end
        [a,b,c,d]=faces{7}.get_rect_points; % Gets the four points that define the face of the cube
        output=lineSegmentIntersect(mat_pos,[a b;d c;a d;b c]); % Check intersection
    end

% This function bleaches positions within the region (face 7)
    function bleach_region()
        face_MTs=faces{5}.get_MTs();
        for l=1:length(face_MTs)
            mt_pos=face_MTs{l}.get_positions;
            output=find_intersections(mt_pos);
            iam=output.intAdjacencyMatrix; % Obtains intersection results
            [x,y]=find(iam'==1); % y vector corresponds to intersection segment indexes
            [~,index1,~]=unique([output.intMatrixX(y,x) output.intMatrixY(y,x)],'rows','stable');
            y=y(index1);
            in=in_region(mt_pos{1}.get_position);
            ind_pos1=in;
            ind_pos2=0;
            for m=1:length(y)
                in=1-in; % Toggles
                if in
                    ind_pos1=y(m);
                else
                    ind_pos2=y(m)+1;
                    bleaches_MT(ind_pos1,ind_pos2,mt_pos);
                end
            end
            if in
                ind_pos2=length(mt_pos);
                bleaches_MT(ind_pos1,ind_pos2,mt_pos);
            end
        end
    end

% This function calculates the amount of new polymer within the region (face 7)
    function total_pol=get_polymer_within_region(debug,i,tt_pol_aux)
        if (debug)
            fid8=fopen(strcat(folder,'/Results-Part8.txt'),'w'); % Open results_part8
        end
        
        face_MTs=faces{5}.get_MTs();
        total_pol=zeros(1,2);
        for l=1:length(face_MTs)
            mt_pos=face_MTs{l}.get_positions;
            output=find_intersections(mt_pos);
            iam=output.intAdjacencyMatrix; % Obtains intersection results
            [x,y]=find(iam'==1); % y vector corresponds to intersection segment indexes
            [~,index1,~]=unique([output.intMatrixX(y,x) output.intMatrixY(y,x)],'rows','stable');
            x=x(index1);
            y=y(index1);
            in=in_region(mt_pos{1}.get_position);
            int_pos1=-1;
            int_pos2=-1;
            ind_pos1=in;
            ind_pos2=0;
            pol=zeros(1,2);
            for m=1:length(y)
                in=1-in; % Toggles
                if in
                    int_pos1=[output.intMatrixX(y(m),x(m)) output.intMatrixY(y(m),x(m))];
                    ind_pos1=y(m); % Intersection index point
                else
                    ind_pos2=y(m)+1; % Intersection index point
                    pol=pol+sum_pol_MT(int_pos1,ind_pos1,ind_pos2,mt_pos);
                end
            end
            if in
                ind_pos2=length(mt_pos);
                pol=pol+sum_pol_MT(int_pos1,ind_pos1,ind_pos2,mt_pos);
            end
            total_pol=total_pol+pol;
            
            % Debugging...
            if (debug)
                if (isempty(y))
                    fprintf(fid8,'%i %6.8f %6.8f %6.8f %6.8f %i %i %i %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %i %i %i %i \n',i,tt_pol_aux,l,length(mt_pos),length(iam),mt_pos{1}.get_position,total_pol,face_MTs{l}.get_mt_length,pol,0,0,0,0);
                end
                for m=1:length(y)
                    fprintf(fid8,'%i %6.8f %6.8f %6.8f %6.8f %i %i %i %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %i %i %6.8f %6.8f \n',i,tt_pol_aux,l,length(mt_pos),length(iam),mt_pos{1}.get_position,total_pol,face_MTs{l}.get_mt_length,pol,x(m),y(m),output.intMatrixX(y(m),x(m)),output.intMatrixY(y(m),x(m)));
                end
            end
            
        end
        
        if (debug)
            fclose(fid8);
        end
        
    end

% This function writes current results to output files
    function write_results(i)
        % Write results_part1
        interm_length=zeros(1,length(MTs));
        interm_interactions=zeros(1,length(MTs));
        interm_orientations=zeros(1,length(MTs));
        for l=1:length(MTs)
            interm_length(l)=MTs{l}.get_mt_length;
            interm_interactions(l)=MTs{l}.get_interactions;
            interm_orientations(l)=MTs{l}.get_plus_end_pos.get_orientation(1);
        end
        y1=[ones(1,l)*i;interm_length;interm_interactions;interm_orientations];
        fprintf(fid1,'%6.8f %6.8f %6.8f %6.8f\n',y1);
        
        % Write results_part2
        avg_pts=phase_times./length(MTs);
        avg_rls=run_lengths./length(MTs);
        y2=[i;fmt;bmt;N;length(remov);[avg_pts(1,:) avg_pts(2,:)]';[avg_rls(1,:) avg_rls(2,:)]'];
        fprintf(fid2,'%6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f\n',y2);
        
        % Write results_part3
        clear vels;
        vels=length_MTs_vars./MTs_vars;
        dlmwrite(strcat(folder,'/Results-Part3.txt'),vels(1,:),'-append');
        
        % Write results_part4
        clear trans;
        trans=face_transitions./length(MTs);
        dlmwrite(strcat(folder,'/Results-Part4.txt'),[trans(1,:) trans(2,:) trans(3,:) trans(4,:) trans(5,:) trans(6,:)],'-append');
        
        % Write results_part5
        areas=zeros(1,6);
        for l=1:6
            dim=(faces{l}.get_dimensions).*1e6;
            areas(l)=dim(1)*dim(2);
        end
        clear nmns;
        clear npls;
        nmns=npolends_minus./(repmat(areas',1,3));
        npls=npolends_plus./(repmat(areas',1,3));
        dlmwrite(strcat(folder,'/Results-Part5.txt'),[nmns(1,:) nmns(2,:) nmns(3,:) nmns(4,:) nmns(5,:) nmns(6,:)],'-append');
        dlmwrite(strcat(folder,'/Results-Part5.txt'),[npls(1,:) npls(2,:) npls(3,:) npls(4,:) npls(5,:) npls(6,:)],'-append');
    end

% This function writes the current MT positions to an output file
    function write_MT_positions(new_MTs)
        % Write results_part5
        fid6=fopen(strcat(folder,'/Results-Part6.txt'),'w');
        for i=1:length(new_MTs)
            positions=new_MTs{i}.get_positions;
            for j=1:length(positions)
                fp=positions{j}.get_frontier_point; % frontier point
                bp=positions{j}.get_bleached; % bleached
                position=positions{j}.get_position;
                fprintf(fid6,'%1.0f %1.0f %6.8f %6.8f ',fp,bp,position);
            end
            fprintf(fid6,'\n');
        end
        fclose(fid6);
    end

% This function produces a matrix of transition values from the cell matrix of transitions
    function res=evaluate(transitions)
        [x,y]=size(transitions);
        res=zeros(x,y);
        for m=1:x
            for n=1:y
                if ~(m==n)
                    res(m,n)=eval(transitions{m,n});
                end
            end
        end
    end

% This function creates the rectangular box
    function [h]=create_frame()
        h=figure;
        for i=1:6
            rect_faces=[faces{i}.get_position faces{i}.get_dimensions];
            rectangle('Position',rect_faces);
            hold on;
        end
    end

% This function creates N microtubules
    function create_MTs(N,vg_p,vs_p,it,probs,nfprobs)
        for l=1:N
            faces_pos=find(rand<cumsum(nfprobs));
            new_face_ID=faces_pos(1); % Selects the polymer's face
            face=faces{new_face_ID}; % Gets the face object
            face_MTs=face.get_MTs(); % Gets the list of MTs from the face
            bound=0;
            if ~isempty(face_MTs) && rand<=probs(1)
                mt1=face_MTs{randi(length(face_MTs))}; % Selects a MT for the microtubule dependent nucleation
                posits=mt1.get_positions; % Gets all the points of the MT
                posits=posits(cellfun(@(x)((x.get_nuc_point) && (~x.get_frontier_point) && (x.get_face_ID==new_face_ID)),posits)>0); % Select positions capable of nucleation
                if ~isempty(posits)
                    bound=1;
                    random_pos=randi(length(posits)); % Selects the position
                    branch_id=branch_id+1; % Sets new branch_id
                    pos1=posits{random_pos}; % Gets a position in the MT
                    pos1.set_nuc_point(0); % Sets nucleation point
                    pos1.set_branch_id(branch_id); % Sets branch point
                    minus_end_pos=pos1.get_position; % Selects the minus end position
                    % Gets the segment end orientation
                    if random_pos==length(posits)
                        ori=pos1.get_orientation(1);
                    else
                        ori=pos1.get_orientation(2);
                    end
                    ori_probs=probs(2:end); % Orientation probabilities (right,left,forward,backward)
                    aux=find(rand<cumsum(ori_probs)); % Selects orientation
                    if aux(1)<3 % right or left
                        ori=ori+(-1)^aux(1)*40+randn.*5;
                    else if aux(1)>3 % backward
                            ori=ori-180;
                        end
                    end
                    orientation=[ori ori]; % Sets the orientation
                end
            end
            if ~bound
                ori=orientations(randi(36,1)); % Selects the polymer's plus end orientation
                orientation=[ori ori]; % Sets the orientation
                dim=face.get_dimensions-[1e-16 1e-16]; % Gets the dimensions of the selected face and makes sure the frontier is not included
                minus_end_pos=face.get_position+[1e-16 1e-16]+arrayfun(@(b)b.*rand,dim); % Selects the minus end position
            end
            
            %minus_end_pos=[0.00003989 0.00004502];
            %orientation=315;
            %new_face_ID=6;
            
            %minus_end_pos=[53e-6 10e-6];
            %orientation=45;
            %new_face_ID=2;
            
            plus_end_pos=minus_end_pos; % Sets the plus end position
            pos2=Posit(minus_end_pos,orientation,new_face_ID,0,0,bound*branch_id); % Creates the minus end position
            pos3=Posit(plus_end_pos,orientation,new_face_ID,0,0,bound*branch_id); % Creates the plus end position
            mt2=MT({pos2,pos3},new_face_ID,bound,it); % Creates the polymer (two positions and a known curvature angle)
            [mt2,extinct]=mt2.update_position(2,[vg_p vs_p],faces); % Updates the MT plus end position
            if ~extinct
                arrayfun(@(x)faces{x}.add_MT(mt2),mt2.get_faces_ID); % updates MT objects in faces cell array
                MTs={MTs{1:end} mt2}; % Adds the selected polymer to the list of polymers
            end
        end
    end

end














