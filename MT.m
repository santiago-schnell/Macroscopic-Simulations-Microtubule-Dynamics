classdef MT < handle
    % Properties
    properties(GetAccess = 'private', SetAccess = 'private')
        positions={}; % MT position points
        tip_state=[2 1]; % MT tip states for both minus (first position) and plus (second position) tip (1-growth, 2-pause, 3-shortening)
        mt_length=0; % MT length
        mt_length_vars=zeros(1,4); % MT length variations
        state_hist=zeros(2,3); % MT state hist for both tip ends
        consec_states=zeros(2,3); % Stores the number of consecutive states
        block_states=zeros(2,3); % Stores blocks of consecutive states
        nblock_states=zeros(2,3); % Stores the number of added blocks
        handles={}; % MT graphic handles
        faces_ID=[]; % MT faces_ID
        bound=0; % MT bound to another MT
        init_it=1; % Stores the iteration at which the MT nucleated
        interactions=0; % Stores the number of interactions
        ind_cat=[0 0]; % Induced catastrophe
        face_transitions=zeros(6,6); % Stores the face transition count of the MT
    end
    
    % This section defines public methods
    methods(Access='public')
        % Constructor
        function this=MT(new_positions,new_faces_ID,new_bound,new_init_it)
            this.positions=new_positions; % Sets MT initial positions
            this.faces_ID=[this.faces_ID new_faces_ID]; % Sets MT face
            this.bound=new_bound; % Sets MT bound state
            this.init_it=new_init_it; % Sets initial iteration
        end
        
        % This function adds a position to the MT
        function add_position(this,new_pos)
            this.positions{end+1}=new_pos;
        end
        
        % This function returns the MT position
        function res=get_positions(this)
            res=this.positions;
        end
        
        % This function returns the MT minus end position
        function res=get_minus_end_pos(this)
            res=this.positions{1};
        end
        
        % This function returns the MT plus end position
        function res=get_plus_end_pos(this)
            res=this.positions{end};
        end
        
        % This function returns the MT tip states
        function res=get_tip_state(this)
            res=this.tip_state;
        end
        
        % This function returns the MT bound state
        function res=get_bound(this)
            res=this.bound;
        end
        
        % This function returns the length of the MT
        function res=get_mt_length(this)
            res=this.mt_length;
        end
        
        % This function returns the variation length of the MT
        function res=get_mt_length_vars(this)
            res=this.mt_length_vars;
        end
        
        % This function return the iteration at which the MT started
        function res=get_init_it(this)
            res=this.init_it;
        end
        
        % This function return the MT number of interactions
        function res=get_interactions(this)
            res=this.interactions;
        end
        
        % This function returns the MT faces ID
        function res=get_faces_ID(this)
            res=this.faces_ID;
        end
              
        % This function returns the face transition count of the MT
        function res=get_face_transitions(this)
            res=this.face_transitions;
        end
        
        % This function returns the MT state history for both ends
        function res=get_state_hist(this)
            res=this.state_hist;
        end
        
        % This function returns the MT run length
        function res=get_run_length(this)
            res=zeros(2,3);
            for k=1:2
                for l=1:3
                    if (this.tip_state(k)==l)
                        block_state=this.block_states(k,l)+this.consec_states(k,l)+1;
                        nblock_state=this.nblock_states(k,l)+1;
                        res(k,l)=block_state/nblock_state;
                    else
                        if this.nblock_states(k,l)==0
                            res(k,l)=0;
                        else
                            res(k,l)=this.block_states(k,l)./this.nblock_states(k,l);
                        end
                    end
                end
            end
        end
        
        % This function initiates the variation length of the MT
        function this=init_mt_length_vars(this)
            this.mt_length_vars=zeros(1,4);
        end
        
        % This function updates the state of the MT
        function this=update_state(this,tfm,tfp)
            matrix=tfm;
            old_tip_state=this.tip_state; % Stores old state
            for k=1:2 % Update both ends of the MT
                if this.ind_cat(k)==1 % Checks if tip is in catastrophe mode
                    this.tip_state(k)=3; % MT is shortening at the affected tip
                else
                    val=rand;
                    aux_comp=1:3;
                    ets=this.tip_state(k); % Get the end tip state of the MT
                    aux_comp(ets)=[]; % Sets possible comparisons for the transition
                    tf1=matrix(ets,aux_comp(1));
                    if val<=tf1 % Check first transition probability
                        this.tip_state(k)=aux_comp(1); % Change state of the tip to first possible transition
                    else
                        if val<=(tf1+matrix(ets,aux_comp(2))) % Check second transition probability
                            this.tip_state(k)=aux_comp(2); % Change state of the tip to second possible transition
                        end
                    end
                end
                matrix=tfp;
            end
            
            % Releases bound microtubule if a transition to shortening occurs
            if this.bound
                if  this.tip_state(1)==3
                    this.bound=0;
                else
                    this.tip_state(1)=2;
                end
            end
            
            % Updates state hist
            for k=1:2
                this.state_hist(k,this.tip_state(k))=this.state_hist(k,this.tip_state(k))+1;
            end
            
            %  Updates run length info
            for k=1:2
                for l=1:3
                    curr_state=this.tip_state(k);
                    if (curr_state==old_tip_state(k))
                        this.consec_states(k,curr_state)=this.consec_states(k,curr_state)+1;
                    else
                        this.block_states(k,curr_state)=this.block_states(k,curr_state)+this.consec_states(k,curr_state)+1;
                        this.nblock_states(k,curr_state)=this.nblock_states(k,curr_state)+1;
                        this.consec_states(k,curr_state)=0;
                    end
                end
            end
        end
                 
        % This function updates the tip position of a MT
        function [this,extinct,rel_nuc_site]=update_position(this,tip,velocs,faces)
            extinct=0;
            rel_nuc_site=0;
            aux_tip_state=this.tip_state;
     
            % Polymerizes or depolymerizes according to the tip state
            if aux_tip_state(tip)==1
                this.polymerize(tip,velocs(1),faces);
            else if aux_tip_state(tip)==3
                    [extinct,rel_nuc_site]=this.depolymerize(tip,velocs(2),faces);
                end
            end 
        end
        
        % This function releases dimers from the MT tip
        function [extinct rel_nuc_site]=depolymerize(this,tip,veloc,faces)
            tt_dist=0;
            extinct=0;
            rel_nuc_site=0;
            op=(-1)^(tip-1);
            lst_size=length(this.positions);
            prev=1+(tip-1)*(lst_size-1);
            next=prev+op;
            remove_faces_ID=this.positions{prev}.get_face_ID;
            intersection=0;
            dist_over=0;
            nrem=0;
            
            % Find the segments to remove by releasing dimers
            while (tt_dist<=veloc && next>=1 && next<=lst_size)
                if (this.positions{prev}.get_frontier_point && this.positions{next}.get_frontier_point)
                    dist=0; % Frontier points have no distance between them
                    remove_faces_ID=[remove_faces_ID this.positions{next}.get_face_ID];
                    nrem=nrem+1;
                else
                    dist=norm(this.positions{prev}.get_position-this.positions{next}.get_position);
                    if dist==0
                        rel_nuc_site=1;
                    end
                end
                tt_dist=tt_dist+dist;
                prev=prev+op;
                next=prev+op;
            end
            % Microtubule depolymerizes to extinction if next is not within the list boundaries
            if tt_dist<=veloc
                extinct=1;
                this.positions={};
                this.mt_length_vars(tip+2)=this.mt_length;
                %display('Just got extinct...');
            else
                % Determines a new end position for the depolymerized segment
                aux_pos=prev-op;
                orientation=this.positions{aux_pos}.get_orientation(3-tip)-(tip-1)*180;
                old_pos=this.positions{aux_pos}.get_position;
                new_pos=old_pos+(veloc-tt_dist+dist).*[cos(orientation*pi/180) sin(orientation*pi/180)];
                
                % Check if the new point intersects with the cube faces
                % (this is done because of precision issues with the calculation of the new position)
                if out_of_bounds(new_pos)
                    face_ID=this.positions{aux_pos}.get_face_ID;
                    face=faces{face_ID};
                    old_pos=this.positions{aux_pos+op}.get_position;
                    [a,b,c,d]=face.get_rect_points; % Gets the four points that define a face of the cube
                    %a=a+[1e-16 1e-16];
                    %b=b+[1e-16 -1e-16];
                    %c=c+[-1e-16 -1e-16];
                    %d=d+[-1e-16 1e-16];
                    out=lineSegmentIntersect([new_pos(1) new_pos(2) old_pos(1) old_pos(2)],[a b;d c;a d;b c]); % Check intersection
                    int_side=find(out.intAdjacencyMatrix==1); % Gets intersection sides
                    if ~isempty(int_side) % if there is an intersection
                        intersection=1;
                        rand_side=randi(length(int_side));
                        side=int_side(rand_side);
                        dist_over=out.intNormalizedDistance1To2(side)*norm(old_pos-new_pos); % Determine the distance from the new point to the intersection point
                        new_pos=[out.intMatrixX(side) out.intMatrixY(side)]; % Gets the intersection point
       
                        % Some intersection points are not correctly stated by the
                        % lineSegmentIntersect function because of precision issues
                        switch side
                            case 1
                                new_pos(1)=a(1);
                            case 2
                                new_pos(1)=d(1);
                            case 3
                                new_pos(2)=a(2);
                            case 4
                                new_pos(2)=b(2);
                        end
                    end
                end
                
                % Initialize remaining position properties
                this.positions{aux_pos}.set_position(new_pos);
                this.positions{aux_pos}.set_frontier_point(0);
                this.positions{aux_pos}.set_nuc_point(1);
                this.positions{aux_pos}.set_branch_id(0);
                   
                % Remove depolymerized faces from the face list
                if tip==1
                    this.faces_ID(1:nrem)=[];
                else
                    this.faces_ID=this.faces_ID(1:end-nrem);
                end
                % Updates transition matrix
                for m=1:length(remove_faces_ID)-1
                    this.face_transitions=update_transitions(remove_faces_ID(m),remove_faces_ID(m+1),this.face_transitions);
                end
                % Removes segments due to depolymerization
                if tip==1
                    erase_pos=1:aux_pos-1;
                else
                    erase_pos=aux_pos+1:length(this.positions);
                end
                
                % Calculate the exact subtracted length to the MT
                difference=(veloc-tt_dist+dist)-norm(old_pos-new_pos);
                if difference>0
                    sub_length=veloc-difference;
                else
                    sub_length=veloc+difference;
                end
                if intersection
                    sub_length=sub_length-dist_over;
                end
                
                % Updates the MT length
                this.mt_length=this.mt_length-sub_length;
                
                % Checks if the MT is extinct
                if this.mt_length<0
                    extinct=1;
                    this.positions={};
                    this.mt_length_vars(tip+2)=this.mt_length;
                else
                    this.mt_length_vars(tip+2)=sub_length;
                    this.positions(erase_pos)=[];
                end                          
            end
        end
        
        % This function adds dimers to the MT tip
        function polymerize(this,tip,veloc,faces)
            pos_obj=[this.get_minus_end_pos this.get_plus_end_pos]; % Gets the minus and plus end object position
            old_pos=pos_obj(tip).get_position; % Sets old tip position
            branch_id=pos_obj(tip).get_branch_id; % Gets branch id
            curvature_angle=1.*randn; % Selects the polymer curvature (modify later)
            prev_orientation=pos_obj(tip).get_orientation(3-tip); % Gets previous orientation
            new_orientation=prev_orientation+curvature_angle; % Determines new orientation at the MT tip
            new_pos_ori=(new_orientation+(2-tip)*180)*pi/180;
            new_pos=old_pos+veloc.*[cos(new_pos_ori) sin(new_pos_ori)]; % Determines the new tip position
            face_ID=pos_obj(tip).get_face_ID; % Gets the current face ID
            face=faces{face_ID}; % Gets the object face
            
            [new_orientation_intersect,ic,new_pos]=check_MT_intersections(this,tip,new_orientation,veloc,old_pos,branch_id,new_pos,face); % Check if the MT intersects with another MT in the current face
            if ~(new_orientation_intersect==new_orientation)
                prev_orientation=new_orientation;
            end
            new_orientation=new_orientation_intersect;
            if ic % If there was intersection leading to an induced catastrophe
                this.state_hist(tip,1)=this.state_hist(tip,1)-1;
                this.state_hist(tip,2)=this.state_hist(tip,2)+1;
            else
                pos_obj(tip).set_orientation(new_orientation,tip); % Sets orientation in the old position
                final_orientation(tip)=0; % Sets orientation in the new position
                final_orientation(3-tip)=new_orientation; % Sets orientation in the new position
                new_side=[];
                intersection=1;
                inter_lst={};
                % We need to test two intersections whenever one
                % intersection with any side of the cube faces occurs
                while intersection
                    % Find intersection if it occurs
                    [intersection,new_side,norm_dist,new_face_ID,initial_point,int_point,new_orientation]=find_intersection(face_ID,faces,old_pos,new_pos,new_orientation,new_side);
                    if intersection
                        dist=norm_dist*norm(old_pos-new_pos);
                        old_pos=initial_point; % Sets new old_pos
                        new_pos_ori=(new_orientation+(2-tip)*180)*pi/180;
                        new_pos=initial_point+dist.*[cos(new_pos_ori) sin(new_pos_ori)]; % Sets new position
                        final_orientation=[new_orientation new_orientation];
                        % Sets new face ID on the list of MT faces
                        if tip==1
                            this.faces_ID=[new_face_ID this.faces_ID];
                        else
                            this.faces_ID=[this.faces_ID new_face_ID];
                        end
                        this.face_transitions=update_transitions(face_ID,new_face_ID,this.face_transitions); % Updates transition matrix
                        int_point_obj=Posit(int_point,[prev_orientation prev_orientation],face_ID,1,0,0); % Tip sets the int point
                        initial_point_obj=Posit(initial_point,final_orientation,new_face_ID,1,0,0); % Tip sets the initial point
                        % Sets the intersection points for the final list of MT points
                        if tip==1
                            points={initial_point_obj int_point_obj};
                            inter_lst={points{1,:} inter_lst{1,:}};
                        else
                            points={int_point_obj initial_point_obj};
                            inter_lst={inter_lst{1,:} points{1,:}};
                        end
                    end
                    face_ID=new_face_ID; % Sets new_face_ID
                end
                        
                % Creates the new final position and updates the list of points
                new_pos_obj=Posit(new_pos,final_orientation,face_ID,0,1,0); % Tip sets the int point
                if tip==1
                    this.positions={new_pos_obj inter_lst{1,:} this.positions{1,:}}; % Adds to the minus end
                else
                    this.positions={this.positions{1,:} inter_lst{1,:} new_pos_obj}; % Adds to the plus end
                end
                
                % Calculate the exact added length to the MT
                add_length=0;
                op=(-1)^(tip-1);
                pos=1+(tip-1)*(length(this.positions)-1);
                for n=1:length(inter_lst)+1
                    if ~(this.positions{pos}.get_frontier_point && this.positions{pos+op}.get_frontier_point)
                        add_length=add_length+norm(this.positions{pos}.get_position-this.positions{pos+op}.get_position);
                    end
                    pos=pos+op;
                end
                this.mt_length=this.mt_length+add_length;
                this.mt_length_vars(tip)=add_length;
             end
        end
        
        % This function displays the MTs
        function display_MT(this)
            % Eliminates current handles
            for l=1:length(this.handles)
                set(this.handles{l},'xdata',[0 0],'ydata',[0 0]);
            end
            this.handles={};
            % Creates new handles
            for l=2:length(this.positions)
                pos1_obj=this.positions{l-1};
                pos2_obj=this.positions{l};
                pos1=pos1_obj.get_position;
                pos2=pos2_obj.get_position;
                if ~(pos1_obj.get_frontier_point && pos2_obj.get_frontier_point)
                    h1=line([pos1(1) pos2(1)],[pos1(2) pos2(2)],'color','r','erase','xor','linewidth',1);
                    this.handles{end+1}=h1;
                end
            end
        end
        
        % This function severes the MT according to the specified probabilities
        function [this,newMT]=severe(this,sfprobs,gh)
            newMT=-1;
            inds=cellfun(@(x)((x.get_branch_id==0) && (~x.get_frontier_point)),this.positions); % Select positions
            inds(1:2)=0;
            inds(end-2:end)=0;
            aux=find(inds==1);
            if ~isempty(aux)
                ind_pos=aux(randi(length(aux)));
                pos=this.positions{ind_pos}; % Gets the selected position in the MT
                if (rand<sfprobs(pos.get_face_ID))
                    % Eliminates handles of the original MT
                    if gh
                        for l=1:length(this.handles)
                            set(this.handles{l},'xdata',[0 0],'ydata',[0 0]);
                        end
                    end
                    new_positions=copy_Posits(this.positions(ind_pos:end));
                    new_faces_ID=unique(cellfun(@(x)(x.get_face_ID),new_positions));
                    newMT=MT(new_positions,new_faces_ID,0,this.init_it);
                    newMT=copy_Props(newMT,this);
                    newMT.tip_state(1)=3;
                    % Makes changes to the original MT
                    this.positions(ind_pos+1:end)=[];
                    this.tip_state(2)=3;
                end
            end
        end
        
    end
end

% ------------- This section sets auxiliary functions --------------

% This function copies designated positions into a new vector
function new_Posits=copy_Posits(positions)
    new_Posits={};
    % Creates the new positions
    for i=1:length(positions)
        new_pos=positions{i}.get_position;
        new_ori=positions{i}.get_ori;
        new_face_ID=positions{i}.get_face_ID;
        new_frontier_point=positions{i}.get_frontier_point;
        new_nuc_point=positions{i}.get_nuc_point;
        new_branch_id=positions{i}.get_branch_id;
        new_pos=Posit(new_pos,new_ori,new_face_ID,new_frontier_point,new_nuc_point,new_branch_id);
        new_Posits{end+1}=new_pos;
    end
end

% This function copies the properties of one MT to another MT
function newMT=copy_Props(newMT,MT)
    newMT.tip_state=MT.tip_state;
    newMT.mt_length_vars=MT.mt_length_vars;
    newMT.state_hist=MT.state_hist;
    newMT.consec_states=MT.consec_states;
    newMT.block_states=MT.block_states;
    newMT.handles=MT.handles;
    newMT.interactions=MT.interactions;
    newMT.ind_cat=MT.ind_cat;
    newMT.face_transitions=MT.face_transitions;
end

% This function checks if the MT intersects with another MT in the current face
function [new_orientation,ic,new_pos]=check_MT_intersections(MT,tip,new_orientation,veloc,old_pos,branch_id,new_pos,face)
ic=0;
face_MTs=face.get_MTs();
MT_coords=[old_pos(1) old_pos(2) new_pos(1) new_pos(2)];
for m=1:length(face_MTs)
    if ~eq(MT,face_MTs{m})
        posits=face_MTs{m}.get_positions;
        other_MT_coords=[];
        for n=1:length(posits)
            pos=posits{n};
            if (~branch_id || branch_id~=pos.get_branch_id)
                other_MT_coords=[other_MT_coords pos.get_position];
            end
        end
        if length(other_MT_coords)>2 % This is only effective when no workers are being used
            [X0,Y0]=intersections(MT_coords(1,1:2:end)',MT_coords(1,2:2:end)',other_MT_coords(1,1:2:end)',other_MT_coords(1,2:2:end)',1);
            if ~isempty(X0)
                %display('I JUST INTERSECTED...')
                aux1=other_MT_coords(1,1:2:end)-X0(end);
                aux2=other_MT_coords(1,2:2:end)-Y0(end);
                dists=zeros(1,length(aux1));
                for n=1:length(aux1)
                    dists(n)=norm([aux1(n) aux2(n)]);
                end
                other_MT_int_position=face_MTs{m}.get_positions{dists==min(dists)};
                [new_orientation,interaction,ic]=act_intersect(tip,face_MTs{m},new_orientation,other_MT_int_position);
                MT.interactions=MT.interactions+interaction;
                MT.ind_cat(tip)=ic;
                break;
            end
        end
    end
end
new_pos_ori=(new_orientation+(2-tip)*180)*pi/180;
new_pos=old_pos+veloc.*[cos(new_pos_ori) sin(new_pos_ori)]; % Determines the new tip position
end

% This function acts on the MT intersection
function [new_orientation,interaction,ic]=act_intersect(tip,other_MT,new_orientation,other_MT_int_pos)
ic=0;
interaction=0;

% Gets orientation of plus end
if ~eq(other_MT_int_pos,other_MT.get_plus_end_pos)
    other_MT_orient_plus_end=other_MT_int_pos.get_orientation(2);
else
    other_MT_orient_plus_end=other_MT_int_pos.get_orientation(1);
end
other_MT_orient_minus_end=other_MT_orient_plus_end-180; % Gets orientation of minus end

angle_diff_p=new_orientation+(2-tip)*180-other_MT_orient_plus_end; % Gets angle difference in relation to the plus end
angle_diff_p=abs(angle_diff_p-360*(round(angle_diff_p/360))); % Converts angle
prob_bundle_p=1-(1./(1.+exp(-(angle_diff_p/5-8)))); % Gets probability of bundle with plus end
val=rand; % Get random number from uniform distribution
if val<prob_bundle_p
    new_orientation=other_MT_orient_plus_end+(2-tip)*180;
    interaction=1;
    %display('JUST BUNDLED 1...')
else
    angle_diff_m=new_orientation+(2-tip)*180-other_MT_orient_minus_end; % Gets angle difference in relation to the minus end
    angle_diff_m=abs(angle_diff_m-360*(round(angle_diff_m/360))); % Converts angle
    prob_bundle_m=1-(1./(1.+exp(-(angle_diff_m/5-8)))); % Gets probability of bundle with minus end
    if val<prob_bundle_m
        new_orientation=other_MT_orient_minus_end+(2-tip)*180;
        interaction=1;
        %display('JUST BUNDLED 2...')
    else
        prob_cat=1./(0.8*sqrt(2*pi))*exp(-0.5*((angle_diff_p-90)/0.8).^2); % Gets catastrophe induced depolymerization probability
        if rand < prob_cat
            %display('Induced catastrophe activated...')
            ic=1; % Activates induced catastrophe
            interaction=1;
        end
    end
end
end

% This function finds the intersection with any side of a cube face
% It also returns info such as the intersection point, initial point
% and new orientation
function [intersection,new_side,norm_dist,new_face_ID,initial_point,int_point,new_orientation]=find_intersection(face_ID,faces,old_pos,new_pos,orientation,new_side)
% Initiate return values
intersection=0;
norm_dist=0;
new_face_ID=face_ID;
initial_point=[0 0];
int_point=[0 0];
new_orientation=orientation;

face=faces{face_ID};
face_position=face.get_position; % Gets position of the current face
[a,b,c,d]=face.get_rect_points; % Gets the four points that define a face of the cube
out=lineSegmentIntersect([new_pos(1) new_pos(2) old_pos(1) old_pos(2)],[a b;d c;a d;b c]); % Check intersection
out.intAdjacencyMatrix(new_side)=0; % Eliminates known intersections (used for the possible second intersection)
int_side=find(out.intAdjacencyMatrix==1); % Gets intersection sides

if ~isempty(int_side) % if there is an intersection
    intersection=1;
    rand_side=randi(length(int_side)); % Randomly chooses one of the sides (if there is more than one)
    side=int_side(rand_side); % Gets the identifier of the side (1 for < x, 2 for > x, 3 < y and 4 for > y)
    face_transitions=face.get_face_transitions; % Gets the transition vector for the current face
    new_face_ID=face_transitions(side); % Gets the face_ID to which we are going to make the transition
    side_transitions=face.get_side_transitions(); % Gets the side transition vector
    new_side=side_transitions(new_face_ID); % Gets the side to which we are going to make the transition
    norm_dist=out.intNormalizedDistance1To2(side); % Determine the normalized distance from the new point to the intersection point
    trans_ID=[face_ID new_face_ID];
    diff_trans_ID=diff(trans_ID);
    corner=length(int_side)-1; % Determines if the intersection point is a corner
    rotation=((abs(diff_trans_ID)>1) && (length(intersect(trans_ID,[1 3]))<2));
    new_face=faces{new_face_ID}; % Gets the face object of the face for which we are making the transition to
    new_face_position=new_face.get_position; % Gets the face position from the face we are making the transition to
    int_point=[out.intMatrixX(side) out.intMatrixY(side)]; % Gets the intersection point
    
    % Some intersection points are not correctly stated by the
    % lineSegmentIntersect function because of precision issues
    switch side
        case 1 
            int_point(1)=a(1);
        case 2 
            int_point(1)=d(1);
        case 3 
            int_point(2)=a(2);
        case 4 
            int_point(2)=b(2);
    end
    
    if ~rotation
        initial_point=int_point; % Determines the initial point of the new face
    else
        initial_point=det_new_pos(face_ID,new_face_ID,new_face_position,int_point,face_position); % Determines the initial point of the new face
    end
    if rotation || corner
        new_orientation=det_new_orientation(corner,new_face,initial_point,face_ID,new_face_ID,orientation); % Determines the orientation on the new face
    end
end
end

% This function determines the new orientation in case of an intersection
function new_orientation=det_new_orientation(corner,new_face,initial_point,face_ID,new_face_ID,orientation)
if corner % Determines the new orientation if intersection is a corner
    [~,b,c,d]=new_face.get_rect_points;
    d1=abs(det([d-c;initial_point-c]))/norm(d-c);
    d2=abs(det([c-b;initial_point-b]))/norm(c-b);
    %[cos(orientation*pi/180)<0 sin(orientation*pi/180)<0];
    res=abs([d1<1e-6 d2<1e-6]-[cos(orientation*pi/180)<0 sin(orientation*pi/180)<0]);
    add_angle=90*sum(res);
    if xor(res(1),res(2)) add_angle=add_angle*diff(res); end
    new_orientation=orientation+add_angle;
else if (length(intersect([face_ID new_face_ID],[3 6]))==2)
        new_orientation=orientation;
    else if (face_ID<3)
            new_orientations=(-1)^(face_ID).*(1:3).*90;
            new_orientation=orientation+new_orientations(new_face_ID-3);
        else
            new_orientations=(-1)^(new_face_ID-1).*(1:3).*90;
            new_orientation=orientation+new_orientations(face_ID-3);
        end
    end
end
end

% This function determines the new initial point in case of an intersection
function initial_point = det_new_pos(face_ID,new_face_ID,new_face_position,int,a)
if (length(intersect([face_ID new_face_ID],[3 6]))==2)
    if face_ID<new_face_ID
        initial_point=[int(1) 60e-6];
    else
        initial_point=[int(1) 0];
    end
else if (face_ID<3)
        trans_x=[0 0 0]+(face_ID-1).*[40 40 40].*1e-6;
        trans_y=[abs(int(1)) 15e-6-int(2) 15e-6-abs(int(1))]+(face_ID-1).*[-40 0 40].*1e-6;
        initial_point=[trans_x(new_face_ID-3) trans_y(new_face_ID-3)+new_face_position(2)];
    else
        aux=int(2)-a(2);
        trans_x=[-aux -15e-6 -15e-6+aux]+(new_face_ID-1).*[40e-6+2*aux 70e-6 70e-6-2*aux];
        trans_y=[15e-6 15e-6-aux 0];
        initial_point=[trans_x(face_ID-3) trans_y(face_ID-3)];
    end
end
end

% This function updates the matrix of face transitions
function face_transitions=update_transitions(face_ID,new_face_ID,face_transitions)
if ~mod(face_transitions(new_face_ID,face_ID),1)
    face_transitions(face_ID,new_face_ID)=face_transitions(face_ID,new_face_ID)+0.5;
else
    face_transitions(new_face_ID,face_ID)=face_transitions(new_face_ID,face_ID)-0.5;
end
end

% This function tells if the new position is out of bounds (used by the
% depolymerization method)
function res=out_of_bounds(new_pos)
res=0;
if  ((new_pos(1)<-15e-6) || (new_pos(1)>55e-6) || (new_pos(2)<0) || (new_pos(2)>60e-6) ||...
     (new_pos(1)<0 && new_pos(2)>15e-6) || (new_pos(1)>40e-6 && new_pos(2)>15e-6))
    res=1;
end
end





















