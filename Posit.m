classdef Posit < handle
    % Properties
    properties(GetAccess = 'private', SetAccess = 'private')
        position=[]; % Position vector
        orientation=[0 0]; % Position orientation
        face_ID=0; % Face to which position belongs to
        nuc_point=1; % States if the point is a nucleation point
        frontier_point=0; % States if the point is a frontier point
        branch_id=0; % States the branch_id for branched nucleations
        bleached=0; % States if the position is bleached (FRAP simulation)
    end
    methods
        % Constructor
        function this=Posit(new_position,new_orientation,new_face_ID,new_frontier_point,new_nuc_point,new_branch_id)
            this.position=new_position; % Sets position
            this.orientation=new_orientation; % Sets orientation
            this.face_ID=new_face_ID; % Sets face_ID
            this.frontier_point=new_frontier_point; % Sets frontier point
            this.nuc_point=new_nuc_point; % Sets nucleation point
            this.branch_id=new_branch_id; % Sets branch id
        end
        
        % This function returns the position
        function res=get_position(this)
            res=this.position;
        end
        
        % This function returns the point orientation
        function res=get_ori(this)
            res=this.orientation;
        end
        
        % This function returns the point orientation at position position
        function res=get_orientation(this,position)
            res=this.orientation(position);
        end
        
        % This function returns the face ID
        function res=get_face_ID(this)
            res=this.face_ID;
        end
        
        % This function returns a boolean stating if the position is a nucleation point
        function res=get_nuc_point(this)
            res=this.nuc_point;
        end
        
        % This function returns a boolean stating if the position is a frontier point
        function res=get_frontier_point(this)
            res=this.frontier_point;
        end
        
        % This function returns the id of the branch
        function res=get_branch_id(this)
            res=this.branch_id;
        end
        
        % This function returns if the position is bleached
        function res=get_bleached(this)
            res=this.bleached;
        end
        
        % This function sets the position
        function set_position(this,new_position)
            this.position=new_position;
        end
        
        % This function sets the orientation
        function set_orientation(this,new_orientation,new_pos)
            this.orientation(new_pos)=new_orientation;
        end
        
        % This function sets the face
        function set_face_ID(this,new_face_ID)
            this.face_ID=new_face_ID;
        end
        
        % This function sets the point as a nucleation point
        function set_nuc_point(this,new_nuc_point)
            this.nuc_point=new_nuc_point;
        end
        
        % This function sets the point as a frontier point
        function set_frontier_point(this,new_frontier_point)
            this.frontier_point=new_frontier_point;
        end
        
        % This function sets the point's branch id
        function set_branch_id(this,new_branch_id)
            this.branch_id=new_branch_id;
        end
        
        % This function sets the bleached property
        function set_bleached(this,new_bleached)
            this.bleached=new_bleached;
        end
    end
end


