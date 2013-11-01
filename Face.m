classdef Face < handle
    % Properties
    properties(GetAccess = 'private', SetAccess = 'private')
        ID=0; % Face identifier (1-6)
        position=[]; % Position of the face
        dimensions=[]; % Dimensions of the face
        face_transitions=[]; % Specifies face transitions
        side_transitions=[]; % Specifies side transitions
        MTs={}; % Stores the MTs currently in this face
    end
    methods
        % Constructor
        function this=Face(new_ID,new_position,new_dimensions,new_face_transitions,new_side_transitions)
            this.ID=new_ID; % Sets ID
            this.position=new_position.*1e-6; % Sets position
            this.dimensions=new_dimensions.*1e-6; % Sets dimensions
            this.face_transitions=new_face_transitions; % Sets face transitions
            this.side_transitions=new_side_transitions; % Sets side transitions
        end
        
        % This function returns the Face ID
        function res=get_ID(this)
            res=this.ID;
        end
        
        % This function returns the Face position
        function res=get_position(this)
            res=this.position;
        end
        
        % This function returns Face dimensions
        function res=get_dimensions(this)
            res=this.dimensions;
        end
        
        % This function returns Face transitions
        function res=get_face_transitions(this)
            res=this.face_transitions;
        end
        
        % This function returns Side transitions
        function res=get_side_transitions(this)
            res=this.side_transitions;
        end
        
        % This function returns the Microtubules
        function res=get_MTs(this)
            res=this.MTs;
        end
        
        % This function sets the Microtubules
        function set_MTs(this,new_MTs)
            this.MTs=new_MTs;
        end
        
        % This function adds a new MT to the Face
        function add_MT(this,new_MT)
            this.MTs{end+1}=new_MT;
        end
        
        % This function removes MTs from the face
        function remove_MTs(this,positions)
            this.MTs(positions)=[];
        end
        
        % This function returns the 4 distinct points of a cube face
        function [a b c d]=get_rect_points(this)
            a=this.position;
            c=a+this.dimensions;
            b=[a(1) c(2)];
            d=[c(1) a(2)];
        end
    end
end