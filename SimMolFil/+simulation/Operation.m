classdef Operation
    %OPERATION A class to handle operations required by a simulation model
    %   Operation is a abstraction of components/input/output/flow control
    %   of a simulation model. Each node (which is presented by Model)
    %   contains a specific form of Operation.
    
    properties (Access = private)
        name
        op_type
        num_args
        parameters
    end
    
    methods
        function obj = Operation(name, op_type, num_args, parameters)
            %OPERATION Construct an instance of this class
            obj.name = name;
            obj.op_type = op_type;
            obj.num_args = num_args;
            obj.parameters = parameters;
        end
        
        % overload subsref, make operation object callable
        % also handle 1d array sub reference
        function m = subsref(obj,S)
            if S.type == "()"
                if length(S.subs) == 1 && isa(S.subs{:},"double")
                    % Operation_obj_array(ii) 
                    % This is just a normal 1d array subsref
                    % use builtin subsref here
                    m = builtin('subsref',obj,S);
                else
                    if length(obj) == 1
                        % Operation_obj(cell_array of Models)
                        % Call a single Operation object
                        % Check arguments numbers
                        if obj.num_args >=0 && obj.num_args ~= length(S.subs)
                            error("%s: need %d arguments but found %d\n", obj.name, obj.num_args, length(S.subs));
                        end
                        % Check argument types
                        if ~all(cellfun(@(c) isa(c,"simulation.Model"),S.subs))
                            error("%s: need simulation.Model arguments\n", obj.name);
                        end
                        % S.subs contains only Model objects now
%                         s = cellfun(@char,S.subs);
%                         disp(sprintf("You called %s with argument: ", obj.name) + join(s,","));
                        if obj.op_type == "Feedback"
                            % Special Type check for Feedback
                            L_operand_type = S.subs{1}.get_op_type();
                            if L_operand_type ~= "Recurrence"
                                error("Feedback need Reccurence type input\n");
                            end
                            % t = Feedback(l,r)
                            % Feedback result in a same Recurrence model l
                            % with a new input replaced by the r operand of
                            % feedback
                            m = S.subs{1}.set_inputs(S.subs{2});
                        else
                            m = simulation.Model(obj,[S.subs{:}]);
                        end
                    end
                end
            else
                % Everything else just use builtin for now
                m = builtin('subsref',obj,S);
            end
        end
        
        function disp(obj)
            disp(char(obj));
        end
        
        function s = char(obj)
            if length(obj) == 1
                % convert a object to string
                s = obj.name+"."+obj.op_type;
                if isempty(fieldnames(obj.parameters))
                    return
                end
                % put in parameters
                params = [];
                parameter_names = fieldnames(obj.parameters);
                for ii = 1:length(parameter_names)
                    param = obj.parameters.(parameter_names{ii});
                    if isa(param,"double") || isa(param,"string")
                        params = [params, parameter_names{ii}+"="+param];
                    else
                        params = [params, parameter_names{ii}+"="+class(param)];
                    end
                    
                end
                params = join(params,",");
                s = sprintf("%s[%s]",s,params);
            else
                % convert 1d array of objects to string
                s = [];
                for ii = 1:length(obj)
                    s = [s,char(obj(ii))];
                end
                s = join(s,newline);
            end
        end
        
    end
end

