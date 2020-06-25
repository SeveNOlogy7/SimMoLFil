classdef Evaluation
    %EVALUATION A class that represnts a evaluatable/ a runnable program
    %   Need to be properly constructed by a Builder
    %   Once built, a Evaluation object e is callable, will output
    %   simulation results for dufferent sets of arguments.
    
    properties (Access = private)
        program         % collection of models in order
        variables       % collection of variables presenting light fields
    end
    
    methods
        function obj = Evaluation(program)
            %EVAL Construct an instance of this class
            obj.program = program;
        end
        
        function e = run(obj, kwargs)
            %RUN actually evaluate the program
            % Models in the program are evaluated one by one
            
            % kwargs should be a struct, it is like keyword argument
            if ~isa(kwargs, "struct")
                error("Evaluation: run needs one struct-type argument\n");
            end
            % check parameters in kwargs.
            config = kwargs.configuration;
            % Initialise light field variables
            nt = config.nt;
            obj.variables = zeros(nt,length(obj.program));
            valid_operations = methods("simulation.builder.Evaluation");
            % Start calculation
            ii = 1;     % Model Pointer
            end_of_program = length(obj.program);
            while (ii <= end_of_program)
                model = obj.program(ii);
                operation = model.get_operation();
                op_type = obj.program(ii).get_op_type();
                inputs_id = obj.program(ii).get_inputs_id();
                args = obj.variables(:,inputs_id);
                if ~isempty(find(ismember(valid_operations,op_type), 1))
                    obj.variables(:,ii) = ...
                        simulation.builder.Evaluation.(op_type)...
                        (model, operation, args, kwargs);
                else
                    error("%s operation not implemented\n", op_type);
                end
                % move Model Pointer
                % Should handle Recurrence operation here
                % or maybe inside Recurrence function?
                ii = ii + 1;
            end
            % return the value of the last light_memory block
            e = obj.variables(:,end);
        end
        
        function m = subsref(obj,S)
            %SUBSREF make Evaluation object callable
            if S.type == "()"
                if length(obj) == 1
                    % Evaluation_obj(struct(inputs,parameters...))
                    % Check arguments numbers
                    if length(S.subs)~=1 && isa(S.subs{1},"struct")
                        error("Evaluation: need one struct-type argument\n");
                    end
                    m = obj.run(S.subs{1});
                    return
                end
            elseif S.type == "."
                m = builtin('subsref',obj,S);
                return
            end
            % Everything else results in error
            error("Evaluation: not implemented types of subsref\n");
        end
    end
    
    % Implementation of concrete function for each Operations
    methods (Static)
        function u = Input(model, operation, args, kwargs)
            % Read input from kwargs of Evaluation
            value = kwargs.(operation.name);
            if isa(value,"double")
                u = value;
            elseif isa(value,"string")
                if value == "sech"
                    config = kwargs.("configuration");
                    u = rand_sech(config.nt, config.time)';
                else
                    error("Evaluation: not implemented kwargs value %s for operation %s\n", value, operation);
                end
            else
                error("Evaluation: not implemented kwargs type %s for operation %s\n", class(value), operation);
            end
        end
        
        function u = Const(model, operation, args, kwargs)
            % Read Const value from operation parameters
            value = operation.parameters.("value");
            if isa(value,"double")
                u = value;
            else
                error("Evaluation: not implemented kwargs type %s for operation %s\n", class(value), operation);
            end
        end
        
        function u = Fiber(model, operation, args, kwargs)
            % Read Fiber component value from operation parameters
            config = kwargs.("configuration");
            fiber = operation.parameters.("fiber").cal_gamma(config.lambda0);
            % use GNLSE solver
            [u,~,~] =  IP_CQEM_FD(args',config.dt,config.dz,fiber,config.f0,config.tol,0,0);
            u = u';
        end
        
        function u = Show(model, operation, args, kwargs)
            u = args;
            figure, plot(abs(u));
        end
    end
end

