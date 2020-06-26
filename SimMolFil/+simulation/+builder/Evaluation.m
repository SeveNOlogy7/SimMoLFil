classdef Evaluation < handle
    %EVALUATION A class that represnts a evaluatable/ a runnable program
    %   Need to be properly constructed by a Builder
    %   Once built, a Evaluation object e is callable, will output
    %   simulation results for dufferent sets of arguments.
    
    properties (Access = private)
        program         % collection of models in Map
        ids             % id of models in order
        variables       % collection of variables presenting light fields
        iter            % iteration counter for recurrence node
        pointer         % index point to id of current model
    end
    
    methods
        function obj = Evaluation(program)
            %EVAL Construct an instance of this class
            ids = zeros(length(program),1);
            for ii = 1:length(program)
                ids(ii) = program(ii).get_id();
            end
            if length(program)>1
                program = cell(program);
            end
            obj.program = containers.Map(ids,program);
            obj.ids = ids;
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
            % Allocate maximum space for variables
            ones_temp = zeros(nt,max(cell2mat(keys(obj.program)))+1);
            obj.variables = complex(ones_temp, ones_temp);
            obj.iter = 0;
            % Start calculation
            obj.pointer = 1;     % Model Pointer
            end_of_program = length(obj.ids);
            while (obj.pointer <= end_of_program)
                id = obj.ids(obj.pointer);
                model = obj.program(id);
                operation = model.get_operation();
                op_type = obj.program(id).get_op_type();
                % Get inputs idx in variables
                args_idx = obj.program(id).get_inputs_id()+1;
                args = obj.variables(:,args_idx);
                % Throw error if op_type not implemented
                disp(model);
                obj.(op_type)(model, operation, args, kwargs);
            end
            % return the value corresponding to the last model
            if isfield(kwargs, "Out_id")
                e = obj.variables(:,kwargs.Out_id);
            else
                e = obj.variables(:,obj.ids(end)+1);
            end
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
    methods
        function Input(obj, model, operation, args, kwargs)
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
            obj.variables(:,obj.ids(obj.pointer)+1) = u;
            obj.pointer = obj.pointer + 1;
        end
        
        function Const(obj, model, operation, args, kwargs)
            % Read Const value from operation parameters
            value = operation.parameters.("value");
            if isa(value,"double")
                u = value;
            else
                error("Evaluation: not implemented kwargs type %s for operation %s\n", class(value), operation);
            end
            obj.variables(:,obj.ids(obj.pointer)+1) = u;
            obj.pointer = obj.pointer + 1;
        end
        
        function Fiber(obj, model, operation, args, kwargs)
            % Read Fiber component value from operation parameters
            config = kwargs.("configuration");
            fiber = operation.parameters.("fiber").cal_gamma(config.lambda0);
            % use GNLSE solver
            [u,~,~] =  IP_CQEM_FD(transpose(args),config.dt,config.dz,fiber,config.f0,config.tol,0,1);
            obj.variables(:,obj.ids(obj.pointer)+1) = transpose(u);
            obj.pointer = obj.pointer + 1;
        end
        
        function Show(obj, model, operation, args, kwargs)
            config = kwargs.("configuration");
            figure(10), 
            plot(config.t, abs(args).^2,'r.-');axis tight;
            pause(0.1);
            obj.variables(:,obj.ids(obj.pointer)+1) = args;
            obj.pointer = obj.pointer + 1;
        end
        
        function Filter(obj, model, operation, args, kwargs)
            config = kwargs.("configuration");
            filter = operation.parameters.("filter");
            u = filter.filter_gauss(args, config.f0, config.df);
            obj.variables(:,obj.ids(obj.pointer)+1) = u;
            obj.pointer = obj.pointer + 1;
        end
        
        function CouplerOut1(obj, model, operation, args, kwargs)
            config = kwargs.("configuration");
            coupler = operation.parameters.("coupler");
            u = coupler.couplerOut1(args(:,1),args(:,2));
            obj.variables(:,obj.ids(obj.pointer)+1) = u;
            obj.pointer = obj.pointer + 1;
        end
        
        function CouplerOut2(obj, model, operation, args, kwargs)
            config = kwargs.("configuration");
            coupler = operation.parameters.("coupler");
            u = coupler.couplerOut2(args(:,1),args(:,2));
            obj.variables(:,obj.ids(obj.pointer)+1) = u;
            obj.pointer = obj.pointer + 1;
        end
        
        function Switch(obj, model, operation, args, kwargs)
            persistent First_run_flag
            condition = operation.parameters.("condition");
            if condition == "First_Run"
                if isempty(First_run_flag)
                    u = args(:,1);
                    First_run_flag = true;
                else
                    u = args(:,2);
                end
            end
            obj.variables(:,obj.ids(obj.pointer)+1) = u;
            obj.pointer = obj.pointer + 1;
        end
        
        function Recurrence(obj, model, operation, args, kwargs)
            config = kwargs.("configuration");
            condition = operation.parameters.("condition");
            if condition == "Run_until_stop"
                if obj.iter == config.N_trip
                    di = length(obj.program); % jump out of the main loop
                else
                    di = 1;
                end
                obj.iter = obj.iter + 1;
            end
            obj.variables(:,obj.ids(obj.pointer)+1) = args;
            obj.pointer = obj.pointer + di;
        end
        
        function Feedback(obj, model, operation, args, kwargs)
            inputs_id = model.get_inputs_id();
            recurrence_model = obj.program(inputs_id(1));
            % Hack into the recurrence model, reconnect input
            obj.program(inputs_id(1)) = recurrence_model.set_inputs(obj.program(inputs_id(2)));
            % Jump to where the recurrence model is
            obj.pointer = find(obj.ids==inputs_id(1));
        end
        
    end
end

