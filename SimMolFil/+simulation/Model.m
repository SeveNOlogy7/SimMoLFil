classdef Model
    %MODEL A class to handle model information as a data flow graph(DFG)
    %   A DFG is a graph that present complex data flow as nodes and
    %   directed edges. Each node has input and output, and performs some
    %   operation. Use DFG to describle a simulation model. This class aims
    %   to construct a DFG presentation of a model using some MATLAB
    %   build-in arithmatic operators.
    %
    %   Technically, a node is a smallest model unit. A full model is a
    %   combination of model units, hence a DFG. For any node in the model,
    %   call statements to display the generating DFG for that node, which
    %   is basicly a sub-model. A DFG is displayed in multilines of texts
    %   each presents a single operation.
    %
    %   Here is an example:
    %       fiber = component.Fiber();  % Define fiber component
    %       SMF1 = simulation.Fiber("SMF1", fiber)  % Operation of node SMF1
    %       SMF2 = simulation.Fiber("SMF2", fiber)  % Operation of node SMF2
    %       SMF3 = simulation.Fiber("SMF3", fiber)  % Operation of node SMF3
    %       In = simulation.Input("In") % A node(Model) that asks for input
    %       Out = In + SMF1 + SMF2 + SMF3;  % An output defined using +
    %       Simulation = Out.statements()   % Print the DFG
    %   This example shows how to constuct a model that cascades 3 sections
    %   of optical fiber.
    %
    %   output of the statements function.
    %       t0 = In.Input()
    %       t1 = SMF1.Fiber[fiber=component.Fiber](t0)
    %       t2 = SMF2.Fiber[fiber=component.Fiber](t1)
    %       t3 = SMF3.Fiber[fiber=component.Fiber](t2)
    %
    %   See unittests for more examples
    
    properties (Access = private)
        operation
        inputs
        id
    end
    
    methods
        function obj = Model(operation, inputs)
            %MODEL Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0
                error("No Default Constructor in Model Class");
            elseif nargin == 1
                obj = simulation.Model.promote(operation);
                if isempty(obj)
                    error("No Compatible Convesion Constructor in Model Class");
                end
                return
            elseif nargin > 2
                error("Too Many parameters for Constructor in Model Class")
            end
            
            % Check Type
            if ~isa(operation,'simulation.Operation')
                error("Type error: 'operation' need type simulation.Operation, get %s",class(operation));
            end
            if ~isempty(inputs) && ~isa(inputs,'simulation.Model')
                error("Type error: 'inputs' need type simulation.Model, get %s",class(inputs));
            end
            
            persistent next_id;
            if isempty(next_id)
                next_id = 0;
            end
            obj.operation = operation;
            obj.inputs = inputs;
            obj.id = next_id;
            next_id = next_id + 1;
        end
        
        function s = statements(obj)
            %STATEMENTS convert an model into readable statement
            %   each statement corresponds to a node in the DFG
            lines = string([]);
            function visitor(o)
                s_temp = char(o);
                % If some inputs are reused, there will be duplicate
                % statements found by dfs, should avoid duplication
                if isempty(find(lines == s_temp, 1))
                    lines = [lines, char(o)];
                end
            end
            obj.depth_first_search({},@visitor);
            s = join(lines,newline);
        end
        
        function varargout = plus(a,b)
            %PLUS Overload operator+
            % This overloads operator+ so that a expression like
            % Model+Operation will perform Operation on a Model
            % output. Ideally an operator like -> looks better.
            % But no operator -> in matlab, use + instead.
            % Model->Operation
            arguments
                a simulation.Model
                b simulation.Operation
            end
            for ii = 1:length(b)
                op = b(ii);
                inputs = cell(a);
                varargout{ii} = op(inputs{:});
            end
        end
        
        function eval = build(obj, builder)
            %BUILD convert an model into evaluation
            lines = string([]);
            function visitor(o)
                s_temp = char(o);
                if isempty(find(lines == s_temp, 1))
                    lines = [lines, char(o)];
                    builder.append(o);
                end
            end
            obj.depth_first_search({}, @visitor);
            eval = builder.build();
        end
        
        % overload cell call
        function c = cell(obj)
            c = cell(1,length(obj));
            for ii = 1:length(obj)
                c{ii} = obj(ii);
            end
        end
        
        function disp(obj)
            disp(char(obj));
        end
        
        function s = char(obj)
            if length(obj) == 1
                args = strings(length(obj.inputs),1);
                for ii = 1:length(obj.inputs)
                    args(ii) = sprintf("t%d", obj.inputs(ii).id);
                end
                if isempty(args)
                    args = "";
                else
                    args = join(args,",");
                end
                s = sprintf("t%d = %s(%s)",obj.id, char(obj.operation),args);
            else
                % convert 1d array of objects to string
                s = strings(length(obj),1);
                for ii = 1:length(obj)
                    s(ii) = char(obj(ii));
                end
                s = join(s,newline);
            end
        end
        
        % getters
        function t = get_operation(obj)
            t = obj.operation;
        end
        function t = get_op_type(obj)
            t = obj.operation.op_type;
        end
        
        function t = get_inputs_id(obj)
            t = zeros(1,length(obj.inputs));
            for ii = 1:length(t)
                t(ii) = obj.inputs(ii).id+1;    % to MATLAB convention
            end
        end
        
        % setters
        function obj = set_inputs(obj, inputs)
            if ~isempty(inputs) && ~isa(inputs,'simulation.Model')
                error("Type error: 'inputs' need type simulation.Model, get %s",class(inputs));
            end
            obj.inputs = inputs;
        end
    end
    
    methods (Access = private)
        function depth_first_search(obj, ids, visitor)
            ids = containers.Map(obj.id, true);
            for ii = 1: length(obj.inputs)
                if isKey(ids,obj.inputs(ii).id)
                    continue
                end
                obj.inputs(ii).depth_first_search(ids, visitor)
            end
            visitor(obj);
        end
    end
    
    methods (Static)
        function p = promote(r)
            if isa(r,"simulation.Model")
                p = r;
            elseif isa(r,"double")
                p = simulation.Const(r);
            end
        end
    end
end

