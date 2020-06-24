classdef Builder < handle
    %BUILDER A builder class that convert a Model into Evaluation
    % append method, add a model to the program that is being built
    % build method, return a Evaluation based on current program
    
    properties (Access = private)
        program
    end
    
    methods
        function obj = Builder()
            %BUILDER Construct an instance of this class
            obj.program = [];
        end
        
        function append(obj, model)
            % APPEDND, add a model to the program that is being built
            obj.program = [obj.program; model];
        end
        
        function eval = build(obj)
            % BUILD, return a Evaluation based on current program
            eval = simulation.builder.Evaluation(obj.program);
        end
    end
end

