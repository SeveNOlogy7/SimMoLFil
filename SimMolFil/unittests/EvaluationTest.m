classdef EvaluationTest < matlab.unittest.TestCase
    %EVALUATIONTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(TestMethodSetup)
        function createFixture(testCase)
            function clearModel()
                clear Model
            end
            testCase.addTeardown(@clearModel);
        end
    end
    
    methods (Test)
        function testSingleFiberEvaluation(testCase)
            % disp, char, + should all work in this test
            fiber = component.Fiber();
            SMF = simulation.Fiber("SMF", fiber);
            In = simulation.Input("In");
            show = simulation.Show(struct());
            Out = show(In + SMF);
            Eval = Out.build(simulation.builder.Builder());
            Eval(struct("configuration",simulation.Configuration()));
            testCase.assertEqual(true,true);
        end
    end
end

