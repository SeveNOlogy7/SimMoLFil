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
        function testArrayInput(testCase)
            In = simulation.Input("In");
            show = simulation.Show(struct());
            Out = show(In);
            Eval = Out.build(simulation.builder.Builder());
            
            config = simulation.Configuration();
            input = rands(config.nt,1);
            result = Eval(struct("configuration",config,...
                "In",input));
            testCase.assertEqual(result,input);
        end
        
        function testSechInput(testCase)
            In = simulation.Input("In");
            show = simulation.Show(struct());
            Out = show(In);
            Eval = Out.build(simulation.builder.Builder());
            
            config = simulation.Configuration();
            input = rand_sech(config.nt, config.time)';
            result = Eval(struct("configuration",config,...
                "In","sech"));
            testCase.assertEqual(result,input);
        end
        
        function testConst(testCase)
            config = simulation.Configuration();
            const = rands(config.nt,1);
            Out = simulation.Const(const);
            Eval = Out.build(simulation.builder.Builder());
            
            result = Eval(struct("configuration",config));
            testCase.assertEqual(result,const);
        end
        
        function testSingleFiberEvaluation(testCase)
            fiber = component.Fiber();
            fiber.L = 0.001;
            SMF = simulation.Fiber("SMF", fiber);
            In = simulation.Input("In");
            show = simulation.Show(struct());
            Out = show(In + SMF);
            Eval = Out.build(simulation.builder.Builder());
            
            config = simulation.Configuration();
            input = rands(config.nt,1);
            result = Eval(struct("configuration",config,...
                "In",input));
            testCase.assertEqual(true,true);
        end
        
        function testMultipleFiberEvaluation(testCase)
            fiber = component.Fiber();
            fiber.L = 0.001;
            createFixture(testCase);
            SMF1 = simulation.Fiber("SMF1", fiber);
            SMF2 = simulation.Fiber("SMF2", fiber);
            SMF3 = simulation.Fiber("SMF3", fiber);
            In = simulation.Input("In");
            Out = In + SMF1 + SMF2 + SMF3;
            Eval = Out.build(simulation.builder.Builder());
            
            config = simulation.Configuration();
            input = rands(config.nt,1);
            result = Eval(struct("configuration",config,...
                "In",input));
            testCase.assertEqual(true,true);
        end
    end
end

