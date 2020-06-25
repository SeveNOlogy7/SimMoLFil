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
        
        function testFilterEvaluation(testCase)
            BPF = simulation.Filter("BPF", component.Filter(1030, 1030, 7));
            In = simulation.Input("In");
            Out = In + BPF;
            Eval = Out.build(simulation.builder.Builder());
            
            config = simulation.Configuration();
            input = rand_sech(config.nt, config.time)';
            result = Eval(struct("configuration",config,...
                "In","sech"));
            testCase.assertEqual(true,true);
        end
        
        function testCouplerEvaluation(testCase)
            % Test Coupler model
            createFixture(testCase);
            rho = 0.3;
            In = simulation.Input("In");
            OC = simulation.Coupler("OC", component.Coupler(rho));
            
            [Out1,Out2] = [In, 0] + OC;
            Eval1 = Out1.build(simulation.builder.Builder());
            Eval2 = Out2.build(simulation.builder.Builder());
            
            config = simulation.Configuration();
            input = rand_sech(config.nt, config.time)';
            result1 = Eval1(struct("configuration",config,...
                "In","sech"));
            result2 = Eval2(struct("configuration",config,...
                "In","sech"));
            
            testCase.assertEqual(result1.*conj(result1),rho*input.^2,"RelTol",1e-6);
            testCase.assertEqual(result2.*conj(result2),(1-rho)*input.^2,"RelTol",1e-6);
        end
    end
end

