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
            config = simulation.Configuration();
            fiber = component.Fiber();
            fiber.L = 0.1;
            fiber = fiber.cal_gamma(config.lambda0);
            
            SMF = simulation.Fiber("SMF", fiber);
            In = simulation.Input("In");
            Out = In + SMF;
            Eval = Out.build(simulation.builder.Builder());
            
            input = rands(config.nt,1);
            result = Eval(struct("configuration",config,...
                "In",input));
            
            compare = IP_CQEM_FD(input',config.dt,config.dz,fiber,config.f0,config.tol,0,1);
            testCase.assertEqual(result,transpose(compare));
        end
        
        function testActiveFiberEvaluation(testCase)
            config = simulation.Configuration();
            
            fiber = component.ActiveFiber();
            fiber.L = 0.001;
            fiber = fiber.cal_gamma(config.lambda0)...
                .set_RepeatFre(fiber.L);
            
            SMF = simulation.Fiber("SMF", fiber);
            In = simulation.Input("In");
            Out = In + SMF;
            Eval = Out.build(simulation.builder.Builder());

            input = rands(config.nt,1);
            result = Eval(struct("configuration",config,...
                "In",input));
            
            compare = IP_CQEM_FD(input',config.dt,config.dz,fiber,config.f0,config.tol,0,1);
            testCase.assertEqual(result,transpose(compare));
        end
        
        function testMultipleFiberEvaluation(testCase)
            config = simulation.Configuration();
            
            fiber = component.Fiber();
            fiber.L = 0.001;
            fiber = fiber.cal_gamma(config.lambda0);
            
            createFixture(testCase);
            SMF1 = simulation.Fiber("SMF1", fiber);
            SMF2 = simulation.Fiber("SMF2", fiber);
            SMF3 = simulation.Fiber("SMF3", fiber);
            In = simulation.Input("In");
            Out = In + SMF1 + SMF2 + SMF3;
            Eval = Out.build(simulation.builder.Builder());
            
            input = rands(config.nt,1);
            result = Eval(struct("configuration",config,...
                "In",input));
            
            u1 = IP_CQEM_FD(input',config.dt,config.dz,fiber,config.f0,config.tol,0,1);
            u2 = IP_CQEM_FD(u1,config.dt,config.dz,fiber,config.f0,config.tol,0,1);
            u3 = IP_CQEM_FD(u2,config.dt,config.dz,fiber,config.f0,config.tol,0,1);
            
            testCase.assertEqual(result,transpose(u3));
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
            In1 = simulation.Input("In1");
            In2 = simulation.Input("In2");
            OC = simulation.Coupler("OC", component.Coupler(rho));
            
            [Out1,Out2] = [In1, In2] + OC;
            Eval1 = Out1.build(simulation.builder.Builder());
            Eval2 = Out2.build(simulation.builder.Builder());
            
            config = simulation.Configuration();
            input1 = rand_sech(config.nt, config.time)';
            input2 = fliplr(rand_sech(config.nt, config.time))';
            result1 = Eval1(struct("configuration",config,...
                "In1",input1,"In2",input2));
            result2 = Eval2(struct("configuration",config,...
                "In1",input1,"In2",input2));
            
            [ u1o,u2o ] = coupler(input1,input2,rho);
            testCase.assertEqual(result1,u1o,"RelTol",1e-6);
            testCase.assertEqual(result2,u2o,"RelTol",1e-6);
        end
        
        function testRecuurenceSwitch(testCase)
            createFixture(testCase);
            
            config = simulation.Configuration();
            
            fiber = component.Fiber();
            fiber.L = 0.001;
            fiber = fiber.cal_gamma(config.lambda0);
            
            SMF = simulation.Fiber("SMF", fiber);
            S = simulation.Switch("S","First_Run");
            
            feedback = simulation.Feedback();
            
            In = simulation.Input("In");
            u = simulation.Recurrence("u",In,"Run_until_stop");
            u1 = S(In,u) + SMF;
            t = feedback(u,u1);
            
            config.N_trip = 10;
            Eval = t.build(simulation.builder.Builder());
            input = rand_sech(config.nt,config.time)';
            result = Eval(struct("configuration",config,...
                "In",input,"Out_id",4));
            
            compare = input';
            for ii = 1:config.N_trip
                compare = IP_CQEM_FD(compare,config.dt,config.dz,fiber,config.f0,config.tol,0,0);
            end
            testCase.assertEqual(result,transpose(compare));
        end
        
        function testFigure8Evaluation(testCase)
            createFixture(testCase);
            
            % Structural Parameters to investigate
            landa_bw = 11; 
            smf5_L = 0.025;
            NOLM_L = 0.0015;
            PsatdBm = 23.5;
            smf4_L = 0.002;
            
            % used default configuration
            config = simulation.Configuration();
            
            % construct components
            coupler_NOLM = component.Coupler(0.3);
            coupler_output = component.Coupler(0.4);
            passive_fiber = component.Fiber("F10125").cal_gamma(config.lambda0);
            active_fiber = component.ActiveFiber("YDF10125", 37.5, PsatdBm,config.lambda0,80).cal_gamma(config.lambda0);
            filter = component.Filter(config.lambda0, landa_bw, 7);
            
            % assign components to models
            SMF1 = simulation.Fiber("SMF1", passive_fiber.set_L(0.0007));
            SMF2 = simulation.Fiber("SMF2", passive_fiber.set_L(NOLM_L-0.0007));
            SMF4 = simulation.Fiber("SMF4", passive_fiber.set_L(smf4_L));
            SMF5 = simulation.Fiber("SMF5", passive_fiber.set_L(smf5_L));
            AF = simulation.Fiber("AF", ...
                active_fiber.set_L(0.0008)...
                .set_RepeatFre(sum([NOLM_L,smf4_L,smf5_L,0.0008])));
            BPF = simulation.Filter("BPF", filter);
            OC = simulation.Coupler("OC", coupler_output);
            Coupler_NOLM = simulation.Coupler("Coupler_NOLM", coupler_NOLM);
            
            % define functional models
            S = simulation.Switch("S","First_Run");
            feedback = simulation.Feedback();
            show = simulation.Show(struct());
            
            % connect models
            % S(In,u)->BPF->SMF4->AMF1->SMF5->NOLM->OC->u
            In = simulation.Input("In");
            
            u = simulation.Recurrence("u",In,"Run_until_stop");
            
            u1 = S(In,u) + BPF + SMF4 + AF + SMF5;
            t = u1;
            [uf, ub] = [u1, 0] + Coupler_NOLM;
            ufo = uf + SMF1 + SMF2;
            ubo = ub + SMF2 + SMF1;
            [ur,ut] = [ubo, ufo] + Coupler_NOLM;
            
            [u2, uout] = [ut, 0] + OC;  % no way to show uout for now
            t = feedback(u,show(u2));
                  
            input = rand_sech(config.nt,config.time)';
            config.N_trip = 20;
            config.tol = 1e-7;
            Eval = t.build(simulation.builder.Builder());
            result = Eval(struct("configuration",config,...
                "In",input));
        end
    end
end

