classdef DFGTest < matlab.unittest.TestCase
    %DFGTEST Summary of this class goes here
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
        function testSingleFiber(testCase)
            % disp, char, + should all work in this test
            fiber = component.Fiber();
            SMF = simulation.Fiber("SMF", fiber)
            In = simulation.Input("In")
            Out = In + SMF;
            statements = Out.statements()
            testCase.verifyEqual(statements,...
                join(["t0 = In.Input()",...
                "t1 = SMF.Fiber[fiber=component.Fiber](t0)"],newline));
        end
        
        function testMultipleFiber(testCase)
            createFixture(testCase);
            fiber = component.Fiber();
            SMF1 = simulation.Fiber("SMF1", fiber);
            SMF2 = simulation.Fiber("SMF2", fiber);
            SMF3 = simulation.Fiber("SMF3", fiber);
            In = simulation.Input("In");
            Out = In + SMF1 + SMF2 + SMF3;
            statements = Out.statements()
            testCase.verifyEqual(statements,...
                join(["t0 = In.Input()",...
                "t1 = SMF1.Fiber[fiber=component.Fiber](t0)",...
                "t2 = SMF2.Fiber[fiber=component.Fiber](t1)",...
                "t3 = SMF3.Fiber[fiber=component.Fiber](t2)"],newline));
        end
        
        function testCoupler(testCase)
            % Test Coupler model
            createFixture(testCase);
            coupler = component.Coupler();
            In1 = simulation.Input("In1");
            In2 = simulation.Input("In2");
            OC = simulation.Coupler("OC", coupler);
            
            % cell() is overloaded here, the warning here can be neglected
            [Out1,Out2] = [In1, In2] + OC;
            
            statements1 = Out1.statements()
            statements2 = Out2.statements()
            testCase.verifyEqual(statements1,...
                join(["t0 = In1.Input()",...
                "t1 = In2.Input()",...
                "t2 = OC.CouplerOut1[coupler=component.Coupler](t0,t1)"],newline));
            testCase.verifyEqual(statements2,...
                join(["t0 = In1.Input()",...
                "t1 = In2.Input()",...
                "t3 = OC.CouplerOut2[coupler=component.Coupler](t0,t1)"],newline));
        end
        
        function testNOLM(testCase)
            % Test Creation of a NOLM model
            createFixture(testCase);
            coupler = component.Coupler();
            fiber = component.Fiber();
            SMF = simulation.Fiber("SMF1", fiber);
            In = simulation.Input("In");
            OC = simulation.Coupler("OC", coupler);
            [u1,u2] = [In, 0] + OC;
            u3 = u1 + SMF;
            u4 = u2 + SMF;
            [~, Out] = [u3, u4] + OC;
            statements = Out.statements()
            testCase.verifyEqual(statements,...
                join(["t0 = In.Input()",...
                "t1 = .Const[value=0]()",...
                "t2 = OC.CouplerOut1[coupler=component.Coupler](t0,t1)",...
                "t4 = SMF1.Fiber[fiber=component.Fiber](t2)",...
                "t3 = OC.CouplerOut2[coupler=component.Coupler](t0,t1)",...
                "t5 = SMF1.Fiber[fiber=component.Fiber](t3)",...
                "t7 = OC.CouplerOut2[coupler=component.Coupler](t4,t5)"],newline));
        end
        
        function testSwitch(testCase)
            createFixture(testCase);
            In1 = simulation.Input("In1");
            In2 = simulation.Input("In2");
            S = simulation.Switch("S","First_Run");
            Out = S(In1,In2);
            statements = Out.statements()
            testCase.verifyEqual(statements,...
                join(["t0 = In1.Input()",...
                "t1 = In2.Input()",...
                "t2 = S.Switch[condition=First_Run](t0,t1)"],newline));
        end
        
        function testShow(testCase)
            createFixture(testCase);
            In = simulation.Input("In");
            show = simulation.Show(struct());
            In = show(In); % It is fine to write like this
            statements1 = In.statements()
            testCase.verifyEqual(statements1,...
                join(["t0 = In.Input()",...
                "t1 = .Show(t0)",],newline));
            
            In = simulation.Input("In");
            Out = show(In); % It is better to write like this
            statements2 = Out.statements()
            testCase.verifyEqual(statements2,...
                join(["t2 = In.Input()",...
                "t3 = .Show(t2)",],newline));
            
            % In->SMF->Show->Out
            In = simulation.Input("In");
            fiber = component.Fiber();
            SMF = simulation.Fiber("SMF", fiber);
            SMF = simulation.Fiber("SMF", fiber);
            Out = show(In+SMF);
            statements3 = Out.statements()
            testCase.verifyEqual(statements3,...
                join(["t4 = In.Input()",...
                "t5 = SMF.Fiber[fiber=component.Fiber](t4)",...
                "t6 = .Show(t5)"],newline));
        end
        
        function testRecuurence(testCase)
            createFixture(testCase);
            
            coupler = component.Coupler();
            fiber = component.Fiber();
            
            SMF = simulation.Fiber("SMF", fiber);
            OC = simulation.Coupler("OC", coupler);
            S = simulation.Switch("S","First_Run");
            
            feedback = simulation.Feedback();
            
            In = simulation.Input("In");
            u = simulation.Recurrence("u",0,"Run_until_stop");
            u1 = S(In,u);
            [u2, ~] = [u1, 0] + OC;
            t = feedback(u,u2);
            
            % use the latest model to see the whole statements
            statements = t.statements()
            testCase.verifyEqual(statements,...
                join(["t0 = In.Input()",...
                "t1 = .Const[value=0]()",...
                "t2 = u.Recurrence[condition=Run_until_stop](t1)",...
                "t3 = S.Switch[condition=First_Run](t0,t2)",...
                "t4 = .Const[value=0]()",...
                "t5 = OC.CouplerOut1[coupler=component.Coupler](t3,t4)",...
                "t2 = u.Recurrence[condition=Run_until_stop](t5)"],newline));
        end
        
        function testFigure8(testCase)
            createFixture(testCase);
            
            coupler_NOLM = component.Coupler();
            coupler_output = component.Coupler();
            passive_fiber = component.Fiber();
            active_fiber = component.Fiber();
            filter = component.Filter();
            
            SMF1 = simulation.Fiber("SMF1", passive_fiber);
            SMF2 = simulation.Fiber("SMF2", passive_fiber);
            SMF3 = simulation.Fiber("SMF3", passive_fiber);
            SMF4 = simulation.Fiber("SMF4", passive_fiber);
            SMF5 = simulation.Fiber("SMF5", passive_fiber);
            AF = simulation.Fiber("AF", active_fiber);
            BPF = simulation.Filter("BPF", filter);
            
            OC = simulation.Coupler("OC", coupler_output);
            Coupler_NOLM = simulation.Coupler("Coupler_NOLM", coupler_NOLM);
            S = simulation.Switch("S","First_Run");
            
            feedback = simulation.Feedback();
            
            In = simulation.Input("In");
            u = simulation.Recurrence("u",0,"Run_until_stop");
            u1 = S(In,u);
            u3 = u1 + SMF1 + AF + SMF2;
            [u4, ~] = [u3, 0] + OC;
            u5 = u4 + SMF3 + BPF + SMF4;
            [u6, u7] = [u5, 0] + Coupler_NOLM;
            u8 = u6 + SMF5;
            u9 = u7 + SMF5;
            [~, u10] = [u8, u9] + Coupler_NOLM;
            t = feedback(u,u10);
            
            % use the latest model to see the whole statements
            statements = t.statements()
            testCase.verifyEqual(true, true)
        end
    end
    
end

