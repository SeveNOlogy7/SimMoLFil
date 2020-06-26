classdef ActiveFiber < component.Fiber
    %ACTIVEFIBER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gssdB
        PsatdBm
        
        lambda_c
        lambda_bw
        fbw
        fc
        
        RepeatFre
    end
    
    methods
        function obj = ActiveFiber(varargin)
            %ACTIVEFIBER Construct an instance of this class
            %   Detailed explanation goes here
            import component.*
            import simulation.*
            switch(nargin)
                case 0
                    sup_varargin = {"F10125"};
                    gssdB = 37.5;
                    PsatdBm = 23.5;
                    lambda_c = 1030;
                    lambda_bw = 80;
                case 1
                    % Fiber(type)
                    if isa(varargin{1},"string")
                        switch (varargin{1})
                            case "YDF10125"
                                sup_varargin = {"F10125"};
                                gssdB = 37.5;
                                PsatdBm = 23.5;
                                lambda_c = 1030;
                                lambda_bw = 80;
                            otherwise
                                error("Fiber: not implemented fiber type %s\n",varargin{1});
                        end
                    else
                        error("Fiber: not implemented argument types\n");
                    end
                case 5
                    % Fiber(type, gssdB, PsatdBm, lambda_c, lambda_bw)
                    if isa(varargin{1},"string")
                        switch (varargin{1})
                            case "YDF10125"
                                sup_varargin = {"F10125"};
                                gssdB = varargin{2};
                                PsatdBm = varargin{3};
                                lambda_c = varargin{4};
                                lambda_bw = varargin{5};
                            otherwise
                                error("Fiber: not implemented fiber type %s\n",varargin{1});
                        end
                    else
                        error("Fiber: not implemented argument types\n");
                    end
                case 7
                    % Fiber(core_radius, n2, betaw, gssdB, PsatdBm, lambda_c, lambda_bw)
                    sup_varargin = varargin{1:3};
                    gssdB = varargin{4};
                    PsatdBm = varargin{5};
                    lambda_c = varargin{6};
                    lambda_bw = varargin{7};
                otherwise
                    error("Fiber: not implemented argument types\n");
            end
            obj = obj@component.Fiber(sup_varargin{:});
            obj.gssdB = gssdB;
            obj.PsatdBm = PsatdBm;
            obj.lambda_c = lambda_c;
            obj.lambda_bw = lambda_bw;
            
            obj.fbw = Constants.c/(obj.lambda_c)^2*obj.lambda_bw;
            obj.fc = Constants.c/obj.lambda_c;
        end
        
        function obj = set_RepeatFre(obj, cavity_length)
            obj.RepeatFre = simulation.Constants.c/2/cavity_length;
        end
    end
end

