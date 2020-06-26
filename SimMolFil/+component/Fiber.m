classdef Fiber
    %FIBER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        core_radius     % um
        Amod            % Mode Field Area (um^2)
        n2              % Kerr coefficient (10^-16*cm^2/W)
        gamma           % W^-1 * km^-1
        alpha = 0    	% Atenuation coefficient (km^-1)
        betaw           % beta coefficients (ps^n/nm)
        L = 0   	% fiber length (km)
        raman = 0   	% 0 = exclude raman effect
        ssp = 0        	% 0 = disable self sharpen effect
        type            % fiber type, to load predefined fiber specs
    end
    
    methods
        function obj = Fiber(varargin)
            %FIBER Construct an instance of this class
            import component.*
            switch(nargin)
                case 0
                    core_radius = 3.3;
                    n2 = 2.6;
                    betaw = [0 0 24.5864 26.1949e-3];
                    obj.type = "PM980";
                case 1
                    % Fiber(type)
                    if isa(varargin{1},"string")
                        switch (varargin{1})
                            case "PM980"
                                core_radius = 3.3;
                                n2 = 2.6;
                                betaw = [0 0 24.5864 26.1949e-3];
                                obj.type = "PM980";
                            case "PM1025"
                                core_radius = 5;
                                n2 = 2.6;
                                betaw = [0 0 20.8634 34.9621e-3];
                                obj.type = "PM1025";
                            case "F10125"
                                core_radius = 5.25;
                                n2 = 2.6;
                                betaw = [0 0 20.1814 36.8057e-3];
                                obj.type = "F10125";
                            otherwise
                                error("Fiber: not implemented fiber type %s\n",varargin{1});
                        end
                    else
                        error("Fiber: not implemented argument types\n");
                    end
                case 3
                    % Fiber(core_radius, n2, betaw)
                    core_radius = varargin{1};
                    n2 = varargin{2};
                    betaw = varargin{3};
                otherwise
                    error("Fiber: not implemented argument types\n");
            end
            obj.core_radius = core_radius;
            obj.Amod = pi*obj.core_radius^2;
            obj.n2 = n2;
            obj.betaw = betaw;
        end
        
        function obj = cal_gamma(obj,lambda0)
            obj.gamma = 2*pi*obj.n2/lambda0/obj.Amod*1e4;
        end
        
        function obj = set_L(obj,L)
            obj.L = L;
        end
    end
end

