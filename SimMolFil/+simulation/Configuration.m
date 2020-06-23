classdef Configuration
    %CONFIGURATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nt
        time
        tol
        N_trip
    end
    
    methods
        function obj = Configuration(varargin)
            %CONFIGURATION Construct an instance of Configuration
            obj.nt = 2^10;
            obj.time = 100;
            obj.tol = 1e-7;
            obj.N_trip = 300;
            switch(length(varargin))
                case 0
                    % pass
                case 1
                    obj.nt = varargin{1};
                case 2
                    obj.nt = varargin{1};
                    obj.time = varargin{2};
                case 3
                    obj.nt = varargin{1};
                    obj.time = varargin{2};
                    obj.tol = varargin{3};
                case 4
                    obj.nt = varargin{1};
                    obj.time = varargin{2};
                    obj.tol = varargin{3};
                    obj.N_trip = varargin{4};
            end
        end
    end
end

