classdef Configuration
    %CONFIGURATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nt = 2^12       % number of points
        time = 120      % width of time window (ps)
        lambda0 = 1030  % pulse central wavelength (nm)
        tol = 1e-7      % tolerance
        N_trip = 300    % number of total runtrips
        
        f0          	% pulse central frequency (THz)
        
        dz = 0.00001  	% Initial longitudinal step (km)
        dt              % time step (ps)
        df              % frequencies separation (THz)
        
        t               % time vector (ps)
        f               % frequencies vector (THz)
        w               % angular frequencies vector (THz)
        lambda          % lambdas vector (nm)
       
    end
    
    methods
        function obj = Configuration(varargin)
            %CONFIGURATION Construct an instance of Configuration
            for ii = 1:length(varargin)
                switch(ii)
                    case 1
                        obj.nt = varargin{ii};
                    case 2
                        obj.time = varargin{ii};
                    case 3
                        obj.lambda0 = varargin{ii};
                    case 4
                        obj.tol = varargin{ii};
                    case 5
                        obj.N_trip = varargin{ii};
                end
            end
            obj.dt = obj.time/obj.nt;
            obj.df = 1/obj.time;
            obj.f0 = simulation.Constants.c/obj.lambda0;
            
            c = simulation.Constants.c;
            obj.t = -obj.time/2:obj.dt:(obj.time/2-obj.dt);
            obj.f=-(obj.nt/2)*obj.df:obj.df:(obj.nt/2-1)*obj.df; 
            obj.lambda = c./(obj.f + c/obj.lambda0); 
            obj.w = 2*pi*obj.f;
        end
    end
end

