classdef Filter
    %FILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lambda_c       	% filter centre wavelength (nm)
        lambda_bw      	% filter bandwidth (nm)
        order          	% filter order
        
        fc              % filter centre frequency (THz)
        f3dB            % filter 3dB bandwidth(THz)
    end
    
    methods
        function obj = Filter(lambda_c, lambda_bw, order)
            %FILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.lambda_c = lambda_c;
            obj.lambda_bw = lambda_bw;
            obj.order = order;
            
            obj.fc = simulation.Constants.c/lambda_c;
            obj.f3dB = simulation.Constants.c/(lambda_c)^2*lambda_bw;
        end
        
        function uo = filter_gauss(obj, ui, f0, df)
            %FILTER_GAUSS Apply gaussian filterting to light field
            % ui: input field amplitude (column vector)
            % fo: central pulse frequency (THz)
            % df: frequencies separation (THz)
            uo = filter_gauss(transpose(ui),obj.f3dB,obj.fc,obj.order,f0,df) ;
        end
    end
    
end

