classdef orthDecomp < handle
	%% properties
    properties (Access=public)
        orthBasis;  % Orthonormalization
        rawBasis;   % Input
        interval;   % Interval
        appFcn;     % Approaching Function
        inner;      % Inner Product
        fcn;        % Function
        err;        % Error
    end
    %% constructor
    methods
        function obj = orthDecomp(rawBasis,fcn,interval)
            % default properties
            if nargin<3 || length(interval)<2
                interval = [0,1];
            end
            if nargin<2
                fcn = @(x)sqrt(x);
            end
            if nargin<1
                syms x;
                rawBasis = [1,x,x^2,x^3];
            end
            % assign values
            obj.rawBasis = rawBasis(:);
            obj.interval = interval;
            obj.inner    = @(fx,fy)int(fx*fy,interval(1),interval(2));
            obj.fcn      = fcn;
        end
        % Input: [sym,anonymousFcn,coordinate]
    end
    %% methods
    methods (Static=true)
        % Orthonormalize
        function orthBasis = orthnorm(obj)
            % initial
            orthBasis = obj.rawBasis;
            inner     = obj.inner;
            len       = length(orthBasis);
            coe       = zeros(1,len);
            % main
            for i=1:len
                coe(i) = 1;
                for j=i-1:-1:1
                    coe(j) = -inner(orthBasis(i),orthBasis(j));
                end
                orthBasis(i) = coe*orthBasis;
                orthBasis(i) = orthBasis(i)/sqrt(inner(orthBasis(i),orthBasis(i)));
            end
            % assignment
            obj.orthBasis = orthBasis;
        end
        % calculate error
        function err = calcErr(obj)
            % initial
            orthBasis = obj.orthBasis;
            inner     = obj.inner;
            len       = length(orthBasis);
            tmp       = zeros(len);
            for i=1:len
                for j=1:i
                    tmp(i,j) = inner(orthBasis(i),orthBasis(j));
                end
            end
            % assignment
            err     = norm(tmp-eye(len));
            obj.err = err;
        end
        % approach
        function fcn = approach(obj)
            % initial
            basis = obj.orthBasis;
            inner = obj.inner;
            fcn   = obj.fcn;
            % main
            len = length(basis);
            coe = zeros(1,len);
            for i=1:len
                coe(i) = inner(fcn,basis(i));
            end
            fcn        = coe*basis;
            obj.appFcn = fcn;
        end
        % visualization
        function err = visualize(obj)
            % initial
            interval = obj.interval;
            appFcn   = matlabFunction(obj.appFcn);
            fcn      = obj.fcn;
            % main
            x = linspace(interval(1),interval(2),1000);
            figure('FileName','orthDecomp','NumberTitle','off');
            plot(x,fcn(x),'-',x,appFcn(x),'.');
            legend('Original','Result');
            % calculate error
            err = sum(abs(fcn(x)-appFcn(x)))/sum(abs(fcn(x)));
            disp(['Error: ',num2str(err*100),'%']);
        end
    end
end













