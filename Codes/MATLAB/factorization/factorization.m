classdef factorization < handle
	properties (Access=public)
		method;
	end
	%% constructor
	methods
		function obj = factorization()
			obj.method = {'lu', 'qr', 'doolittle', 'cholesky', 'crout'};
		end
	end
	%% methods
	methods (Static=true)
        % lu
		function [L,U] = lu(matrix)
			[L,U] = lu(matrix);
        end
        % qr
        function [Q,R] = qr(matrix)
            [Q,R] = qr(matrix);
        end
        % doolittle
        function [L, U] = doolittle(matrix)
            L = doolittle(matrix);
            if nargout>1
                U = tril(L, -1) + eye(size(L));
                L = triu(L);
            end
        end
        % cholesky
        function [L,R] = cholesky(matrix)
            L = chol(matrix, 'lower');
            if nargout>1
                R = L';
            end
        end
        % crout
        function [L,U] = crout(matrix)
            [L,U] = crout(diag(matrix), diag(matrix,1), diag(matrix,-1));
        end
	end
end