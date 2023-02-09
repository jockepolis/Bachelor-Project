function [N, val, err] = mcintegral(D, N)
% [N, val, err] = mcintegral(D, N);
% Function solving a multidimensional integral using Monte Carlo.
% Input: D -  The domain. Each row in the matrix D represents the range
%             of the corresponding variable. For example, a 3-dimensional
%             integral solved on the unit cube is given by
%             D = [0 1; 0 1; 0 1];
%             Default is a 6-dimensional integral in the range -5 to 5 in
%             each dimension
%        N -  A vector with the numbers of random numbers. For example, if
%             N = [100; 1000; 2000] three integral values will be computed,
%             based on 100, 1000 and 2000 random numbers, respectively.
%             The default is N = [100; 500; 1000; 5000; 10000; 1e5];
% Output: N   - The number of random numbers used (same as the input value
%               N)
%         val - the resulting integral value, a vector where val(i)
%               is the integral computed in N(i) random numbers
%         err - the error in the result (the standard deviation)

% Andreas Hellander, 2009. Updated by Stefan Pålsson 2018.

% Number of replications to base the error estimation on
M = 150;

% Default number of pseudorandom numbers
if (nargin < 2) || (isempty(N))
    N = [100; 500; 1000; 5000; 10000; 1e5];
end
N = N(:);  % Make N a column vector if it isn't

% Default 6D and range [-5 5] in each dimension
if (nargin == 0) || isempty(D)
    D = [-5 5;-5 5;-5 5;-5 5;-5 5; -5 5];
end

% The dimension and range
[~, columns] = size(D);
if columns ~= 2
    error('There must be 2 numbers (lower and upper limit) in each row in D');
end

dim = numel(D)/2;
disp(' ------');
disp(['A ',num2str(dim),'-dimensional integral is solved'])
disp(['for random numbers: ' num2str(N')]);

% Plot if 1D or 2D
if dim == 1
    fplot(@fnorm, [D(1) D(2)]);
end
if dim == 2
    NN = 100;
    x = linspace(D(1,1), D(1,2), NN);
    y = linspace(D(2,1), D(2,2), NN);
    [X,Y] = meshgrid(x,y);
    Z = fnorm([X(:)';Y(:)']);
    Z = reshape(Z,NN,NN);
    surf(X,Y,Z);
    shading interp;
end
drawnow;

% Calculate the integral
j=1;
val = zeros(length(N),1);
err = val;
for i=1:length(N)
    [val(j),err(j)]=mcint(@fnorm, D, N(i), M);
    j=j+1;
end

disp('Done!')
disp(' ------');

end % end main

% ----
% Internal functions
% ----

function y = fnorm(x)
% y = fnorm(x)
% Function defining the density of the d-dim standard normal distribution.
% Input: x - d*N-matrix. Each column of x is a point where the function 
% is evaluated.  
%

%Andreas Hellander, 2009.

   [d,~]=size(x);
   C = eye(d);
   y = 1/(2*pi)^(d/2)*exp(-0.5*sum(x.*(C*x),1));
end



function [val,err] = mcint(f,D,N,M) 
% [val,err] = mcint(@fun,D,N,M)  
% Function that solves an integral using Monte Carlo.    
% Input:  fun - function defining the integrand. 
%               fun will be evaluated in dim*N points, where dim is the
%               dimension of the problem. 
%         D   - the domain. Each row in the matrix D represents the range 
%               of the corresponding variable. For example, a 3-dimensional
%               integral solved on the unit cube would be given by             
%               D = [0 1; 0 1; 0 1];
%         N   - the number of points in each realization.
%         M   - the number of repetitions used for error estimation. 
%               (Recommendation, M = 30+).
%               Total number of points used is thus M*N.  
%
% Output:  val - the resulting integral value
%          err - the error in the result (the standard deviation)

% Andreas Hellander, 2009.
    
    V   = cumprod(D(:,2)-D(:,1));
    dim = numel(V);
    V   = V(end);
    
    r=zeros(dim,N);
    
    for j=1:M
        r(:,:) = repmat(D(:,1),1,N)+rand(dim,N).*repmat(D(:,2)-D(:,1),1,N);
        I(j)=V*mean(f(r));
    end
    
    val = mean(I);
    err = std(I);
end
