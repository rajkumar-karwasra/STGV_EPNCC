function result = RollingGuidanceFilter(I, sigma_s, sigma_r, ...
     iteration,GaussianPrecision)

if ~exist('iteration','var')
    iteration = 4;
end

if ~exist('sigma_s','var')
    sigma_s = 3;
end

if ~exist('sigma_r','var')
    sigma_r = 0.1;
end

if ~exist('GaussianPrecision','var')
    GaussianPrecision = 0.08;
end
 
[a b c] = size(I);

GaussianWindow = WindowBlock(sigma_s, GaussianPrecision);

N = ( size(GaussianWindow,1) -1)/2;

J_plus = I;

%Guide = imguidedfilter(I);
Guide = imguidedfilter(I,...
    'NeighborhoodSize',[9 9],...
    'DegreeOfSmoothing',0.0005);

I = ExpandBorder(I,N);

J = ExpandBorder(Guide,N);


for l = 1:iteration
    for k = 1 : c
        for j = 1+N : b+N
            for i = 1+N : a+N
                RangeWeight = exp(-(J(i-N:i+N,j-N:j+N,k)-J(i,j,k)).^2/(2*sigma_r^2));

LocalGuide = J(i-N:i+N,j-N:j+N,k);

%[Gx,Gy] = gradient(LocalGuide);
[Gx,Gy] = imgradientxy(LocalGuide,'sobel');
%EdgeWeight = exp(-2*sqrt(Gx.^2 + Gy.^2));
edge = sqrt(Gx.^2 + Gy.^2);

edge = edge/(max(edge(:))+eps);

EdgeWeight = exp(-0.6*edge);
H = GaussianWindow .* RangeWeight .* EdgeWeight;
               K_p = sum(H(:));

if K_p < eps
    K_p = eps;
end

J_plus(i-N,j-N,k) = sum(sum(H .* I(i-N:i+N,j-N:j+N,k))) / K_p;

            end
        end
    end

    J = ExpandBorder(J_plus, N);

end   % <-- closes the for l = 1:iteration loop

Residual = I(N+1:end-N,N+1:end-N,:) - J_plus;

grayR = rgb2gray(J_plus);

[Gx,Gy] = imgradientxy(grayR);

edgeMap = sqrt(Gx.^2 + Gy.^2);

edgeMap = edgeMap ./ (max(edgeMap(:))+eps);
Weight = exp(-2 * edgeMap);

%Weight = exp(-2*edgeMap);

Weight = repmat(Weight,[1 1 size(I,3)]);

noiseLevel = std(Residual(:));

alpha = max(0.05,min(0.15,0.5*noiseLevel));

result = J_plus + alpha .* Weight .* Residual;

result = max(min(result,1),0);

end    % <-- closes RollingGuidanceFilter function

%==================================================================

function GaussianWindow = WindowBlock(sigma_s, GaussianPrecision)
%function GaussianWindow = WindowBlock(sigma_s, GaussianPrecision)
%product a Block of 2D Gaussian.
%input    sigma_s: the Gaussian standard deviation; GaussianPrecision: contorl the size of Block.
%output   GaussianWindow: a square matrix of 2D Gaussian.

%right below
pq = bsxfun(@plus, ([0:sigma_s*3].^2)', [0:sigma_s*3].^2); 

% gaussian distribution
pqrb = exp(-pq/2/sigma_s^2); 

% element that is less than GaussianPrecision ar equal zero
pqrb = pqrb .* (pqrb>GaussianPrecision); 

% remove all zero column
pqrb(:, pqrb(1,:)==0) = []; 

% remove all zero row
pqrb(pqrb(:,1)==0, :) = []; 

%left below
pqlb = fliplr(pqrb); 

%right upper
pqru = flipud(pqrb); 

%left upper
pqlu = fliplr(pqru); 

GaussianWindow = [pqlu(:, 1:end-1)     pqru;
                  pqlb(2:end, 1:end-1) pqrb(2:end, :)];
  
end

%==================================================================

function imageExpand = ExpandBorder(image, N)
%function imageExpand = ExpandBorder(image, N)
%is On the edge of the pixels extend outward N pixels.
%input    image: the origic image; N: the number of expansion need.
%output   imageExpand: generate a new image

imageExpand = [repmat(image(1,:,:), [N,1,1]) ; image ; repmat(image(end,:,:), [N,1,1])];
imageExpand = [repmat(imageExpand(:,1,:), [1,N,1])  imageExpand  repmat(imageExpand(:,end,:), [1,N,1])];

end
