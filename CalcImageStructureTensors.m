function [ tT ] = CalcImageStructureTensors( mInputImage, sharpnessLevel, anisotropyLevel, gradientSmoothness, tensorSmoothness )
% ----------------------------------------------------------------------------------------------- %
% [ tT ] = CalcImageStructureTensors( mInputImage, sharpnessLevel, anisotropyLevel, gradientSmoothness, tensorSmoothness )
%   Calculates the Image Structure Tensor
% Input:
%   - mInputImage           -   Input image.
%                               Matrix, 1 Channels, Floating Point, [0, 1]
%   - sharpnessLevel        -   Local Window Radius.
%                               Scalar, Floating Point, {1, 2, ..., 10}.
%   - anisotropyLevel           -   Local Window Gaussian Kernel STD.
%                               Scalar, Floating Point [0.1, 20].
%   - gradientSmoothness       -   Search Window Radius.
%                               Scalar, Floating Point, {1, 2, ..., 10}.
%   - tensorSmoothness            -   Weights STD Factor.
%                               Scalar, Floating Point [0.1, 20].
% Output:
%   - mOutputImage          -   Input image.
%                               Matrix, 1 Channels, Floating Point, [0, 1]
% Remarks:
%   1.  Prefixes:
%       -   't' - Tensor.
%       -   'm' - Matrix.
%       -   'v' - Vector.
%   2.  Cl
% TODO:
%   1.  aa
%   Release Notes:
%   -   1.0.000     27/10/2014  Or Yair
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

tTensors = zeros([size(mInputImage), 3]);

sharpnessLevel = max(sharpnessLevel, 1e-5);
	
power1 = sharpnessLevel / 2;
power2 = power1    / (1e-7 + 1 - anisotropyLevel);

%%
% Gaussian Blurring Kernel
vGradientSmoother    = fspecial('gaussian', [ceil(3 * gradientSmoothness),...
	                           ceil(3 * gradientSmoothness)], gradientSmoothness);
mInputImageSmoothed = imfilter(mInputImage, vGradientSmoother, 'same', 'replicate');

% Notmalize ?? %

%%
vDerivativeKernel = 0.5 * [-1, 0, 1];
mIx = imfilter(mInputImageSmoothed, vDerivativeKernel,  'replicate', 'same');
mIy = imfilter(mInputImageSmoothed, vDerivativeKernel', 'replicate', 'same');

mIxx = mIx .* mIx;
mIyy = mIy .* mIy;
mIxy = mIx .* mIy; % Ixy = Iyx

%%
% Gaussian Blurring Kernel
vTensorSmoother  = fspecial('gaussian', [ceil(3 * tensorSmoothness), ceil(3 * tensorSmoothness)], tensorSmoothness);
tTensors(:,:,1) = imfilter(mIxx, vTensorSmoother, 'same', 'replicate');
tTensors(:,:,2) = imfilter(mIxy, vTensorSmoother, 'same', 'replicate');
tTensors(:,:,3) = imfilter(mIyy, vTensorSmoother, 'same', 'replicate');

%%
% Using Jacobi Method
Tau = ( tTensors(:,:,3) - tTensors(:,:,1) ) ./ ( 2 * tTensors(:, :, 2) );
t   = sign(Tau) ./ ( abs(Tau) + sqrt( 1 + (Tau .* Tau) ) );
c   = 1 ./ sqrt( 1 + (t .* t) );
s   = c .* t;

tThetaA = cat(3, c, -s); 
tThetaB = cat(3, s, c);  

mLambdaA = tTensors(:, :, 1) - ( t .* tTensors(:, :, 2) ); 
mLambdaB = tTensors(:, :, 3) + ( t .* tTensors(:, :, 2) ); 

mMinLambdaMask = mLambdaA < mLambdaB;
tMinLambdaMask = cat(3, mMinLambdaMask, mMinLambdaMask);

mLambda1 = mMinLambdaMask .* mLambdaA + (1 - mMinLambdaMask) .* mLambdaB;
mLambda2 = mMinLambdaMask .* mLambdaB + (1 - mMinLambdaMask) .* mLambdaA;

tTheta1 = tMinLambdaMask .* tThetaA + (1 - tMinLambdaMask) .* tThetaB;
tTheta2 = tMinLambdaMask .* tThetaB + (1 - tMinLambdaMask) .* tThetaA;
 
% Step 6 - T Matrix
% n1 = (1 + mLambda1 + mLambda2) .^ (-power1 / 2);
% n2 = (1 + mLambda1 + mLambda2) .^ (-power2 / 2);

n1 = (1 + mLambda1 + mLambda2) .^ (-power1);
n2 = (1 + mLambda1 + mLambda2) .^ (-power2);

tT(:, :, 1) = n1 .* (tTheta1(:, :, 1) .* tTheta1(:, :, 1)) + ...
                    n2 .* (tTheta2(:, :, 1) .* tTheta2(:, :, 1));
		 
tT(:, :, 2) = n1 .* (tTheta1(:, :, 1) .* tTheta1(:, :, 2)) + ...
                    n2 .* (tTheta2(:, :, 1) .* tTheta2(:, :, 2));
		 
tT(:, :, 3) = n1 .* (tTheta1(:, :, 2) .* tTheta1(:, :, 2)) + ...
                    n2 .* (tTheta2(:, :, 2) .* tTheta2(:, :, 2));
                
end
