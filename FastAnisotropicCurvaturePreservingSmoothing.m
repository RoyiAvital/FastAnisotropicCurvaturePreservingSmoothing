function [ mOutputImage ] = FastAnisotropicCurvaturePreservingSmoothing( mInputImage, ...
    smoothingAmplitude, sharpnessLevel, anisotropyLevel, gradientSmoothness, tensorSmoothness, stepSize )
% ----------------------------------------------------------------------------------------------- %
% [ mOutputImage ] = FastAnisotropicCurvaturePreservingSmoothing( mInputImage, ...
%    smoothingAmplitude, spatialPrecision, gradientSmoothness, tensorSmoothness )
%   Applies the Fast Anisotropic Curvature Preserving Smoothing on an Input Image
% Input:
%   - mInputImage           -   Input image.
%                               Matrix, 1 Channels, Floating Point, [0, 1]
%   - smoothingAmplitude        -   Local Window Radius.
%                               Scalar, Floating Point, {1, 2, ..., 10}.
%   - sharpnessLevel           -   Local Window Gaussian Kernel STD.
%                               Scalar, Floating Point [0.1, 20].
%   - anisotropyLevel       -   Search Window Radius.
%                               Scalar, Floating Point, {1, 2, ..., 10}.
%   - gradientSmoothness            -   Weights STD Factor.
%                               Scalar, Floating Point [0.1, 20].
%   - tensorSmoothness            -   Weights STD Factor.
%                               Scalar, Floating Point [0.1, 20].
%   - stepSize            -   Weights STD Factor.
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

ANGLE_TO_RAD_FACTOR = pi / 180;

gaussianKernelPrecision = 2;
vAngles     = [0, 30, 60, 90, 120, 150] * ANGLE_TO_RAD_FACTOR;
numAngles   = length(vAngles);

tTensors     = CalcImageStructureTensors(mInputImage, sharpnessLevel, anisotropyLevel, gradientSmoothness, tensorSmoothness);

mOutputImage = zeros( size(mInputImage) );


for fieldAngle = vAngles
	tTensorField(:, :, 1) = (tTensors(:, :, 1) * cos(fieldAngle)) + (tTensors(:, :, 2) * sin(fieldAngle));
	tTensorField(:, :, 2) = (tTensors(:, :, 2) * cos(fieldAngle)) + (tTensors(:, :, 3) * sin(fieldAngle));
	
	tTensorField(:, :, 3)   = sqrt((tTensorField(:, :, 1) .* tTensorField(:, :, 1)) + (tTensorField(:, :, 2) .* tTensorField(:, :, 2)));
	dln                     = stepSize ./ tTensorField(:, :, 3);
	
	tTensorField(:, :, 1) = tTensorField(:, :, 1) .* dln;
	tTensorField(:, :, 2) = tTensorField(:, :, 2) .* dln;
	
	mFsigma  = tTensorField(:, :, 3) * (sqrt(2 * smoothingAmplitude));
	mLength = gaussianKernelPrecision * mFsigma;
		
    mOutputImage = mOutputImage + ApplyImageLineIntegralConvolution(mInputImage, tTensorField, mLength, stepSize);
    
end

mOutputImage = mOutputImage / numAngles;


end

