% FastAnisotropicCurvaturePreservingSmoothingDemo
% ----------------------------------------------------------------------------------------------- %
%   Demo script of Fast Anisotropic Curvature Preserving Smoothing
% Remarks:
%   1.  This code is maintained in:
%       *   https://github.com/RoyiAvital/Fast-Anisotropic-Curvature-Preserving-Smoothing
%       *   "Add Matlab File Exchange Link".
%   2.  Prefixes:
%       -   't' - Tensor.
%       -   'm' - Matrix.
%       -   'v' - Vector.
%   3.  To run this demo just click "Run".
% Known Issues:
%   1.  Images must be in the range [0, 255].
%   2.  High amplitude values creates artifacts in bright areas in the
%       image.
% TODO:
%   1.  Move from images at the range [0, 255] to the range [0, 1].
%   2.  Add support for more advanced interpolation mwrhods for the "LIC"
%       pahse.
%   3.  Insert teh angle discretezation level as an input parameter.
%   4.  Add "numIterations" parameter.
%   5.  Add the option to use Gaussian Weight Window for the LIC phase.
%   6.  Add option to enable multi threaded operation (By choice).
%   7.  GPU mode of operation.
%   Release Notes:
%   -   1.0.000     27/10/2014  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

mInputImage = double(imread('Lena512.png'));

% Only single channel image are supported, for color image, run on each
% channel.
mInputImage = mInputImage(:, :, 1);

numRows = size(mInputImage, 1);
numCols = size(mInputImage, 2);

% Additive White Gaussian Noise parameters
mAwgnMean   = 0;
mAwgnStd    = 7;

mInputImageNoisy = mInputImage + ((mAwgnMean .* ones(numRows, numCols)) + (mAwgnStd .* randn(numRows, numCols)));

% Smoothing parameters
smoothingAmplitude  = 60;
sharpnessLevel      = 0.7;
anisotropyLevel     = 0.6;
gradientSmoothness  = 0.6;
tensorSmoothness    = 1.1;
stepSize            = 0.8;

% Running the algorithm
hTimerSatrt = tic();
[ mOutputImage ] = FastAnisotropicCurvaturePreservingSmoothing(mInputImage, ...
    smoothingAmplitude, sharpnessLevel, anisotropyLevel, gradientSmoothness, tensorSmoothness, stepSize);
runTime = toc(hTimerSatrt);

disp(['FastAnisotropicCurvaturePreservingSmoothing Run Time - ', num2str(runTime), ' [Sec]']);

figure();
imshow(mInputImage, [0, 255]);

figure();
imshow(mInputImageNoisy, [0, 255]);

figure();
imshow(mOutputImage, [0, 255]);



