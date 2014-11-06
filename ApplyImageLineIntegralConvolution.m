function [ mOutputImage ] = ApplyImageLineIntegralConvolution(mInputImage, tW, mLength, stepSize)
% ----------------------------------------------------------------------------------------------- %
% [ mOutputImage ] = ApplyImageLineIntegralConvolution(mInputImage, tW, mLength, stepSize)
%   Applies the Line Integral Convolution (LIC) on an image
% Input:
%   - mInputImage   -   Input image.
%                       Matrix, 1 Channels, Floating Point, [0, 1]
%   - tW            -   Vector Field Tensor.
%                       Tensor, Floating Point, [].
%   - mLength       -   ****
%                       Scalar, Floating Point [].
%   - stepSize      -   Discrete Integration Step Size.
%                       Scalar, Floating Point, [].
% Output:
%   - mOutputImage  -   Input image.
%                       Matrix, 1 Channels, Floating Point, [0, 1]
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

mOutputImage = zeros( size(mInputImage) );

	for x = 1 : size(mInputImage, 2)
		for y = 1 : size(mInputImage, 1)
			val     = 0;
			l       = 0;
			counter = 0;
			X       = x;
			Y       = y;
			while l <= mLength(y,x)
            
				cx = round(X);
				cy = round(Y);
				if (cx > size(tW, 2) || cy > size(tW, 1) || cx <= 0 || cy <= 0)
					break;
				end
				
				val = val + mInputImage(cy, cx);
				
				X = X + tW(cy, cx, 1);
				Y = Y + tW(cy, cx, 2);
				l = l + stepSize;
				
				counter = counter + 1;
            end
			
			mOutputImage(y, x) = mOutputImage(y, x) + val/counter;
		end
	end
end