function r=mucalib_fun(p,mucal)
%
% This version is for simultaneous calibration of the 10-Be and 26-Al attenuation lengths, the
% erosion rate, and the muon production parameters for 10-Be and 26-Al.  
%
% Based on corefun1026.m from CRONUScalc: https://bitbucket.org/cronusearth/cronus-calc
% Jakob Heyman - 2016-2019 (jakob.heyman@gu.se)
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% The parameters are:
%   p(1)             erosion rate
%   p(2)             attenuation length for 10-Be
%   p(3)             attenuation length for 26-Al
%   p(4)             fstar10 (scaled by 1.0e-3)
%   p(5)             sigma010 (scaled by 1.0e-30)
%   p(6)             fstar26 (scaled by 1.0e-3)
%   p(7)             sigma026 (scaled by 1.0e-30)

% Deal with the case that any parameter is negative.
if any(p<0)
    r=1.0e30*ones(mucal.nresids,1);
    return
end

tvzv = p(1).*1E-3 .* mucal.tv; % depth vector for tv (g/cm2)
tvzm10 = bsxfun(@plus,tvzv',mucal.Nz10'); % depth matrix for tv10 (one col per depth sample)
tvzm26 = bsxfun(@plus,tvzv',mucal.Nz26'); % depth matrix for tv26 (one col per depth sample)

% sp production
Psp10rm = exp(-tvzm10./p(2)); % spallation 10 depth prod ratio
Psp26rm = exp(-tvzm26./p(3)); % spallation 26 depth prod ratio
Psp10m = bsxfun(@times,Psp10rm,mucal.Psp10'); % spallation 10Be production matrix
Psp26m = bsxfun(@times,Psp26rm,mucal.Psp26'); % spallation 26Al production matrix

% 10Be muon prod
Pfast10m = interp1(mucal.dz,mucal.Pfast10d,tvzm10,'pchip') .* p(5).*1E-31; % calc Pfast10 prod matr
Pneg10m = interp1(mucal.dz,mucal.Pneg10d,tvzm10,'pchip') .* p(4).*1E-3; % calc Pneg10 prod matr

% 26Al muon prod
Pfast26m = interp1(mucal.dz,mucal.Pfast26d,tvzm26,'pchip') .* p(7).*1E-30; % calc Pfast26 prod matr
Pneg26m = interp1(mucal.dz,mucal.Pneg26d,tvzm26,'pchip') .* p(6).*1E-2; % calc Pneg26 prod matr

% full production plus decay
Pfull10m = Psp10m + Pfast10m + Pneg10m; % full P10 matrix
Pfull26m = Psp26m + Pfast26m + Pneg26m; % full P26 matrix
Pfull10lm = bsxfun(@times,Pfull10m,mucal.dcf10'); % full P10 including decay matrix
Pfull26lm = bsxfun(@times,Pfull26m,mucal.dcf26'); % full P26 including decay matrix

% calculate N by integration
N10v = trapz(mucal.tv',Pfull10lm); % calculated N10
N26v = trapz(mucal.tv',Pfull26lm); % calculated N26

% calculate residuals
r10 = (N10v-mucal.N10')./mucal.N10unc';
r26 = (N26v-mucal.N26')./mucal.N26unc';
r = [r10 r26]';
