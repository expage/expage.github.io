function r=mucalib14_fun(p,mucal)
%
% This version is for simultaneous calibration of the 10-Be and 26-Al attenuation lengths, the
% erosion rate, and the muon production parameters for 10-Be and 26-Al.  
%
% Based on corefun1026.m from CRONUScalc: https://bitbucket.org/cronusearth/cronus-calc
% Jakob Heyman - 2022 (jakob.heyman@gu.se)
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% The parameters are:
%   p(1)             erosion rate
%   p(2)             attenuation length for 14-C
%   p(3)             fstar14 (scaled by 1.0e-1)
%   p(4)             sigma014 (scaled by 1.0e-31)

% Deal with the case that any parameter is negative.
if any(p<0)
    r=1.0e30*ones(mucal.nresids,1);
    return
end

tvzv = p(1).*1E-3 .* mucal.tv; % depth vector for tv (g/cm2)
tvzm14 = bsxfun(@plus,tvzv',mucal.Nz14'); % depth matrix for tv14 (one col per depth sample)

% sp production
Psp14rm = exp(-tvzm14./p(2)); % spallation 14 depth prod ratio
Psp14m = bsxfun(@times,Psp14rm,mucal.Psp14'); % spallation 14C production matrix

% 14C muon prod
Pfast14m = interp1(mucal.dz,mucal.Pfast14d,tvzm14,'pchip') .* p(4).*1E-31; % calc Pfast14 prod matr
Pneg14m = interp1(mucal.dz,mucal.Pneg14d,tvzm14,'pchip') .* p(3).*1E-1; % calc Pneg14 prod matr

% full production plus decay
Pfull14m = Psp14m + Pfast14m + Pneg14m; % full P14 matrix
Pfull14lm = bsxfun(@times,Pfull14m,mucal.dcf14'); % full P14 including decay matrix

% calculate N by integration
N14v = trapz(mucal.tv',Pfull14lm); % calculated N14

% calculate residuals
r14 = (N14v-mucal.N14')./mucal.N14unc';
r = [r14]';
