function J=mucalib14_jac(p,mucal)
%
% This version is for simultaneous calibration of the 10-Be and 26-Al attenuation lengths, the
% erosion rate, and the muon production parameters for 10-Be and 26-Al.  
%
% Based on corejac1026.m from CRONUScalc: https://bitbucket.org/cronusearth/cronus-calc
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

% Make space for the J matrix.
J=zeros(mucal.nresids,mucal.npars);

% Compute the derivative with respect to p(1).
r0=mucalib14_fun(p,mucal);
p1=p;
p1(1)=p(1)*1.001;
r1=mucalib14_fun(p1,mucal);
J(:,1)=(r1-r0)/(0.001*p(1));

% Compute the derivative with respect to p(2).
p1=p;
p1(2)=p(2)*1.001;
r1=mucalib14_fun(p1,mucal);
J(:,2)=(r1-r0)/(0.001*p(2));

% Compute the derivative with respect to p(3).
p1=p;
p1(3)=p(3)*1.001;
r1=mucalib14_fun(p1,mucal);
J(:,3)=(r1-r0)/(0.001*p(3));

% Compute the derivative with respect to p(4).
p1=p;
p1(4)=p(4)*1.001;
r1=mucalib14_fun(p1,mucal);
J(:,4)=(r1-r0)/(0.001*p(4));
