function J=mucalib_jac(p,mucal)
%
% This version is for simultaneous calibration of the 10-Be and 26-Al attenuation lengths, the
% erosion rate, and the muon production parameters for 10-Be and 26-Al.  
%
% Based on corejac1026.m from CRONUScalc: https://bitbucket.org/cronusearth/cronus-calc
% Jakob Heyman - 2016-2018 (jakob.heyman@gu.se)
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

% Make space for the J matrix.
J=zeros(mucal.nresids,mucal.npars);

% Compute the derivative with respect to p(1).
r0=mucalib_fun(p,mucal);
p1=p;
p1(1)=p(1)*1.001;
r1=mucalib_fun(p1,mucal);
J(:,1)=(r1-r0)/(0.001*p(1));

% Compute the derivative with respect to p(2).
p1=p;
p1(2)=p(2)*1.001;
r1=mucalib_fun(p1,mucal);
J(:,2)=(r1-r0)/(0.001*p(2));

% Compute the derivative with respect to p(3).
p1=p;
p1(3)=p(3)*1.001;
r1=mucalib_fun(p1,mucal);
J(:,3)=(r1-r0)/(0.001*p(3));

% Compute the derivative with respect to p(4).
p1=p;
p1(4)=p(4)*1.001;
r1=mucalib_fun(p1,mucal);
J(:,4)=(r1-r0)/(0.001*p(4));

% Compute the derivative with respect to p(5).
p1=p;
p1(5)=p(5)*1.001;
r1=mucalib_fun(p1,mucal);
J(:,5)=(r1-r0)/(0.001*p(5));

% Compute the derivative with respect to p(6).
p1=p;
p1(6)=p(6)*1.001;
r1=mucalib_fun(p1,mucal);
J(:,6)=(r1-r0)/(0.001*p(6));

% Compute the derivative with respect to p(7).
p1=p;
p1(7)=p(7)*1.001;
r1=mucalib_fun(p1,mucal);
J(:,7)=(r1-r0)/(0.001*p(7));
