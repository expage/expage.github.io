function LSDspal_out = LSDspal(h,Rc,SPhi,w,nucl10,nucl26,consts)

% Implements the Lifton Sato Dunai scaling scheme for spallation.
%
% Syntax: scalingfactor = LSDspal(h,Rc,SPhi,w,nuclide);
%
% Where:
%   h = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   SPhi = solar modulation potntial (Phi, see source paper)
%   w = fractional water content of ground (nondimensional)
%
% Vectorized. Send in scalars or vectors of common length. 
%

% Modified by Jakob Heyman (jakob.heyman@gu.se) 2015-2018
% from code written by Nat Lifton 2013, Purdue University
% nlifton@purdue.edu
% Based on code by Greg Balco -- Berkeley Geochronology Lab
% balcs@bgc.org
% April, 2007
%
% Copyright 2001-2013, University of Washington, Purdue University
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).

EthRef = consts.ethfluxRef;
ThRef = consts.thfluxRef;

% Select reference values for nuclide of interest or flux
if nucl10 == 1
    BeRef = consts.P10nRef + consts.P10pRef;
end

if nucl26 == 1
    AlRef = consts.P26nRef + consts.P26pRef;
end

%if nuclide == 3
%    HeRef = consts.P3nRef + consts.P3pRef;
%elseif nuclide == 10
%    BeRef = consts.P10nRef + consts.P10pRef;
%elseif nuclide == 14
%    CRef = consts.P14nRef + consts.P14pRef;
%elseif nuclide == 26
%    AlRef = consts.P26nRef + consts.P26pRef;
%else
%    SpRef = consts.nfluxRef + consts.pfluxRef;
%    % Sato et al. (2008) Reference hadron flux integral >1 MeV
%end

% Site nucleon fluxes
NSite = Neutrons(h,Rc,SPhi,w,consts,nucl10,nucl26);
PSite = Protons(h,Rc,SPhi,consts,nucl10,nucl26);

%Nuclide-specific scaling factors as f(Rc)
if nucl10 == 1
    Site.sp10 = (NSite.P10n + PSite.P10p)./BeRef;
end

if nucl26 == 1
    Site.sp26 = (NSite.P26n + PSite.P26p)./AlRef;
end

%if nuclide == 3
%    Site.He = (NSite.P3n + PSite.P3p)./HeRef;
%elseif nuclide == 10
%    Site.Be = (NSite.P10n + PSite.P10p)./BeRef;
%elseif nuclide == 14
%    Site.C = (NSite.P14n + PSite.P14p)./CRef;
%elseif nuclide == 26
%    Site.Al = (NSite.P26n + PSite.P26p)./AlRef;
%else %Total nucleon flux scaling factors as f(Rc)
%    Site.sp = ((NSite.nflux + PSite.pflux))./SpRef; % Sato et al. (2008) Reference hadron flux
% integral >1 MeV
%end

%Site.E = NSite.E;%Nucleon flux energy bins
%Site.eth = ethflux./EthRef; %Epithermal neutron flux scaling factor as f(Rc)
%Site.th = thflux./ThRef;%Thermal neutron flux scaling factor as f(Rc)

LSDspal_out = Site;
