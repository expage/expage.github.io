function out = mucalib_Pmu(z,h,Rc,SPhi,nucl10,nucl26,consts)

% Calculates the production rate scaling of Be-10 and Al-26 by muons as a function of depth below
% the surface z (g/cm2) and site atmospheric pressure h (hPa), cutoff rigidity Rc (GV), and solar
% modulation parameter SPhi (MV).
%
% This code does not include fstar and sigma0 and it is made for calibration of these parameters
% 
% syntax out = mucalib_Pmu(z,h,Rc,SPhi,nucl10,nucl26,consts)
% 
% out.P_fast10   without sigma0!!
% out.P_neg10    without fstar!!
% out.P_fast26   without sigma0!!
% out.P_neg26    without fstar!!
%
% This uses the scheme in Heisinger and others (2002, 2 papers), scaled using the omnidirectional
% muon fluxes generated by PARMA, as presented in Sato et al. (2008). It then converts those scaled
% fluxes at the surface to vertical fluxes following Heisinger et al. (2002). The rest of the
% calculation generally parallels the procedure presented in the hard-copy documentation for the
% function "P_mu_total.m" as part of Balco et al., 2008.
%
% Production rate parameterization with fstar and sigma0 is based on the CRONUScalc calculator
% (Marrero et al. 2016; Phillips et al. 2016) calibration based on Antarctica depth profile.
% 
% Modified by Jakob Heyman (jakob.heyman@gu.se) 2015-2018
%
% Modified by Nat Lifton -- Purdue University 
% nlifton@purdue.edu
% March 2011
% from code written by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% March, 2006
%
% Copyright 2001-2013, University of Washington and Purdue University
%
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).

% what version is this?
ver = '201902';

% remember what direction the z vector came in
in_size = size(z);

% standardize direction
if size(z,1) > 1; z = z';end;

% figure the atmospheric depth in g/cm2
H = (1013.25 - h).*1.019716;
Href = 1013.25;

% find the omnidirectional flux at the site
mflux = Muons(h,Rc,SPhi);%Generates omnidirectional muon flux at site from Sato et al. (2008) model
mfluxRef = consts.mfluxRef;

phi_site = (mflux.neg(1,:) + mflux.pos(1,:));
phiRef = (mfluxRef.neg + mfluxRef.pos);

% changed from P_mu_totalLSD.m
% linear extrapolation of mflux.E, phi_site, and phiRef (see E2R below)
phi_site_d = (phi_site(1)-phi_site(2))/(mflux.E(1)-mflux.E(2));
phiRef_d = (phiRef(1)-phiRef(2))/(mflux.E(1)-mflux.E(2));
mflux.E = [5.710 mflux.E];
phi_site = [(mflux.E(1)-mflux.E(2))*phi_site_d+phi_site(1) phi_site];
phiRef = [(mflux.E(1)-mflux.E(2))*phiRef_d+phiRef(1) phiRef];

% find the vertical flux at SLHL
a = 258.5*(100.^2.66);
b = 75*(100.^1.66);
phi_vert_slhl = (a./((z+21000).*(((z+1000).^1.66) + b))).*exp(-5.5e-6 .* z);

% The above expression is only good to 2e5 g/cm2. We don't ever consider production
% below that depth. The full-depth scheme appears in the comments below.
% ------ begin full-depth flux equations -------
%phiz_1 = (a./((z+21000).*(((z+1000).^1.66) + b))).*exp(-5.5e-6 .* z);
%phiz_2 = 1.82e-6.*((121100./z).^2).*exp(-z./121100) + 2.84e-13;
%out(find(z<200000)) = phiz_1(find(z<200000));
%out(find(z>=200000)) = phiz_2(find(z>=200000));
% ------ end full-depth flux equations -------

% Convert E to Range
Temp = E2R(mflux.E); 
RTemp = Temp.R;

% find the stopping rate of vertical muons at site
SFmu = phi_site./phiRef;

% Prevent depths less than the minimum range in E2R to be used below
%z(z < min(RTemp)) = min(RTemp);
ztemp = z; % changed from P_mu_totalLSD.m to reduce near-surface production artifacts
ztemp(ztemp < min(RTemp)) = min(RTemp);

% Find scaling factors appropriate for energies associated with stopping muons at depths z
%Rz = interp1(RTemp,SFmu,z);
Rz = interp1(RTemp,SFmu,ztemp); % changed to reduce near-surface production artifacts

% Set upper limit to stopping range to test comparability with measurements
%StopLimit = 10;
% find the stopping rate of vertical muons at site
% find all ranges <10 g/cm2
%stopindex = find(RTemp<StopLimit,1,'last');
%SFmuslow = sum(phi_site(1:stopindex))./sum(phiRef(1:stopindex));
SFmuslow = phi_site(1)./phiRef(1); % changed to reduce near-surface production artifacts
Rz(Rz>SFmuslow) = SFmuslow;

RzSpline = spline(RTemp, SFmu);

% find the stopping rate of vertical muons at the site, scaled from SLHL
% this is done in a subfunction Rv0, because it gets integrated below.
R_vert_slhl = Rv0(z);
R_vert_site = R_vert_slhl.*Rz;

% find the flux of vertical muons at site
for a = 1:length(z);
	if z(a)/10 == round(z(a)/10);
        fprintf(1,'%.0f ',z(a));
    end;
    % integrate
    % ends at 200,001 g/cm2 to avoid being asked for an zero
    % range of integration -- 
    % get integration tolerance -- want relative tolerance around 1 part in 10^4. 
    tol = phi_vert_slhl(a) * 1e-4;
    [temp,fcnt] = quad(@(x) Rv0(x).*ppval(RzSpline,x),z(a),(2e5+1),tol);
    % second variable assignment here to preserve fcnt if needed
    phi_vert_site(a) = temp;
end;
   
% invariant flux at 2e5 g/cm2 depth - constant of integration
% calculated using commented-out formula above
phi_200k = (a./((2e5+21000).*(((2e5+1000).^1.66) + b))).*exp(-5.5e-6 .* 2e5);
phi_vert_site = phi_vert_site + phi_200k;

% find the total flux of muons at site

% angular distribution exponent
nofz = 3.21 - 0.297.*log((z+H)./100 + 42) + 1.21e-5.*(z+H);
% derivative of same
dndz = (-0.297./100)./((z+H)./100 + 42) + 1.21e-5;

phi_temp = phi_vert_site.*2.*pi./(nofz+1);

% that was in muons/cm2/s
% convert to muons/cm2/yr
phi = phi_temp*60*60*24*365;

% find the total stopping rate of muons at site
R_temp = (2.*pi./(nofz+1)).*R_vert_site ... 
    - phi_vert_site.*(-2.*pi.*((nofz+1).^-2)).*dndz;
    
% that was in total muons/g/s
% convert to negative muons/g/yr

R = R_temp*0.44*60*60*24*365;

% Attenuation lengths
LambdaMu =(Href-h)./(log(phi_site)-log(phiRef));

% Now calculate the production rates. 

% Depth-dependent parts of the fast muon reaction cross-section

% Per John Stone, personal communication 2011 - see text 
aalpha = 1.0;
Beta = 1.0;

Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*z)) + 50.7.*(1-exp(-5.05e-7.*z));

% internally defined constants
%if nuclide == 14
%	Natoms = consts.Natoms14;
%	k_negpartial = consts.k_negpartial14;
%	sigma0sp = consts.sigma0_14sp;
%	sigma0nu = consts.sigma0_14nu;
%	fstar_sp = consts.fstar14sp;
%	fstar_nu = consts.fstar14nu;
%elseif nuclide == 26
%	Natoms = consts.Natoms26;
%	k_negpartial = consts.k_negpartial26;
%	sigma0sp = consts.sigma0_26sp;
%	sigma0nu = consts.sigma0_26nu;
%	fstar_sp = consts.fstar26sp;
%	fstar_nu = consts.fstar26nu;
%else %defaults to 10Be
%	Natoms = consts.Natoms10;
%	k_negpartial = consts.k_negpartial10;
%	sigma0sp = consts.sigma0_10sp;
%	sigma0nu = consts.sigma0_10nu;
%	fstar_sp = consts.fstar10sp;
%	fstar_nu = consts.fstar10nu;
%end

if nucl10 == 1
	Natoms10 = consts.Natoms10;
	k_negpartial10 = consts.k_negpartial10;
%	sigma0nu10 = consts.sigma0_10nu;
%	fstar_nu10 = consts.fstar10nu;
	
	% fast muon production
	P_fast10 = phi.*Beta.*(Ebar.^aalpha).*Natoms10; % without sigma0nu10 !!
	
	% negative muon capture
	P_neg10 = R.*k_negpartial10; % without fstar_nu10 !!
end

if nucl26 == 1
	Natoms26 = consts.Natoms26;
	k_negpartial26 = consts.k_negpartial26;
%	sigma0nu26 = consts.sigma0_26nu;
%	fstar_nu26 = consts.fstar26nu;
	
	% fast muon production
	P_fast26 = phi.*Beta.*(Ebar.^aalpha).*Natoms26; % without sigma0nu26 !!
	
	% negative muon capture
	P_neg26 = R.*k_negpartial26; % without fstar_nu26 !!
end

% return
if nucl10 == 1
	out.P_fast10 = P_fast10 .* mflux.pint ./ mflux.pint(1);
	out.P_neg10 = P_neg10 .* mflux.nint ./ mflux.nint(1);
end
if nucl26 == 1
	out.P_fast26 = P_fast26 .* mflux.pint ./ mflux.pint(1);
	out.P_neg26 = P_neg26 .* mflux.nint ./ mflux.nint(1);
end

% -------------------------------------------------------------------------

function out = Rv0(z)

% this subfunction returns the stopping rate of vertically traveling muons
% as a function of depth z at sea level and high latitude.

a = exp(-5.5e-6.*z);
b = z + 21000;
c = (z + 1000).^1.66 + 1.567e5;
dadz = -5.5e-6 .* exp(-5.5e-6.*z);
dbdz = 1;
dcdz = 1.66.*(z + 1000).^0.66;

out = -5.401e7 .* (b.*c.*dadz - a.*(c.*dbdz + b.*dcdz))./(b.^2 .* c.^2);

% full depth calculation appears in comments below
%R_1 = -5.401e7 .* (b.*c.*dadz - a.*(c.*dbdz + b.*dcdz))./(b.^2 .* c.^2);
%f = (121100./z).^2;
%g = exp(-z./121100);
%dfdz = (-2.*(121100.^2))./(z.^3);
%dgdz = -exp(-z./121100)./121100;
%R_2 = -1.82e-6.*(g.*dfdz + f.*dgdz);
%out(find(z<200000)) = R_1(find(z<200000));
%out(find(z>=200000)) = R_2(find(z>=200000));

% -------------------------------------------------------------------------

function out = E2R(x)

% this subfunction returns the range and energy loss values for
% muons of energy E in MeV

% define range/energy/energy loss relation
% table for muons in standard rock
% http://pdg.lbl.gov/2010/AtomicNuclearProperties/ Table 281

data = [5.710 0.1 8.162
    1.0e1 8.400e-1 6.619
    1.4e1 1.530e0 5.180
    2.0e1 2.854e0 4.057
    3.0e1 5.687e0 3.157
    4.0e1 9.133e0 2.702
    8.0e1 2.675e1 2.029
    1.0e2 3.695e1 1.904
    1.4e2 5.878e1 1.779
    2.0e2 9.331e1 1.710
    3.0e2 1.523e2 1.688
    4.0e2 2.114e2 1.698
    8.0e2 4.418e2 1.775
    1.0e3 5.534e2 1.808
    1.4e3 7.712e2 1.862
    2.0e3 1.088e3 1.922
    3.0e3 1.599e3 1.990
    4.0e3 2.095e3 2.038
    8.0e3 3.998e3 2.152
    1.0e4 4.920e3 2.188
    1.4e4 6.724e3 2.244
    2.0e4 9.360e3 2.306
    3.0e4 1.362e4 2.383
    4.0e4 1.776e4 2.447
    8.0e4 3.343e4 2.654
    1.0e5 4.084e4 2.747
    1.4e5 5.495e4 2.925
    2.0e5 7.459e4 3.187
    3.0e5 1.040e5 3.611
    4.0e5 1.302e5 4.037
    8.0e5 2.129e5 5.748];

% units are range in g cm-2 (column 2)
% energy in MeV (column 1)
% Total energy loss/g/cm2 in MeV cm2/g(column 3)

% deal with zero situation

%too_low = find(x < 10);
too_low = find(x < 5.710); % changed to match extrapolation of data
x(too_low) = ones(size(too_low));

% obtain ranges
% use log-linear interpolation

% out = exp(spline(log(data(:,1)),log(data(:,2)),log(x')))';
out.R = exp(interpolate(log(data(:,1)),log(data(:,2)),log(x')))';
out.Eloss = exp(interpolate(log(data(:,1)),log(data(:,3)),log(x')))';
