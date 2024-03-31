function out = make_consts_expage();

% This function creates and saves a structure with relevant constants and external data for the
% expage calculator.
% Based on make_al_be_consts_v22.m (Balco et al. 2008) and make_consts_LSD.m Lifton et al. (2014),
% created and modified by Greg Balco, Brent Goehring, Nat Lifton.
%
% Paleomagnetic records updated to Lifton (2016) data: PMag_Sep12.mat changed to PavonPmag.mat.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2018-2024 (jakob.heyman@gu.se)

consts.version = '202403';
consts.prepdate = fix(clock);

% Be-10 decay constant -- new Euro value
consts.l10 = -log(0.5)./1.387e6; % Chmeleff/Korschinek value
dldt = -log(0.5).*(1.387e6^-2);
consts.l10unc = sqrt((dldt.*0.012e6)^2); % Chmeleff/Korschinek value

% Al-26 decay constant -- value compatible with Nishiizumi standards
% lambda = 9.83e-7 --> t(1/2) = 0.705e6 yr
% See Nishiizumi (2004) for details.
consts.l26 = 9.83e-7; 
consts.l26unc = 2.5e-8;

% C-14 decay constant -- t(1/2) = 5700 ± 30 yr (Kutschera 2013)
consts.l14 = -log(0.5)./5700;
dldt = -log(0.5).*(5700^-2);
consts.l14unc = sqrt((dldt.*30)^2);

% Effective attenuation length for spallation in rock
% Commonly accepted value: 160 g/cm2 (Gosse and Phillips 2001)
% Only used for simple age calculation - else rawattenuationlength is used.
consts.Lsp = 160;

% Be-10 standardization info.
% Standards comparison/conversion lookup table. A zero placeholder is allowed.
consts.std10 = {'07KNSTD','KNSTD','NIST_Certified','NIST_30000','NIST_30200','NIST_30300',...
    'NIST_30500','NIST_30600','NIST_27900','LLNL31000','LLNL10000','LLNL3000','LLNL1000',...
    'LLNL300','S555','S2007','BEST433','S555N','S2007N','BEST433N','STD11','0','-'};
consts.std10_cf = [1 0.9042 1.0425 0.9313 0.9251 0.9221 0.9157 0.9130 1 0.8761 0.9042 0.8644 ...
    0.9313 0.8562 0.9124 0.9124 0.9124 1 1 1 1 1 1]';

% Same for Al-26. A zero placeholder is allowed.
consts.std26 = {'KNSTD','ZAL94','ZAL94N','SMAL11','ASTER','Z92-0222','0','-'};
consts.std26_cf = [1 0.9134 1 1.021 1.021 1 1 1]';

% Reference production rates at SLHL for spallation. Letter codes sp and nu refer to the LSD
% spallation scaling and the LSD nuclide-specific scaling of Lifton et al. (2014), respectively.
% 10Be production rates are referenced to 07KNSTD.
% 26Al production rates are referenced to KNSTD.

% Be-10 production rates - global expage-202403 ref prod rate
consts.Pref10 = 4.06;
consts.Pref10unc = 0.25;
consts.Pref10iso = 4.08;
consts.Pref10isounc = 0.26;
consts.P10_ref_sp = 4; % not used and not calibrated!
consts.delP10_ref_sp = 0.3; % not used and not calibrated!

% Al-26 production rates - global expage-202403 ref prod rate
consts.Pref26 = 28.31;
consts.Pref26unc = 2.22;
consts.Pref26iso = 28.19;
consts.Pref26isounc = 2.35;
consts.P26_ref_sp = 28; % not used and not calibrated!
consts.delP26_ref_sp = 3; % not used and not calibrated!

% C-14 production rates - global expage-202403 ref prod rate
consts.Pref14 = 12.81;
consts.Pref14unc = 1.25;
consts.Pref14iso = 12.76;
consts.Pref14isounc = 1.16;

% Muon interaction cross-sections. All follow Heisinger (2002a,b).
consts.Natoms10 = 2.006e22;
consts.Natoms26 = 1.003e22;
consts.Natoms14 = 2.006e22;
consts.Natoms3 = 2.006e22; % not used in present version

% sigma0 and fstar calibrated for the expage calculator using mucalib.m (mucalib14.m for 14C)
% (ref prodrate expage-202403)
consts.sigma0_10nu = 1.792e-31;
consts.delsigma0_10nu = 0.399e-31;
consts.fstar10nu = 1.587e-3;
consts.delfstar10nu = 0.532e-3;
consts.sigma0_26nu = 2.651e-30;
consts.delsigma0_26nu = 1.218e-30;
consts.fstar26nu = 1.030e-2;
consts.delfstar26nu = 0.635e-2;
consts.sigma0_14nu = 0.792e-31;
consts.delsigma0_14nu = 125.111e-31;
consts.fstar14nu = 1.355e-1;
consts.delfstar14nu = 1.040e-1;

% sigma0 and fstar from CRONUScalc calculator (Marrero et al. 2016; Phillips et al. 2016)
% not used in present version
consts.sigma0_10sp = 0.252E-30;
consts.sigma0_26sp = 4.10E-30;
consts.delsigma0_10sp = 0.015e-30;
consts.delsigma0_26sp = 0.73e-30;
consts.sigma0_14sp = 8.79321E-30;
consts.fstar10sp = 1.89e-3;
consts.fstar26sp = 11.7e-3;
consts.delfstar10sp = 0.18e-3;
consts.delfstar26sp = 2.6e-3;
consts.fstar14sp = 0.137;

% k_negpartial from CRONUScalc calculator (Marrero et al. 2016; Phillips et al. 2016)
consts.k_negpartial10 = (0.704 * 0.1828)./1.106;
consts.k_negpartial26 = 0.296 * 0.6559;
consts.k_negpartial14 = 0.704 * 0.1828;
% delk_negpartial from CRONUScalc calculator (Marrero et al. 2016; Phillips et al. 2016)
consts.delk_negpartial10 = (0.704 * 0.1828 * 0.0003)./1.106;
consts.delk_negpartial26 = 0.296 * 0.6559 * 0.002;
consts.delk_negpartial14 = 0.704 * 0.1828 * 0.0011;

% New Spallogenic Nuclide Production Cross-Sections (n & p) from Bob Reedy 9/2010
load XSectsReedyAll;

consts.O16nxBe10 = O16nxBe10;
consts.O16pxBe10 = O16pxBe10;
consts.SinxBe10 = SinxBe10;
consts.SipxBe10 = SipxBe10;
consts.O16nn2pC14 = O16nn2pC14;
consts.O16pxC14 = O16pxC14;
consts.SinxC14 = SinxC14;
consts.SipxC14 = SipxC14;
consts.Aln2nAl26 = Aln2nAl26;
consts.AlppnAl26 = AlppnAl26;
consts.SinxAl26 = SinxAl26;
consts.SipxAl26 = SipxAl26;
consts.OnxHe3T = OnxHe3T;
consts.OpxHe3T = OpxHe3T;
consts.SinxHe3T = SinxHe3T;
consts.SipxHe3T = SipxHe3T;

% Paleomagnetic records for use in time-dependent production rate schemes
% Load the magnetic field data - Lifton (2016) version:
% 2010-1950: Rc from DGRFs
% 0-14 ka:   SHA.DIF.14k
% 14-75 ka:  GLOPIS-75
% >75 ka:    PADM2M
load PavonPmag;
% Includes updated SPhi values from Usoskin et al, 2011, Journal of Geophysical Research, v. 116,
% no. A2. End value of tM = 1e7 to enable to run in MATLAB 2012a and later. 

% Relative dipole moment and time vector
consts.M = GLOPADMM0; 
consts.tM = tM; 

% These start at 14000 yr -- time slices are 100-yr from 0 to 3000, 200-yr from 3000 to 75000,
% 1000-yr from 75-ka to 800-ka, 2000-yr from 800-ka to 2-Ma, and 10-Ma at the end. The interval
% spacing is adjusted from that of Lifton et al. (2014) to improve computational speed while
% remaining faithful to the temporal variability in each model.

% Cutoff rigidity blocks for past 14000 yr. 
% PavonRc has lat x lon x time blocks of Rc values for the past 14000 years.
% Cutoff rigidity obtained by Nat Lifton using trajectory tracing from the SHA.DIF.14k magnetic
% field model of Pavón-Carrasco et al. (2014).  
consts.PavonRc = PavonRc; % data block
consts.lat_Rc = lat_Rc; % lat and lon indices for Rc data block
consts.lon_Rc = lon_Rc;
consts.tRc = tRc; % time vector for Rc data block

% Solar variability from Usoskin et al. 2011
% 0-11260 yr

% Per Tatsuhiko Sato, personal communication, 2013, convert annually averaged Usoskin et al. (2011)
% solar modulation potential to Sato Force Field Potential due to different assumed Local
% Interstellar Spectrum and other factors
SPhi = 1.1381076.*SPhi - 1.2738468e-4.*SPhi.^2;

consts.SPhi = SPhi;
consts.SPhiInf = trapz(tRc(1:78),SPhi)/tRc(78); % time-averaged SPhi

load Reference;
% Reference values for scaling via Sato et al. (2008) spectra
consts.E = E;
consts.P3nRef = P3nRef; %3He neutron reference production in SiO2
consts.P3pRef = P3pRef; %3He proton reference production in SiO2
consts.P10nRef = P10nRef; %10Be neutron reference production in SiO2
consts.P10pRef = P10pRef; %10Be proton reference production in SiO2
consts.P14nRef = P14nRef; %14C neutron reference production in SiO2
consts.P14pRef = P14pRef; %14C proton reference production in SiO2
consts.P26nRef = P26nRef; %26Al neutron reference production in SiO2
consts.P26pRef = P26pRef; %26Al proton reference production in SiO2
consts.nfluxRef = nfluxRef; %integral reference neutron flux
consts.pfluxRef = pfluxRef; %integral reference proton flux
consts.ethfluxRef = ethfluxRef; %integral reference epithermal flux
consts.thfluxRef = thfluxRef; %integral reference thermal neutron flux
consts.mfluxRef = mfluxRef; % reference muon flux components

% Finish up
save consts_expage consts
disp(['consts_expage v. ' consts.version ' saved']);
