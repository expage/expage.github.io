function out = P_sp_expage(h,Rc,s,w,consts,nucl10,nucl26,nucl14)

% Calculates spallation scaling factors as a function of atmospheric pressure h (hPa), cutoff
% rigidity Rc (GV), solar modulation parameter s (MV), and water content w (0-1).
% Function based on Neutrons.m, Protons.m, and LSDscaling.m from Lifton et al. (2014).
% Sato et al. (2008) Analytical Function Approximation (PARMA)
% Implemented in MATLAB by Nat Lifton, 2013
% Purdue University, nlifton@purdue.edu
% modified by Jakob Heyman (jakob.heyman@gu.se) 2016-2023

% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).

% what version is this?
ver = '202306';

% make h,Rc,s row vectors
h = h(:)';
Rc = Rc(:)';
s = s(:)';

% common Neutron and Proton parameters ===============================
% extend h if shorter than Rc
if length(h)<length(Rc);
    h(end+1:length(Rc)) = h(1);
end;
x = h.*1.019716; % Convert pressure (hPa) to atm depth (g/cm2)

% E = logspace(-8,5,1000);
% E = [1.1295 11.295 112.95 1129.5 11295];
E = logspace(0,5.3010,200);

% Flatten low rigidities.
lowRc = find(Rc < 1.0);
Rc(lowRc) = 1.0 + zeros(size(lowRc));

smin = 400; %units of MV
smax = 1200; %units of MV

NP.P10n = zeros(1,length(Rc));
NP.P10p = zeros(1,length(Rc));

% Neutron parameters =================================================
% Integrated neutron flux <15 MeV
NEt = 2.5e-8; % Thermal Neutron Energy in MeV

Na6 = 1.8882e-4;
Na7 = 4.4791e-1;
Na8 = 1.4361e-3;
Na12 = 1.4109e-2;

Nb11min = 2.5702e1;
Nb11max = -6.9221;
Nb12min = -5.0931e-1;
Nb12max = 1.1336;
Nb13min= 7.4650;
Nb13max = 2.6961e1;
Nb14min = 1.2313e1;
Nb14max = 1.1746e1;
Nb15min = 1.0498;
Nb15max = 2.6171;
Nb21min = 6.5143e-3;
Nb21max = 5.3211e-3;
Nb22min = 3.3511e-5;
Nb22max = 8.4899e-5;
Nb23min = 9.4415e-4;
Nb23max = 2.0704e-3;
Nb24min = 1.2088e1;
Nb24max = 1.1714e1;
Nb25min = 2.7782;
Nb25max = 3.8051;
Nb31min = 9.8364e-1;
Nb31max = 9.7536e-1;
Nb32min = 1.4964e-4;
Nb32max = 6.4733e-4;
Nb33min = -7.3249e-1;
Nb33max = -2.2750e-1;
Nb34min = -1.4381;
Nb34max = 2.0713;
Nb35min = 2.7448;
Nb35max = 2.1689;
Nb41min = 8.8681e-3;
Nb41max = 9.1435e-3;
Nb42min = -4.3322e-5;
Nb42max = -6.4855e-5;
Nb43min = 1.7293e-2;
Nb43max = 5.8179e-3;
Nb44min = -1.0836;
Nb44max = 1.0168;
Nb45min = 2.6602;
Nb45max = 2.4504;

Nb121 = 9.31e-1;
Nb122 = 3.70e-2;
Nb123 = -2.02;
Nb124 = 2.12;
Nb125 = 5.34;
Nb131 = 6.67e-4;
Nb132 = -1.19e-5;
Nb133 = 1.00e-4;
Nb134 = 1.45;
Nb135 = 4.29;

% Basic Spectrum
Nb51 = 9.7337e-4;
Nb52 = -9.6635e-5;
Nb53 = 1.2121e-2;
Nb54 = 7.1726;
Nb55 = 1.4601;
Nb91 = 5.7199e2;
Nb92 = 7.1293;
Nb93 = -1.0703e2;
Nb94 = 1.8538;
Nb95 = 1.2142;
Nb101 = 6.8552e-4;
Nb102 = 2.7136e-5;
Nb103 = 5.7823e-4;
Nb104 = 8.8534;
Nb105 = 3.6417;
Nb111 = -5.0800e-1;
Nb112 = 1.4783e-1;
Nb113 = 1.0068;
Nb114 = 9.1556;
Nb115 = 1.6369;

Nc1 = 2.3555e-1; % lethargy^-1
Nc2 = 2.3779; % MeV
Nc3 = 7.2597e-1;
Nc5 = 1.2391e2; % MeV
Nc6 = 2.2318; % MeV
Nc7 = 1.0791e-3; % lethargy^-1
Nc8 = 3.6435e-12; % MeV
Nc9 = 1.6595;
Nc10 = 8.4782e-8; % MeV
Nc11 = 1.5054;

% Ground-Level Spectrum
Nh31 = -2.5184e1;
Nh32 = 2.7298;
Nh33 = 7.1526e-2;
Nh51 = 3.4790e-1;
Nh52 = 3.3493;
Nh53 = -1.5744;

Ng1 = -0.023499;
Ng2 = -0.012938;
Ng3 = 10.^(Nh31 + Nh32./(w + Nh33));
Ng4 = 9.6889e-1;
Ng5 = Nh51 + Nh52.*w + Nh53.*(w.^2);

NfG = 10.^(Ng1 + Ng2.*log10(E./Ng3).*(1-tanh(Ng4.*log10(E./Ng5))));

% Thermal Neutron Spectrum
Nh61 = 1.1800e-1;
Nh62 = 1.4438e-1;
Nh63 = 3.8733;
Nh64 = 6.5298e-1;
Nh65 = 4.2752e1;

Ng6 = (Nh61 + Nh62.*exp(-Nh63.*w))./(1 + Nh64.*exp(-Nh65.*w));

NPhiT = Ng6.*((E./NEt).^2).*exp(-E./NEt);
% end Neutron parameters =============================================


% Proton parameters ==================================================
PA = 1;
PZ = 1;
PEp = 938.27; % Rest mass of a proton
PU = (4-1.675).*pi.*PA./PZ.*1e-7; % Unit conversion factor

% Primary spectrum
Pa1 = 2.1153;
Pa2 = 4.4511e-1;
Pa3 = 1.0064e-2;
Pa4 = 3.9564e-2;
Pa5 = 2.9236;
Pa6 = 2.7076;
Pa7 = 1.2663e4;
Pa8 = 4.8288e3;
Pa9 = 3.2822e4;
Pa10 = 7.4378e3;
Pa11 = 3.4643;
Pa12 = 1.6752;
Pa13 = 1.3691;
Pa14 = 2.0665;
Pa15 = 1.0833e2;
Pa16 = 2.3013e3;

% Secondary Spectrum
Pc11 = 1.2560;
Pc12 = 3.2260e-3;
Pc13 = -4.8077e-6;
Pc14 = 2.2825e-9;
Pc21 = 4.3783e-1;
Pc22 = -5.5759e-4;
Pc23 = 7.8388e-7;
Pc24 = -3.8671e-10;
Pc31 = 1.8102e-4;
Pc32 = -5.1754e-7;
Pc33 = 7.5876e-10;
Pc34 = -3.8220e-13;
Pc41 = 1.7065;
Pc42 = 7.1608e-4;
Pc43 = -9.3220e-7;
Pc44 = 5.2665e-10;

Pb1 = Pc11 + Pc12.*x + Pc13.*x.^2 + Pc14.*x.^3;
Pb2 = Pc21 + Pc22.*x + Pc23.*x.^2 + Pc24.*x.^3;
Pb3 = Pc31 + Pc32.*x + Pc33.*x.^2 + Pc34.*x.^3;
Pb4 = Pc41 + Pc42.*x + Pc43.*x.^2 + Pc44.*x.^3;

Ph11min = 2.4354e-3;
Ph11max = 2.5450e-3;
Ph12min = -6.0339e-5;
Ph12max = -7.1807e-5;
Ph13min= 2.1951e-3;
Ph13max = 1.4580e-3;
Ph14min = 6.6767;
Ph14max = 6.9150;
Ph15min = 9.3228e-1;
Ph15max = 9.9366e-1;
Ph21min = 7.7872e-3;
Ph21max = 7.6828e-3;
Ph22min = -9.5771e-6;
Ph22max = -2.4119e-6;
Ph23min = 6.2229e-4;
Ph23max = 6.6411e-4;
Ph24min = 7.7842;
Ph24max = 7.7461;
Ph25min = 1.8502;
Ph25max = 1.9431;
Ph31min = 9.6340e-1;
Ph31max = 9.7353e-1;
Ph32min = 1.5974e-3;
Ph32max = 1.0577e-3;
Ph33min = -7.1179e-2;
Ph33max = -2.1383e-2;
Ph34min = 2.2320;
Ph34max = 3.0058;
Ph35min = 7.8800e-1;
Ph35max = 9.1845e-1;
Ph41min = 7.8132e-3;
Ph41max = 7.3482e-3;
Ph42min = 9.7085e-11;
Ph42max = 2.5598e-5;
Ph43min = 8.2392e-4;
Ph43max = 1.2457e-3;
Ph44min = 8.5138;
Ph44max = 8.1896;
Ph45min = 2.3125;
Ph45max = 2.9368;

Ph51 = 1.9100e-1;
Ph52 = 7.0300e-2;
Ph53 = -6.4500e-1;
Ph54 = 2.0300;
Ph55 = 1.3000;
Ph61 = 5.7100e-4;
Ph62 = 6.1300e-6;
Ph63 = 5.4700e-4;
Ph64 = 1.1100;
Ph65 = 8.3700e-1;
% end Proton parameters ==============================================

% Make sure the clip index is consistent with the definition of E above
clipindex = find(E <= 1, 1, 'last' );

%for a = 1:length(Rc)
% initial Neutron calculations =====================================================================
Na1min = Nb11min + Nb12min.*Rc + Nb13min./(1 + exp((Rc - Nb14min)./Nb15min));
Na1max = Nb11max + Nb12max.*Rc + Nb13max./(1 + exp((Rc - Nb14max)./Nb15max));
Na2min = Nb21min + Nb22min.*Rc + Nb23min./(1 + exp((Rc - Nb24min)./Nb25min));
Na2max = Nb21max + Nb22max.*Rc + Nb23max./(1 + exp((Rc - Nb24max)./Nb25max));
Na3min = Nb31min + Nb32min.*Rc + Nb33min./(1 + exp((Rc - Nb34min)./Nb35min));
Na3max = Nb31max + Nb32max.*Rc + Nb33max./(1 + exp((Rc - Nb34max)./Nb35max));
Na4min = Nb41min + Nb42min.*Rc + Nb43min./(1 + exp((Rc - Nb44min)./Nb45min));
Na4max = Nb41max + Nb42max.*Rc + Nb43max./(1 + exp((Rc - Nb44max)./Nb45max));

Na5 = Nb51 + Nb52.*Rc + Nb53./(1 + exp((Rc - Nb54)./Nb55));
Na9 = Nb91 + Nb92.*Rc + Nb93./(1 + exp((Rc - Nb94)./Nb95));
Na10 = Nb101 + Nb102.*Rc + Nb103./(1 + exp((Rc - Nb104)./Nb105));
Na11 = Nb111 + Nb112.*Rc + Nb113./(1 + exp((Rc - Nb114)./Nb115));

Nb5 = Nb121 + Nb122.*Rc + Nb123./(1 + exp((Rc - Nb124)./Nb125));
Nb6 = Nb131 + Nb132.*Rc + Nb133./(1 + exp((Rc - Nb134)./Nb135));

Nc4 = Na5 + Na6.*x./(1 + Na7.*exp(Na8.*x)); % lethargy^-1
Nc12 = Na9.*(exp(-Na10.*x) + Na11.*exp(-Na12.*x)); % MeV

NPhiLmin = Na1min.*(exp(-Na2min.*x) - Na3min.*exp(-Na4min.*x)); %length of Rc
NPhiLmax = Na1max.*(exp(-Na2max.*x) - Na3max.*exp(-Na4max.*x)); %length of Rc

Nf3 = Nb5 + Nb6.*x;
Nf2 = (NPhiLmin - NPhiLmax)./(smin.^Nf3 - smax.^Nf3);
Nf1 = NPhiLmin - Nf2.*smin.^Nf3;

NPhiL = Nf1 + Nf2.*s.^Nf3;
% end initial Neutron calculations =================================================================

% initial Proton calculations ======================================================================
Pg1min = Ph11min + Ph12min.*Rc + Ph13min./(1 + exp((Rc - Ph14min)./Ph15min));
Pg1max = Ph11max + Ph12max.*Rc + Ph13max./(1 + exp((Rc - Ph14max)./Ph15max));
Pg2min = Ph21min + Ph22min.*Rc + Ph23min./(1 + exp((Rc - Ph24min)./Ph25min));
Pg2max = Ph21max + Ph22max.*Rc + Ph23max./(1 + exp((Rc - Ph24max)./Ph25max));
Pg3min = Ph31min + Ph32min.*Rc + Ph33min./(1 + exp((Rc - Ph34min)./Ph35min));
Pg3max = Ph31max + Ph32max.*Rc + Ph33max./(1 + exp((Rc - Ph34max)./Ph35max));
Pg4min = Ph41min + Ph42min.*Rc + Ph43min./(1 + exp((Rc - Ph44min)./Ph45min));
Pg4max = Ph41max + Ph42max.*Rc + Ph43max./(1 + exp((Rc - Ph44max)./Ph45max));

PphiPmin = Pg1min.*(exp(-Pg2min.*x) - Pg3min.*exp(-Pg4min.*x)); %length of Rc
PphiPmax = Pg1max.*(exp(-Pg2max.*x) - Pg3max.*exp(-Pg4max.*x)); %length of Rc

Pg5 = Ph51 + Ph52.*Rc + Ph53./(1 + exp((Rc - Ph54)./Ph55));
Pg6 = Ph61 + Ph62.*Rc + Ph63./(1 + exp((Rc - Ph64)./Ph65));

Pf3 = Pg5 + Pg6.*x;
Pf2 = (PphiPmin - PphiPmax)./(smin.^Pf3 - smax.^Pf3);
Pf1 = PphiPmin - Pf2.*smin.^Pf3;

PphiP = Pf1 + Pf2.*s.^Pf3;

PEc = (sqrt((1000.*Rc.*PZ).^2 + PEp.^2) - PEp)./PA;
PEs = Pa13.*(PEc - Pa14.*x);
PEs1 = max(Pa15,PEs);
PEs2 = max(Pa16,PEs);
% end initial Proton calculations ==================================================================

% final calculations ===============================================================================
for a = 1:length(Rc)
    % Neutron calculations =========================================================================
    NPhiB = (Nc1.*(E./Nc2).^Nc3).*exp(-E./Nc2) + Nc4(a).*exp((-(log10(E) - log10(Nc5)).^2) ./ ...
        (2.*(log10(Nc6)).^2)) + Nc7.*log10(E./Nc8).*(1 + tanh(Nc9.*log10(E./Nc10))) .* ...
        (1 - tanh(Nc11.*log10(E./Nc12(a))));
    
    NPhiG = NPhiL(a).*(NPhiB.*NfG + NPhiT);
    NPhiGMev = NPhiG./E;
    % end Neutron calculations =====================================================================
    
    % Proton calculations ==========================================================================
    PEtoa = E + Pa1.*x(a);
    PRtoa = 0.001.*sqrt((PA.*PEtoa).^2 + 2.*PA.*PEp.*PEtoa)./PZ;

    PElis = PEtoa + s(a).*PZ./PA;
    PBeta = sqrt(1-(PEp./(PEp + PElis.*PA)).^2); % Particle speed relative to light
    PRlis = 0.001.*sqrt((PA.*PElis).^2 + 2.*PA.*PEp.*PElis)./PZ;
    PC = Pa7 + Pa8./(1 + exp((PElis - Pa9)./Pa10));

    PphiTOA = (PC.*(PBeta.^Pa5)./(PRlis.^Pa6)).*(PRtoa./PRlis).^2;
    PphiPri = (PU./PBeta).*PphiTOA.*(Pa2.*exp(-Pa3.*x(a)) + (1 - Pa2).*exp(-Pa4.*x(a)));
    
    PphiSec = (PphiP(a).*Pb1(a).*E.^Pb2(a))./(1 + Pb3(a).*E.^Pb4(a));
    
    PphiPtot = PphiPri.*(tanh(Pa11.*(E./PEs1(a) - 1)) + 1)./2 + ...
        PphiSec.*(tanh(Pa12.*(1 - E./PEs2(a))) + 1)./2;
    % end Proton calculations ======================================================================
    
    if nucl10 == 1
		NP.P10n(a) = (trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
            consts.O16nxBe10(clipindex:end)) + trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
            consts.SinxBe10(clipindex:end)./2)).*consts.Natoms10.*1e-27.*3.1536e7;
        NP.P10p(a) = (trapz(E,PphiPtot.*consts.O16pxBe10) + ...
		    trapz(E,PphiPtot.*consts.SipxBe10./2)).*consts.Natoms10.*1e-27.*3.1536e7;
	end
	
	if nucl26 == 1
		NP.P26n(a) = trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
            consts.SinxAl26(clipindex:end)).*consts.Natoms26.*1e-27.*3.1536e7;
        NP.P26p(a) = trapz(E,PphiPtot.*consts.SipxAl26).*consts.Natoms26.*1e-27.*3.1536e7; 
	end

    if nucl14 == 1
        NP.P14n(a) = (trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
            consts.O16nn2pC14(clipindex:end))+ trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
            consts.SinxC14(clipindex:end)./2)).*consts.Natoms14.*1e-27.*3.1536e7;
        NP.P14p(a) = (trapz(E,PphiPtot.*consts.O16pxC14)+ ...
            trapz(E,PphiPtot.*consts.SipxC14./2)).*consts.Natoms14.*1e-27.*3.1536e7;
    end;
	
%    if nuclide == 3
%        NP.P3n(a) = (trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
%            consts.OnxHe3T(clipindex:end)) + trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
%            consts.SinxHe3T(clipindex:end)./2)).*consts.Natoms3.*1e-27.*3.1536e7;
%        NP.P3p(a) = (trapz(E,PphiPtot.*consts.OpxHe3T) + ...
%            trapz(E,PphiPtot.*consts.SipxHe3T./2)).*consts.Natoms3.*1e-27.*3.1536e7;    
%    elseif nuclide == 10
%        NP.P10n(a) = (trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
%            consts.O16nxBe10(clipindex:end)) + trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
%            consts.SinxBe10(clipindex:end)./2)).*consts.Natoms10.*1e-27.*3.1536e7; 
%        NP.P10p(a) = (trapz(E,PphiPtot.*consts.O16pxBe10) + ...
%            trapz(E,PphiPtot.*consts.SipxBe10./2)).*consts.Natoms10.*1e-27.*3.1536e7;
%    elseif nuclide == 14
%        NP.P14n(a) = (trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
%            consts.O16nn2pC14(clipindex:end))+ trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
%            consts.SinxC14(clipindex:end)./2)).*consts.Natoms14.*1e-27.*3.1536e7;
%        NP.P14p(a) = (trapz(E,PphiPtot.*consts.O16pxC14)+ ...
%            trapz(E,PphiPtot.*consts.SipxC14./2)).*consts.Natoms14.*1e-27.*3.1536e7;
%    elseif nuclide == 26
%        NP.P26n(a) = trapz(E(clipindex:end),NPhiGMev(clipindex:end).*...
%            consts.SinxAl26(clipindex:end)).*consts.Natoms26.*1e-27.*3.1536e7;
%        NP.P26p(a) = trapz(E,PphiPtot.*consts.SipxAl26).*consts.Natoms26.*1e-27.*3.1536e7; 
%    else
%        NP.nflux(a) = trapz(E(clipindex:end),NPhiGMev(clipindex:end));
%        NP.pflux(a) = trapz(E(clipindex:end),phiPtot(clipindex:end));
%    end
end

% Select reference values for nuclide of interest or flux and calculate scaling factor
if nucl10 == 1
    BeRef = consts.P10nRef + consts.P10pRef;
    out.sp10 = (NP.P10n + NP.P10p)./BeRef;
end
if nucl26 == 1
    AlRef = consts.P26nRef + consts.P26pRef;
    out.sp26 = (NP.P26n + NP.P26p)./AlRef;
end
if nucl14 == 1
    CRef = consts.P14nRef + consts.P14pRef;
    out.sp14 = (NP.P14n + NP.P14p)./CRef;
end
