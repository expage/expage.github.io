function burial();

% expage 10Be and 26Al exposure/burial duration calculator.
% Read and fix input, use get_1026_burial to calculate burial and exposure durations, and fix and
% write output. If plotting is on, display P-normalized 26/10 banana plot.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

tic();

% plotting? (1 = yes) ==============================================================================
expline = 1; % plot simple exposure line and burial lines
eline = 1; % plot erosion end point line
points = 1; % plot sample points
sigma1line = 1; % plot sample 1 sigma line
sigma2line = 0; % plot sample 2 sigma line
probcontours = 0; % plot probability contour lines
% ==================================================================================================

mc = 1E5; % Monte Carlo size

% write normalized conc data? (for later plotting)
normNout = 1;

% What version is this?
ver = '201806';

% read input file
[samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
    samplein.rho,samplein.othercorr,samplein.E,samplein.N10,samplein.delN10,samplein.be_stds,...
    samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr] = ...
    textread('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n','commentstyle','matlab');

% run and load expage constants
make_consts_expage;
load consts_expage;

% convert 10Be concentrations according to standards
for i = 1:numel(samplein.N10);
    be_mult(i,1) = consts.be_stds_cfs(strcmp(samplein.be_stds(i),consts.be_stds_names));
end;
samplein.N10 = samplein.N10 .* be_mult;
samplein.delN10 = samplein.delN10 .* be_mult;

% convert 26Al concentrations according to standards
for i = 1:numel(samplein.N26);
    al_mult(i,1) = consts.al_stds_cfs(strcmp(samplein.al_stds(i),consts.al_stds_names));
end;
samplein.N26 = samplein.N26 .* al_mult;
samplein.delN26 = samplein.delN26 .* al_mult;

% fix longitude values
samplein.long(find(samplein.long < 0)) = samplein.long(find(samplein.long < 0)) + 360;

% fix sample pressure
std_v = strcmp(samplein.aa,'std');
ant_v = strcmp(samplein.aa,'ant');
pre_v = strcmp(samplein.aa,'pre');
samplein.pressure(std_v) = ERA40atm(samplein.lat(std_v),samplein.long(std_v),samplein.elv(std_v));
samplein.pressure(ant_v) = antatm(samplein.elv(ant_v));
samplein.pressure(pre_v) = samplein.elv(pre_v);

% fix erosion rate unit (mm/ka -> cm/yr)
samplein.E = samplein.E .* 1E-4;

% decay constants
l10 = consts.l10;
l26 = consts.l26;

% declare pl parameters
pl.N10n = []; pl.delN10n = []; pl.N26n = []; pl.delN26n = [];
pl.P10sp = []; pl.P26sp = []; pl.P_mu10 = []; pl.P_mu26 = [];
pl.Lsp = []; pl.RcEst = []; pl.SPhiAv = []; pl.thick_rho = []; pl.atm = []; pl.othercorr = [];

% fix for output
output(1,1) = {'sample'};
if sum(samplein.N10+samplein.N26) > 0;
    output(1,2:5) = {'exposure(yr)','exp-unc(yr)','burial(yr)','bur-unc(yr)'};
    if normNout > 0; % if writing normalized conc data
        output(1,6:9) = {'N10norm','N10unc','N26norm','N26unc'};
    end;
end;

% pick out samples one by one
for i = 1:numel(samplein.lat);
    sample.sample_name = samplein.sample_name(i,:);
    sample.lat = samplein.lat(i);
    sample.long = samplein.long(i);
    sample.elv = samplein.elv(i);
    sample.aa = samplein.aa(i,:);
    sample.thick = samplein.thick(i);
    sample.rho = samplein.rho(i);
    sample.othercorr = samplein.othercorr(i);
    sample.E = samplein.E(i);
    sample.N10 = samplein.N10(i);
    sample.delN10 = samplein.delN10(i);
    sample.be_stds = samplein.be_stds(i,:);
    sample.N26 = samplein.N26(i);
    sample.delN26 = samplein.delN26(i);
    sample.al_stds = samplein.al_stds(i,:);
    sample.samplingyr = samplein.samplingyr(i);
    sample.pressure = samplein.pressure(i);
    
    % write sample name to output
    output(i+1,1) = sample.sample_name;
    
    % check that there is 10Be and 26Al data
    if sample.N10 .* sample.N26 <= 0;
        continue;
    end;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample_name{1});
    
    % do calculations
    results = get_1026_burial(sample,consts,mc);
    
    % write to output
    output(i+1,2) = {num2str(results.expdur,'%.0f')};
    output(i+1,3) = {num2str(results.del_expdur,'%.0f')};
    output(i+1,4) = {num2str(results.burial,'%.0f')};
    output(i+1,5) = {num2str(results.del_burial,'%.0f')};
    
    if normNout > 0; % if writing normalized conc data
        output(i+1,6) = {num2str(results.N10norm,'%.6f')};
        output(i+1,7) = {num2str(results.N10norm_del,'%.6f')};
        output(i+1,8) = {num2str(results.N26norm,'%.6f')};
        output(i+1,9) = {num2str(results.N26norm_del,'%.6f')};
    end;
    
    % display exposure and fix if saturated
    if results.del_expdur >= 0; % if initial exposure within bounds
        fprintf(1,' \texposure = %s ± %s yr',output{i+1,2},output{i+1,3});
    else;
        output(i,2) = strcat('>',num2str(results.expdur,'%.0f'));
        output(i,3) = '-';
        fprintf(1,' \texposure: %s yr (saturated 10Be)',output{i+1,2});
    end;
    
    % display burial and fix if too much/too negative
    if results.del_burial >= 0; % if burial within bounds
        % display results
        fprintf(1,' \tburial = %s ± %s yr\n',output{i+1,4},output{i+1,5});
    else;
        if results.del_burial == -1; % if too much burial
            burchar = '>';
        elseif results.del_burial == -2; % if too much negative burial
            burchar = '<';
        end;
        output(i,4) = strcat(burchar,num2str(results.burial,'%.0f'));
        output(i,5) = '-';
        fprintf(1,' \tburial: %s yr\n',output{i+1,4});
    end;
    
    % for plotting
    pl.N10n(end+1,1) = results.N10norm;
    pl.delN10n(end+1,1) = results.N10norm_del;
    pl.N26n(end+1,1) = results.N26norm;
    pl.delN26n(end+1,1) = results.N26norm_del;
    pl.P10sp(end+1,1) = results.P10sp;
    pl.P26sp(end+1,1) = results.P26sp;
    pl.P_mu10(end+1,1) = results.P_mu10;
    pl.P_mu26(end+1,1) = results.P_mu26;
    pl.Lsp(end+1,1) = results.Lsp;
    pl.RcEst(end+1,1) = results.RcEst;
    pl.thick_rho(end+1,1) = sample.thick * sample.rho;
    pl.atm(end+1,1) = sample.pressure;
    pl.othercorr(end+1,1) = sample.othercorr;
    
    clear sample;
end;

% fix and save output =============================================================
if sum(samplein.N10 + samplein.N26) > 0;
    % fix output string
    outstr = '%s';
    for j = 1:size(output,2)-1;
        outstr = strcat(outstr,'\t%s');
    end;
    outstr = strcat(outstr,'\n');

    % fill empty cells with '-'
    nullidx = cellfun(@isempty,output);
    output(nullidx) = {'-'};
    
    % write out-expage.txt
    out = fopen('out-burial.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    fclose(out);
end;
% end output block ================================================================

% plotting ========================================================================
if (expline + eline + points + sigma1line + probcontours) > 0;
    plotting(expline,eline,points,sigma1line,sigma2line,probcontours,pl,consts);
end;
% end plotting ====================================================================

toc()
% end burial function ==============================================================================


% subfunction get_1026_burial ======================================================================
function results = get_1026_burial(sample,consts,mc);
% This function calculates the simple albe exposure and burial durations and packages the results.

% Select appropriate values for nuclide of interest
% Atoms/g measurement
N10 = sample.N10; delN10 = sample.delN10;
N26 = sample.N26; delN26 = sample.delN26;
% Ref production rates
Pref10 = consts.P10_ref_nu; delPref10 = consts.delP10_ref_nu;
Pref26 = consts.P26_ref_nu; delPref26 = consts.delP26_ref_nu;
% Decay constants
l10 = consts.l10;
l26 = consts.l26;

% fix Rc SPhi and tv
LSDfix = LSD_fix(sample.lat,sample.long,1E7,-1,consts);

% interpolate Lsp (Sato, 2008; Marrero et al., 2016)
Lsp = rawattenuationlength(sample.pressure,LSDfix.RcEst);

% 1a. Thickness scaling factor.
if sample.thick > 0;
    thickSF = (Lsp./(sample.rho.*sample.thick)).*(1 - exp(((-1.*sample.rho.*sample.thick)./Lsp)));
else 
    thickSF = 1;
end;

% nu production scaling
LSDnu = LSDspal(sample.pressure,LSDfix.RcEst,consts.SPhiInf,LSDfix.w,1,1,consts);

% randomized ref prod rates for Monte Carlo simulation
Pref10mc = normrnd(Pref10,delPref10,[mc,1]);
Pref26mc = normrnd(Pref26,delPref26,[mc,1]);

% spallation prod rates
P_nu10 = LSDnu.Be.*Pref10.*thickSF.*sample.othercorr;
P_nu26 = LSDnu.Al.*Pref26.*thickSF.*sample.othercorr;
P_nu10mc = LSDnu.Be.*Pref10mc.*thickSF.*sample.othercorr;
P_nu26mc = LSDnu.Al.*Pref26mc.*thickSF.*sample.othercorr;

% muon production
P_mu = P_mu_LSD((sample.thick.*sample.rho./2),sample.pressure,LSDfix.RcEst,consts.SPhiInf,1,1,...
    consts,'no');
P_mu10 = P_mu.Be.*sample.othercorr;
P_mu26 = P_mu.Al.*sample.othercorr;

% summed prod rates
P10 = P_nu10 + P_mu.Be;
P26 = P_nu26 + P_mu.Al;
P10mc = P_nu10mc + P_mu.Be;
P26mc = P_nu26mc + P_mu.Al;

% normalized N (sent out in result for plotting)
N10norm = N10./(P_nu10 + P_mu10);
N10norm_del = delN10./(P_nu10 + P_mu10);
N26norm = N26./(P_nu26 + P_mu26);
N26norm_del = delN26./(P_nu26 + P_mu26);

% randomized 10Be and 26Al concentrations for mc simulation
N10mc = normrnd(N10,delN10,[mc,1]);
N26mc = normrnd(N26,delN26,[mc,1]);

% test burial from -10M to 10M with 1M yr step
tv1 = (-10:1:10).*1e6;
burtest1 = 0;
burtest1mc = zeros(mc,1);
for i = 1:21;
    % calculate N before burial
    N10t1 = N10./exp(-l10.*tv1(i));
    N26t1 = N26./exp(-l26.*tv1(i));
    N10t1mc = N10mc./exp(-l10.*tv1(i));
    N26t1mc = N26mc./exp(-l26.*tv1(i));
    
    % check and exchange saturated N10 and N26
    N10t1 = N10t1.*(N10t1 <= P10./l10.*(1-exp(-l10.*1e7))) + ...
        (N10t1 > P10./l10.*(1-exp(-l10.*1e7))) .* P10./l10.*(1-exp(-l10.*1e7));
    N26t1 = N26t1.*(N26t1 <= P26./l26.*(1-exp(-l26.*1e7))) + ...
        (N26t1 > P26./l26.*(1-exp(-l26.*1e7))) .* P26./l26.*(1-exp(-l26.*1e7));
    N10t1mc = N10t1mc.*(N10t1mc <= P10mc./l10.*(1-exp(-l10.*1e7))) + ...
        (N10t1mc > P10mc./l10.*(1-exp(-l10.*1e7))) .* P10mc./l10.*(1-exp(-l10.*1e7));
    N26t1mc = N26t1mc.*(N26t1mc <= P26mc./l26.*(1-exp(-l26.*1e7))) + ...
        (N26t1mc > P26mc./l26.*(1-exp(-l26.*1e7))) .* P26mc./l26.*(1-exp(-l26.*1e7));
    
    % calculate simple exposure N26 assuming simple exposure N10t1
    N26t2 = P26./l26.*(1-exp(-l26.*log(1-N10t1.*l10./P10)./-l10));
    N26t2mc = P26mc./l26.*(1-exp(-l26.*log(1-N10t1mc.*l10./P10mc)./-l10));
    
    % test if N26 is too high (too long burial)
    burtest1 = burtest1 + (N26t1 < N26t2);
    burtest1mc = burtest1mc + (N26t1mc < N26t2mc);
end;

tv2 = 0;
tv2mc = zeros(mc,1);
burtest2 = zeros(1,6);
burtest2mc = zeros(mc,6);
for i = 1:6;
    tv2 = (burtest1-11).*1e6 + burtest2(1).*1e5 + burtest2(2).*1e4 + burtest2(3).*1e3 + ...
        burtest2(4).*100 + burtest2(5).*10 + burtest2(6); % time vector to start with
    tv2mc = (burtest1mc-11).*1e6 + burtest2mc(:,1).*1e5 + burtest2mc(:,2).*1e4 + ...
        burtest2mc(:,3).*1e3 + burtest2mc(:,4).*100 + burtest2mc(:,5).*10 + burtest2mc(:,6);
    
    for j = 1:10;
        % calculate N before burial
        N10t1 = N10./exp(-l10.*tv2);
        N26t1 = N26./exp(-l26.*tv2);
        N10t1mc = N10mc./exp(-l10.*tv2mc);
        N26t1mc = N26mc./exp(-l26.*tv2mc);
        
        % check and exchange saturated N10 and N26
        N10t1 = N10t1.*(N10t1 <= P10./l10.*(1-exp(-l10.*1e7))) + ...
            (N10t1 > P10./l10.*(1-exp(-l10.*1e7))) .* P10./l10.*(1-exp(-l10.*1e7));
        N26t1 = N26t1.*(N26t1 <= P26./l26.*(1-exp(-l26.*1e7))) + ...
            (N26t1 > P26./l26.*(1-exp(-l26.*1e7))) .* P26./l26.*(1-exp(-l26.*1e7));
        N10t1mc = N10t1mc.*(N10t1mc <= P10mc./l10.*(1-exp(-l10.*1e7))) + ...
            (N10t1mc > P10mc./l10.*(1-exp(-l10.*1e7))) .* P10mc./l10.*(1-exp(-l10.*1e7));
        N26t1mc = N26t1mc.*(N26t1mc <= P26mc./l26.*(1-exp(-l26.*1e7))) + ...
            (N26t1mc > P26mc./l26.*(1-exp(-l26.*1e7))) .* P26mc./l26.*(1-exp(-l26.*1e7));
        
        % calculate simple exposure N26 assuming simple exposure N10t1
        N26t2 = P26./l26.*(1-exp(-l26.*log(1-N10t1.*l10./P10)./-l10));
        N26t2mc = P26mc./l26.*(1-exp(-l26.*log(1-N10t1mc.*l10./P10mc)./-l10));
        
        % test if N26 is too high (too long burial)
        burtest2(1,i) = burtest2(1,i) + (N26t1 < N26t2);
        burtest2mc(:,i) = burtest2mc(:,i) + (N26t1mc < N26t2mc);
        
        % add time step to tv2
        tv2 = tv2 + 10^(6-i);
        tv2mc = tv2mc + 10^(6-i);
    end;
    burtest2(1,i) = burtest2(1,i)-1;
    burtest2mc(:,i) = burtest2mc(:,i)-1;
    clear tv2;
    clear tv2mc;
end;

% calculate burial ages
burialt = (burtest1-11).*1e6 + burtest2(1).*1e5 + burtest2(2).*1e4 + burtest2(3).*1e3 + ...
    burtest2(4).*100 + burtest2(5).*10 + burtest2(6);
burialtmc = (burtest1mc-11).*1e6 + burtest2mc(:,1).*1e5 + burtest2mc(:,2).*1e4 + ...
    burtest2mc(:,3).*1e3 + burtest2mc(:,4).*100 + burtest2mc(:,5).*10 + burtest2mc(:,6);

% calculate N at start of burial
burialN = N10./exp(-l10.*burialt);
burialNmc = N10mc./exp(-l10.*burialtmc);
% check and exchange saturated N10
burialN = burialN.*(burialN <= P10./l10.*(1-exp(-l10.*1e7))) + ...
    (burialN > P10./l10.*(1-exp(-l10.*1e7))) .* P10./l10.*(1-exp(-l10.*1e7));
burialNmc = burialNmc.*(burialNmc <= P10mc./l10.*(1-exp(-l10.*1e7))) + ...
    (burialNmc > P10mc./l10.*(1-exp(-l10.*1e7))) .* P10mc./l10.*(1-exp(-l10.*1e7));
% calculate simple exposure durations
expt = abs(log(1 - burialN.*l10./P10)./-l10);
exptmc = log(1 - burialNmc.*l10./P10mc)./-l10;

% calculate mean and stdev (mc simulation)
burialmc = mean(burialtmc);
del_burialmc = std(burialtmc);
expmc = mean(exptmc);
del_expmc = std(exptmc);

% fix for saturated samples
if round(expt)+del_expmc >= 1E7;
    del_expmc = -1; % signals 10Be saturation
end;

% max burial: 10 Ma
if burialt+del_burialmc >= 1E7;
    del_burialmc = -1; % signals too much burial
    if burialt > 1E7;
        burialt = 1E7;
    end;
end;

% max negative burial: 10 Ma
if burialt-del_burialmc <= -1E7;
    del_burialmc = -2; % signals too much negative burial
    if burialt < -1E7;
        burialt = -1E7;
    end;
end;

% 5. Results structure assignment
results.burial = burialt;
results.del_burial = del_burialmc;
results.expdur = expt;
results.del_expdur = del_expmc;
results.N10norm = N10norm;
results.N10norm_del = N10norm_del;
results.N26norm = N26norm;
results.N26norm_del = N26norm_del;
results.P10sp = P_nu10;
results.P26sp = P_nu26;
results.P_mu10 = P_mu10;
results.P_mu26 = P_mu26;
results.Lsp = Lsp;
results.RcEst = LSDfix.RcEst;
results.SPhiInf = consts.SPhiInf;
% end subfunction get_1026_burial ==================================================================


% plotting =========================================================================================
function plotting(expline,eline,points,sigma1line,sigma2line,probcontours,pl,consts);
% plots banana
set(gcf,'visible','off'); % hide plots
hold on;
box on;

% simple exposure line,
if expline > 0;
    % create data for the simple-exposure line including ratio uncertainties
    tempt = logspace(0,7,100);
    be = (1/consts.l10)*(1-exp(-consts.l10*tempt));
    al = (1/consts.l26)*(1-exp(-consts.l26*tempt));
    al_be = al./be;
    
    albe_P_del = sqrt((consts.delP10_ref_nu./consts.P10_ref_nu)^2 + ...
        (consts.delP26_ref_nu./consts.P26_ref_nu)^2);
    al_be_high = al_be.*(1 + albe_P_del);
    al_be_low = al_be.*(1 - albe_P_del);
    
    % plot ratio uncertainty area and simple exposure line
    patch([be flip(be)],[al_be_high flip(al_be_low)],[0.9 0.9 0.9],'EdgeColor','none');
    semilogx(be,al_be,'color','black');
    
    % make data for simple burial lines
    bu_be50 = be.*exp(-consts.l10*500000);
    bu_al50 = al.*exp(-consts.l26*500000);
    bu_be100 = be.*exp(-consts.l10*1000000);
    bu_al100 = al.*exp(-consts.l26*1000000);
    bu_be150 = be.*exp(-consts.l10*1500000);
    bu_al150 = al.*exp(-consts.l26*1500000);
    bu_be200 = be.*exp(-consts.l10*2000000);
    bu_al200 = al.*exp(-consts.l26*2000000);
    bu_be250 = be.*exp(-consts.l10*2500000);
    bu_al250 = al.*exp(-consts.l26*2500000);
    bu_be300 = be.*exp(-consts.l10*3000000);
    bu_al300 = al.*exp(-consts.l26*3000000);
    
    % plot simple burial lines
    semilogx(bu_be50,bu_al50./bu_be50,'b');
    semilogx(bu_be100,bu_al100./bu_be100,'b');
    semilogx(bu_be150,bu_al150./bu_be150,'b');
    semilogx(bu_be200,bu_al200./bu_be200,'b');
    semilogx(bu_be250,bu_al250./bu_be250,'b');
    semilogx(bu_be300,bu_al300./bu_be300,'b');
    
    % make data for burial pathways
    burv = (0:1E5:1E7);
    path_be10 = (1/consts.l10)*(1-exp(-consts.l10*1E4)).*exp(-consts.l10.*burv);
    path_be100 = (1/consts.l10)*(1-exp(-consts.l10*1E5)).*exp(-consts.l10.*burv);
    path_be1000 = (1/consts.l10)*(1-exp(-consts.l10*1E6)).*exp(-consts.l10.*burv);
    path_be10000 = (1/consts.l10)*(1-exp(-consts.l10*1E7)).*exp(-consts.l10.*burv);
    path_al10 = (1/consts.l26)*(1-exp(-consts.l26*1E4)).*exp(-consts.l26.*burv);
    path_al100 = (1/consts.l26)*(1-exp(-consts.l26*1E5)).*exp(-consts.l26.*burv);
    path_al1000 = (1/consts.l26)*(1-exp(-consts.l26*1E6)).*exp(-consts.l26.*burv);
    path_al10000 = (1/consts.l26)*(1-exp(-consts.l26*1E7)).*exp(-consts.l26.*burv);
    
    % plot burial pathways
    semilogx(path_be10,path_al10./path_be10,'linestyle','--','color','blue');
    semilogx(path_be100,path_al100./path_be100,'linestyle','--','color','blue');
    semilogx(path_be1000,path_al1000./path_be1000,'linestyle','--','color','blue');
    semilogx(path_be10000,path_al10000./path_be10000,'linestyle','--','color','blue');
end;

% erosion line
if eline > 0;
    fprintf(1,'calculating erosion end-point line...');
    % Precompute P_mu(z) to ~200,000 g/cm2 
    % start at the average sample mid-depth.
    z_sp = [0 logspace(0,5.3,100)];
    z_mu = z_sp+(mean(pl.thick_rho)./2);
    P_mu_z = P_mu_LSD(z_mu,mean(pl.atm),mean(pl.RcEst),consts.SPhiInf,1,1,consts,'no');
    P_mu_z10 = P_mu_z.Be.*mean(pl.othercorr);
    P_mu_z26 = P_mu_z.Al.*mean(pl.othercorr);
    P10sp_z = mean(pl.P10sp).*exp(-z_sp./mean(pl.Lsp));
    P26sp_z = mean(pl.P26sp).*exp(-z_sp./mean(pl.Lsp));
    
    tempe = logspace(-5,1,100);
    for j = 1:numel(tempe);
        tv = z_sp./tempe(j);
        dcf10 = exp(-tv.*consts.l10);
        dcf26 = exp(-tv.*consts.l26);
        N10tv = cumtrapz(tv,(P10sp_z.*dcf10 + P_mu_z10.*dcf10));
        N26tv = cumtrapz(tv,(P26sp_z.*dcf26 + P_mu_z26.*dcf26));
        bee(j+1,1) = N10tv(end)./(mean(pl.P10sp) + mean(pl.P_mu10));
        ale(j+1,1) = N26tv(end)./(mean(pl.P26sp) + mean(pl.P_mu26));
    end;
    bee(1,1) = be(end);
    ale(1,1) = al(end);
    ale_bee = ale./bee;
    
    semilogx(bee,ale_bee,'linestyle','--','color','black');
    fprintf(1,' done!\n');
end;

% probability contours
if probcontours > 0;
    % Estimate range and create mesh
    Rv = pl.N26n./pl.N10n;
    delRv = sqrt((pl.delN26n./pl.N26n).^2 + (pl.delN10n./pl.N10n).^2);
    
    xmin = min(pl.N10n - 4.*pl.delN10n);
    xmax = max(pl.N10n + 4.*pl.delN10n);
    Rmin = min(Rv.*(1 - 4.*delRv));
    Rmax = max(Rv.*(1 + 4.*delRv));
    
    [xa,ya] = meshgrid(xmin:0.1 * mean(pl.delN10n):xmax,Rmin:0.1.*mean(Rv).*mean(delRv):Rmax);
    
    Proba = zeros(size(xa));
    
    for j = 1:numel(pl.N10n(:,1));
        % calculate PDF
        Proba1 = xa.*exp(-0.5.*(((ya.*xa-pl.N26n(j))./pl.delN26n(j)).^2+((xa-pl.N10n(j))./...
            pl.delN10n(j)).^2));
        Proba = Proba + Proba1;
    end;
    
    contour(xa,ya,Proba);
end;

% sigma lines
if sigma1line+sigma2line > 0;
    % loop for each sample with al/be data
    for j = 1:numel(pl.N10n);
        % Estimate range and create mesh
        R = (pl.N26n(j)/pl.N10n(j));
        delR = sqrt((pl.delN26n(j) / pl.N26n(j))^2 + (pl.delN10n(j) / pl.N10n(j))^2);
        
        [x,y] = meshgrid((pl.N10n(j) - 4*pl.delN10n(j)):(0.1*pl.delN10n(j)):(pl.N10n(j) +...
            4*pl.delN10n(j)),(R*(1 - 4*delR)):(0.1*R*delR):(R*(1 + 4*delR)));
        
        % calculate PDF
        Prob = x.*exp(-0.5.*((((y.*x) - pl.N26n(j))./pl.delN26n(j)).^2 + ...
            ((x - pl.N10n(j))./pl.delN10n(j)).^2));
        
        % Now find the 68% probability contour.
        % Normalize to volume = 1
        normP = Prob ./ sum(sum(Prob));
        
        % Now we need to figure out cumulative probabilities.
        % multiply by 10000 to achieve manageable number of values, round to get integers:
        normP = normP * 10000;
        intP = round(normP);
        
        for a = 1:max(max(intP));
            cumprob(a) = sum(intP(find(intP >= a)))/10000;
        end;
        
        probs = 1:max(max(intP));
        sigma1 = find(abs(cumprob - 0.6827) == min(abs(cumprob - 0.6827)));
        sigma2 = find(abs(cumprob - 0.9545) == min(abs(cumprob - 0.9545)));
        
        % weed out cases where adjacent probs are same -- rounding error
        if length(sigma1) ~= 1;
            sigma1 = min(sigma1);
        end;
        if length(sigma2) ~= 1;
            sigma2 = min(sigma2);
        end;
        
        % Now draw the contours.
        cmat = contourc(x(1,:),y(:,1),normP,[sigma1 sigma2]);
        
        if sigma1line > 0;
            % Sometimes contourc returns several contours for one level - grid size issue?
            % This is spurious, so plot only the major one.
            contourStarts = find(cmat(1,:) == sigma1);
            contourSizes = cmat(2,contourStarts);
            contourToPlot = find(contourSizes == max(contourSizes));
            
            x1 = cmat(1,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + ...
                contourSizes(contourToPlot)));
            y1 = cmat(2,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + ...
                contourSizes(contourToPlot)));
            
            % plot 1 sigma uncertainty
            sd1 = plot(x1,y1,'color','red');
        end;
        
        if sigma2line > 0;
            % for sigma2 lines
            contourStarts = find(cmat(1,:) == sigma2);
            contourSizes = cmat(2,contourStarts);
            contourToPlot = find(contourSizes == max(contourSizes));
            
            x2 = cmat(1,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + ...
                contourSizes(contourToPlot)));
            y2 = cmat(2,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + ...
                contourSizes(contourToPlot)));
            
            % plot 2 sigma uncertainty
            sd2 = plot(x2,y2,'color','red');
        end;
    end;
end;

% sample points
if points > 0;
    for j = 1:numel(pl.N10n);
        plot(pl.N10n(j),pl.N26n(j)/pl.N10n(j),'.','color','red','markersize',15);
    end;
end;

% fix and display plot
% fix axis
axis([1000 3000000 0.2 1.2]);
xlabel('[^{10}Be]*');
ylabel('[^{26}Al]*/[^{10}Be]*');

set(gca,'layer','top'); % plot axis on top
set(gca,'XScale','log'); % fix for matlab
set(gcf,'visible','on'); % display plot
hold off;
% end subfunction plotting =========================================================================
