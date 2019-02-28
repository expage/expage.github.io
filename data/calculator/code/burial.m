function burial();

% expage 10Be and 26Al exposure/burial duration calculator.
% Read and fix input, use get_exposure_burial to calculate exposure and burial durations and/or
% get_erosion_burial to calculate erosion rate and burial duration, and fix and write output. If
% plotting is on, display P-normalized 26/10 banana plot.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

clear;
close all;
tic();

% make choices here ================================================================================
% exposure + burial and/or erosion + burial? (1 = yes)
exposure_burial = 1; % calculate burial assuming prior surface exposure
erosion_burial = 0; % calculate burial assuming prior constant erosion

% calculate external and/or internal uncertainties? (1 = yes)
extunc = 0;
intunc = 1;
% Monte Carlo size for uncertainty estimate
mc = 1E5;

% plotting? (1 = yes)
pl.expline = 1; % plot simple exposure line (and burial paths/lines)
pl.eline = 1; % plot erosion end point line (and burial paths/lines)
pl.points = 1; % plot sample points
pl.sigma1line = 1; % plot sample 1 sigma line
pl.sigma2line = 0; % plot sample 2 sigma line
pl.probcontours = 0; % plot probability contour lines

% write normalized conc data? (for later plotting)
normNout = 1;
% ==================================================================================================

% What version is this?
ver = '201902';

% count number of input columns in line 1 of input
inid = fopen('input.txt','r');
line1 = fgets(inid);
line1 = regexprep(line1,' +',' ');
numcols = numel(strfind(line1,' ')) + numel(strfind(line1,sprintf('\t',''))) + 1;
fclose(inid);

% read input file
if numcols == 16; % if erosion rate in input
    [samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
        samplein.rho,samplein.othercorr,samplein.E,samplein.N10,samplein.delN10,samplein.be_stds,...
        samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr] = ...
        textread('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n','commentstyle',...
        'matlab');
else; % has to be 15 columns!
    [samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
        samplein.rho,samplein.othercorr,samplein.N10,samplein.delN10,samplein.be_stds,...
        samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr] = ...
        textread('input.txt','%s %n %n %n %s %n %n %n %n %n %s %n %n %s %n','commentstyle',...
        'matlab');
end;

% run and load expage constants
make_consts_expage;
load consts_expage;

% display external/internal uncertainty info
if extunc == 1 && intunc == 1;
    fprintf(1,'external and (internal) uncertainties\n');
elseif extunc == 1;
    fprintf(1,'external uncertainties\n');
elseif intunc == 1;
    fprintf(1,'internal uncertainties\n');
end;

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

% constants
Pref10 = consts.P10_ref_nu; delPref10 = consts.delP10_ref_nu;
Pref26 = consts.P26_ref_nu; delPref26 = consts.delP26_ref_nu;
% decay constants
l10 = consts.l10; dell10 = consts.dell10;
l26 = consts.l26; dell26 = consts.dell26;

% declare pl parameters
pl.N10n = []; pl.delN10n = []; pl.N26n = []; pl.delN26n = [];
pl.P10sp = []; pl.P26sp = []; pl.P_mu10 = []; pl.P_mu26 = []; pl.P10z = []; pl.P26z = [];
pl.Lsp = []; pl.RcEst = []; pl.SPhiAv = []; pl.thick_rho = []; pl.atm = []; pl.othercorr = [];

% fix for output
output(1,1) = {'sample'};
outn(1) = 1;
if sum(samplein.N10+samplein.N26) > 0;
    if exposure_burial == 1;
        output(1,end+1) = {'exposure(yr)'}; outn(1) = max(outn)+1;
        if extunc == 1; output(1,end+1) = {'exp-uncext(yr)'}; outn(2) = max(outn)+1; end;
        if intunc == 1; output(1,end+1) = {'exp-uncint(yr)'}; outn(3) = max(outn)+1; end;
        output(1,end+1) = {'burial(yr)'}; outn(4) = max(outn)+1;
        if extunc == 1; output(1,end+1) = {'bur-uncext(yr)'}; outn(5) = max(outn)+1; end;
        if intunc == 1; output(1,end+1) = {'bur-uncint(yr)'}; outn(6) = max(outn)+1; end;
    end;
    if erosion_burial == 1;
        output(1,end+1) = {'erosion(mm/ka)'}; outn(7) = max(outn)+1;
        if extunc == 1; output(1,end+1) = {'ero-uncext(mm/ka)'}; outn(8) = max(outn)+1; end;
        if intunc == 1; output(1,end+1) = {'ero-uncint(mm/ka)'}; outn(9) = max(outn)+1; end;
        output(1,end+1) = {'burial(yr)'}; outn(10) = max(outn)+1;
        if extunc == 1; output(1,end+1) = {'bur-uncext(yr)'}; outn(11) = max(outn)+1; end;
        if intunc == 1; output(1,end+1) = {'bur-uncint(yr)'}; outn(12) = max(outn)+1; end;
    end;
    if normNout == 1; % if writing normalized conc data
        output(1,end+1) = {'N10norm'}; outn(13) = max(outn)+1;
        output(1,end+1) = {'N10unc'}; outn(14) = max(outn)+1;
        output(1,end+1) = {'N26norm'}; outn(15) = max(outn)+1;
        output(1,end+1) = {'N26unc'}; outn(16) = max(outn)+1;
    end;
end;

% pick out samples one by one
for i = 1:numel(samplein.lat);
    sample.sample_name = samplein.sample_name(i);
    sample.lat = samplein.lat(i);
    sample.long = samplein.long(i);
    sample.pressure = samplein.pressure(i);
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
    
    % fix Rc SPhi and tv
    LSDfix = LSD_fix(sample.lat,sample.long,1E7,-1,consts);
    
    % interpolate Lsp (Sato, 2008; Marrero et al., 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,LSDfix.RcEst);
    
    % thickness scaling factor.
    if sample.thick > 0;
        thickSF = (sample.Lsp./(sample.rho.*sample.thick)).* ...
            (1-exp(((-1.*sample.rho.*sample.thick)./sample.Lsp)));
    else 
        thickSF = 1;
    end;
    
    % if calculating uncertainties
    if extunc+intunc >= 1;
        % randomized 10Be and 26Al concentrations for mc simulation
        sample.N10mc = normrnd(sample.N10,sample.delN10,[1 mc]);
        sample.N26mc = normrnd(sample.N26,sample.delN26,[1 mc]);
        % too low N not allowed!
        sample.N10mc(sample.N10mc<1) = 1;
        sample.N26mc(sample.N26mc<1) = 1;
        % randomized decay constants for mc simulation
        sample.l10mc = normrnd(l10,dell10,[1 mc]);
        sample.l26mc = normrnd(l26,dell26,[1 mc]);
    end;
    
    % spallation production scaling
    Psp = LSDspal(sample.pressure,LSDfix.RcEst,consts.SPhiInf,LSDfix.w,1,1,consts);
    
    % randomized ref prod rates for Monte Carlo simulation
    Pref10mc = normrnd(Pref10,delPref10,[1 mc]);
    Pref26mc = normrnd(Pref26,delPref26,[1 mc]);
    
    % spallation prod rates
    sample.Psp10 = Psp.sp10.*Pref10.*thickSF.*sample.othercorr;
    sample.Psp26 = Psp.sp26.*Pref26.*thickSF.*sample.othercorr;
    sample.Psp10mc = Psp.sp10.*Pref10mc.*thickSF.*sample.othercorr;
    sample.Psp26mc = Psp.sp26.*Pref26mc.*thickSF.*sample.othercorr;
    
    % fix shielding depth vector z (g/cm2/yr)
    if erosion_burial == 1;
        sample.z = [0 logspace(0,5.3,100)]';
    else;
        sample.z = 0;
    end;
    zmu = sample.z + sample.thick.*sample.rho./2; % add half sample depth
    
    % muon production
    P_mu = P_mu_expage(zmu,sample.pressure,LSDfix.RcEst,consts.SPhiInf,1,1,consts,'no');
    sample.Pmu10 = P_mu.mu10'.*sample.othercorr;
    sample.Pmu26 = P_mu.mu26'.*sample.othercorr;
    
    % normalized N (sent out in results for plotting)
    N10norm = sample.N10./(sample.Psp10+sample.Pmu10(1));
    N10norm_del = sample.delN10./(sample.Psp10+sample.Pmu10(1));
    N26norm = sample.N26./(sample.Psp26+sample.Pmu26(1));
    N26norm_del = sample.delN26./(sample.Psp26+sample.Pmu26(1));
    
    % calculate exposure + burial
    if exposure_burial == 1;
        expbur = get_exposure_burial(sample,consts,mc,extunc,intunc);
        
        % write to output
        output(i+1,outn(1)) = {num2str(expbur.exposure,'%.0f')};
        output(i+1,outn(4)) = {num2str(expbur.burial,'%.0f')};
        if extunc == 1;
            output(i+1,outn(2)) = {num2str(expbur.exposure_uncext,'%.0f')};
            output(i+1,outn(5)) = {num2str(expbur.burial_uncext,'%.0f')};
        end;
        if intunc == 1;
            output(i+1,outn(3)) = {num2str(expbur.exposure_uncint,'%.0f')};
            output(i+1,outn(6)) = {num2str(expbur.burial_uncint,'%.0f')};
        end;
        
        % display exposure/burial and fix output if no solution
        if erosion_burial == 1; fprintf(1,'\n'); end; % separate lines if both exposure/erosion
        if expbur.exposure >= 0; % if initial exposure within bounds
            if extunc == 1 && intunc == 1;
                fprintf(1,' \texposure = %s ± %s (%s) yr',output{i+1,outn(1)},...
                    output{i+1,outn(2)},output{i+1,outn(3)});
                fprintf(1,' \tburial = %s ± %s (%s) yr',output{i+1,outn(4)},output{i+1,outn(5)},...
                    output{i+1,outn(6)});
            elseif extunc == 1;
                fprintf(1,' \texposure = %s ± %s yr',output{i+1,outn(1)},output{i+1,outn(2)});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outn(4)},output{i+1,outn(5)});
            elseif intunc == 1;
                fprintf(1,' \texposure = %s ± %s yr',output{i+1,outn(1)},output{i+1,outn(3)});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outn(4)},output{i+1,outn(6)});
            else;
                fprintf(1,' \texposure = %s yr',output{i+1,outn(1)});
                fprintf(1,' \tburial = %s yr',output{i+1,outn(4)});
            end;
        else;
            fprintf(1,' \tno solution!');
        end;
    end;
    
    % calculate erosion + burial
    if erosion_burial == 1;
        erobur = get_erosion_burial(sample,consts,mc,extunc,intunc);
        
        % write to output
        if erobur.erosion>0.1;
            output(i+1,outn(7)) = {num2str(erobur.erosion,'%.2f')};
            if extunc == 1;
                output(i+1,outn(8)) = {num2str(erobur.erosion_uncext,'%.2f')};
            end;
            if intunc == 1;
                output(i+1,outn(9)) = {num2str(erobur.erosion_uncint,'%.2f')};
            end;
        else;
            output(i+1,outn(7)) = {num2str(erobur.erosion,'%.3f')};
            if extunc == 1;
                output(i+1,outn(8)) = {num2str(erobur.erosion_uncext,'%.3f')};
            end;
            if intunc == 1;
                output(i+1,outn(9)) = {num2str(erobur.erosion_uncint,'%.3f')};
            end;
        end;
        output(i+1,outn(10)) = {num2str(erobur.burial,'%.0f')};
        if extunc == 1;
            output(i+1,outn(11)) = {num2str(erobur.burial_uncext,'%.0f')};
        end;
        if intunc == 1;
            output(i+1,outn(12)) = {num2str(erobur.burial_uncint,'%.0f')};
        end;
        
        % display erosion/burial and fix output if no solution
        if exposure_burial == 1; fprintf(1,'\n'); end; % separate lines if both exposure/erosion
        if erobur.erosion >= 0; % if initial erosion within bounds
            if extunc == 1 && intunc == 1;
                fprintf(1,' \terosion = %s ± %s (%s) mm/ka',output{i+1,outn(7)},...
                    output{i+1,outn(8)},output{i+1,outn(9)});
                fprintf(1,' \tburial = %s ± %s (%s) yr',output{i+1,outn(10)},...
                    output{i+1,outn(11)},output{i+1,outn(12)});
            elseif extunc == 1;
                fprintf(1,' \terosion = %s ± %s mm/ka',output{i+1,outn(7)},output{i+1,outn(8)});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outn(10)},output{i+1,outn(11)});
            elseif intunc == 1;
                fprintf(1,' \terosion = %s ± %s mm/ka',output{i+1,outn(7)},output{i+1,outn(9)});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outn(10)},output{i+1,outn(12)});
            else;
                fprintf(1,' \terosion = %s mm/ka',output{i+1,outn(7)});
                fprintf(1,' \tburial = %s yr',output{i+1,outn(10)});
            end;
        else;
            fprintf(1,' \tno solution!');
        end;
    end;
    
    % if writing normalized conc data
    if normNout > 0;
        output(i+1,outn(13)) = {num2str(N10norm,'%.6f')};
        output(i+1,outn(14)) = {num2str(N10norm_del,'%.6f')};
        output(i+1,outn(15)) = {num2str(N26norm,'%.6f')};
        output(i+1,outn(16)) = {num2str(N26norm_del,'%.6f')};
    end;
    
    % new line
    fprintf(1,'\n');
    
    % for plotting
    pl.N10n(end+1,1) = N10norm;
    pl.delN10n(end+1,1) = N10norm_del;
    pl.N26n(end+1,1) = N26norm;
    pl.delN26n(end+1,1) = N26norm_del;
    pl.P10sp(end+1,1) = sample.Psp10;
    pl.P26sp(end+1,1) = sample.Psp26;
    pl.P_mu10(end+1,1) = sample.Pmu10(1);
    pl.P_mu26(end+1,1) = sample.Pmu26(1);
    pl.Lsp(end+1,1) = sample.Lsp;
    pl.RcEst(end+1,1) = LSDfix.RcEst;
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
if (pl.expline + pl.eline + pl.points + pl.sigma1line + pl.probcontours) > 0;
    pl.exposure_burial = exposure_burial;
    pl.erosion_burial = erosion_burial;
    plotting(pl,consts);
end;
% end plotting ====================================================================

toc()
% end burial function ==============================================================================


% subfunction get_exposure_burial ==================================================================
function out = get_exposure_burial(sample,consts,mc,extunc,intunc);
% This function calculates the simple albe exposure and burial durations

% number of elements in exposure vector
tN = 5;
% limit for while loop
burlim = 10;

% atoms/g measurement
N10 = sample.N10;
N26 = sample.N26;
N10mc = repmat(sample.N10mc,tN,1);
N26mc = repmat(sample.N26mc,tN,1);
% decay constants
l10 = consts.l10;
l26 = consts.l26;
% summed prod rates
P10 = sample.Psp10 + sample.Pmu10(1);
P26 = sample.Psp26 + sample.Pmu26(1);
% min and max t to start from
mint = 1;
maxt = 1E7;
mintmc(1:mc) = mint;
maxtmc(1:mc) = maxt;

% loop for central point
burtest = burlim+1; % parameter for while loop
while min(abs(burtest))>burlim;
    % fix time vector
    expv = linspace(mint,maxt,tN)';
    
    % calculate N after expv
    N10exp = P10./l10.*(1-exp(-l10.*expv));
    N26exp = P26./l26.*(1-exp(-l26.*expv));
    
    % calculate burial needed to reach N10exp/N26exp
    bur10 = log(N10./N10exp)./-l10;
    bur26 = log(N26./N26exp)./-l26;
    
    % calculate difference in burial duration for N10 and N26
    burtest = bur26-bur10;
    
    % if solution not possible
    if max(burtest)<0 | min(burtest)>0;
        out.exposure = -1; % signals no solution
        out.burial = -1;
        break;
    end;
    
    % remove negative values
    burtestpos = burtest.*(burtest>=0) + (burtest<0).*max(abs(burtest));
    
    % find index for min burtest > 0
    [m,idx] = min(burtestpos);
    
    % set new min and max e
    mint = expv(idx);
    maxt = expv(min(idx+1,tN));
end;
if min(abs(burtest))<=burlim;
    % pick out exposure and burial
    out.exposure = interp1(burtest,expv,0);
    out.burial = interp1(burtest,bur10,0);
end;

if extunc+intunc >= 1;
    N10mc = repmat(sample.N10mc,tN,1);
    N26mc = repmat(sample.N26mc,tN,1);
end;

% if calculating external uncertainties
if extunc == 1;
    % if solution exist
    if out.exposure >= 0;
        % decay constants
        l10mc = repmat(sample.l10mc,tN,1);
        l26mc = repmat(sample.l26mc,tN,1);
        % summed prod rates
        P10mc = repmat(sample.Psp10mc + sample.Pmu10(1),tN,1);
        P26mc = repmat(sample.Psp26mc + sample.Pmu26(1),tN,1);
        
        % calculations in mc subfunction
        [out.exposure_uncext out.burial_uncext] = expbur_mc(N10mc,N26mc,l10mc,l26mc,P10mc,P26mc,...
            mintmc,maxtmc,mc,tN,burlim,out.exposure,out.burial);
    else;
        out.exposure_uncext = -1;
        out.burial_uncext = -1;
    end;
end;

% if calculating internal uncertainties
if intunc == 1;
    % if solution exist
    if out.exposure >= 0;
        % calculations in mc subfunction
        [out.exposure_uncint out.burial_uncint] = expbur_mc(N10mc,N26mc,l10,l26,P10,P26,...
            mintmc,maxtmc,mc,tN,burlim,out.exposure,out.burial);
    else;
        out.exposure_uncint = -1;
        out.burial_uncint = -1;
    end;
end;
% end subfunction get_exposure_burial ==============================================================


% subfunction expbur_mc ============================================================================
function [exposure_unc burial_unc] = ...
    expbur_mc(N10mc,N26mc,l10mc,l26mc,P10mc,P26mc,mintmc,maxtmc,mc,tN,burlim,exp_mid,bur_mid);
% this function calculates the mc uncertainties for exposure+burial
burtestmc = burlim+1; % parameter for while loop
loopN = 1; % for first round in while loop
while mean(min(abs(burtestmc)))>burlim;
    % fix time matrix
    %expvmc = linspace(mintmc,maxtmc,tN)'; % does not work in Matlab...
    expvmc = repmat((0:tN-1)',size(mintmc)).*repmat((maxtmc-mintmc)./(tN-1),tN,1) + ...
        repmat(mintmc,tN,1);
    
    % calculate N after expv
    N10expmc = P10mc./l10mc.*(1-exp(-l10mc.*expvmc));
    N26expmc = P26mc./l26mc.*(1-exp(-l26mc.*expvmc));
    
    % calculate burial needed to reach N10exp/N26exp
    bur10mc = log(N10mc./N10expmc)./-l10mc;
    bur26mc = log(N26mc./N26expmc)./-l26mc;
    
    % calculate difference in burial duration for N10 and N26
    burtestmc = bur26mc-bur10mc;
    
    % only in first round
    if loopN == 1;
        % remove mc points without solution
        mcrm = find((max(burtestmc)<0)|(min(burtestmc)>0));
        mc = mc-numel(mcrm);
        burtestmc(:,mcrm) = [];
        bur10mc(:,mcrm) = []; bur26mc(:,mcrm) = [];
        N10mc(:,mcrm) = []; N26mc(:,mcrm) = [];
        expvmc(:,mcrm) = [];
        if numel(l10mc) > 1; % if external uncertainties
            P10mc(:,mcrm) = []; P26mc(:,mcrm) = [];
            l10mc(:,mcrm) = []; l26mc(:,mcrm) = [];
        end;
        loopN = 2;
    end;
    
    % remove negative values
    burtestmcpos = burtestmc.*(burtestmc>=0) + (burtestmc<0).*max(max(abs(burtestmc)));
    
    % find index for min burtest > 0
    [mmc,idxmc] = min(burtestmcpos);
    idxmc1 = sub2ind(size(burtestmc),idxmc,(1:1:mc));
    idxmc2 = sub2ind(size(burtestmc),min(idxmc+1,tN),(1:1:mc));
    
    % set new min and max e
    mintmc = expvmc(idxmc1);
    maxtmc = expvmc(idxmc2);
end;

% pick out mc exposure and burial
exposuremc = -burtestmc(idxmc2)./(-burtestmc(idxmc2)+burtestmc(idxmc1)).*...
    (expvmc(idxmc1)-expvmc(idxmc2))+expvmc(idxmc2);
burialmc = -burtestmc(idxmc2)./(-burtestmc(idxmc2)+burtestmc(idxmc1)).*...
    (bur10mc(idxmc1)-bur10mc(idxmc2))+bur10mc(idxmc2);
% uncertainty based on mc deviation from exposure and burial
exposure_unc = sqrt(sum((exposuremc-exp_mid).^2)./mc);
burial_unc = sqrt(sum((burialmc-bur_mid).^2)./mc);
% end subfunction expbur_mc ========================================================================


% subfunction get_erosion_burial ===================================================================
function out = get_erosion_burial(sample,consts,mc,extunc,intunc);
% this function calculates the simple albe erosion rate and burial duration

% number of elements in erosion vector
eN = 5;
% limit for while loop (yr)
burlim = 10;

% atoms/g measurement
N10 = sample.N10;
N26 = sample.N26;
N10mc = repmat(sample.N10mc,eN,1);
N26mc = repmat(sample.N26mc,eN,1);
% decay constants
l10 = consts.l10;
l26 = consts.l26;
% fix prod rates
Psp10 = sample.Psp10.*exp(-sample.z./sample.Lsp);
Psp26 = sample.Psp26.*exp(-sample.z./sample.Lsp);
Pmu10 = sample.Pmu10;
Pmu26 = sample.Pmu26;
P10 = Psp10 + Pmu10;
P26 = Psp26 + Pmu26;
% min and max e to start from (g/cm2/yr)
mine = 0;
maxe = 10;
minemc(1:mc) = 0;
maxemc(1:mc) = 10;

% loop for central point
burtest = burlim+1; % parameter for while loop
while min(abs(burtest))>burlim;
    % fix erosion vectors
    if mine > 0;
        erov = logspace(log10(maxe),log10(mine),eN);
    else;
        erov = [logspace(log10(maxe),log10(maxe.*1E-4),eN-1) 0];
    end;
    
    % fix for zero erosion
    erovfix = (erov>0).*erov + (erov==0);
    
    % calculate N with continuous erosion rate erov
    tv = repmat(sample.z,size(erovfix))./repmat(erovfix,size(sample.z));
    idx0 = find(erov==0);
    tv(:,idx0) = repmat(linspace(0,1E7,numel(sample.z))',size(idx0));
    P10fix = repmat(P10,size(erov));
    P26fix = repmat(P26,size(erov));
    P10fix(:,idx0) = repmat(P10(1),numel(sample.z),numel(idx0));
    P26fix(:,idx0) = repmat(P26(1),numel(sample.z),numel(idx0));
    dcf10 = exp(-tv.*l10);
    dcf26 = exp(-tv.*l26);
    N10ero = trapz_m(tv,P10fix.*dcf10);
    N26ero = trapz_m(tv,P26fix.*dcf26);
    
    % calculate burial needed to reach N10ero/N26ero
    bur10 = log(N10./N10ero)./-l10;
    bur26 = log(N26./N26ero)./-l26;
    
    % calculate difference in burial duration for N10 and N26
    burtest = bur26-bur10;
    
    % if solution not possible
    if max(burtest)<0 | min(burtest)>0;
        out.erosion = -1; % signals no solution
        out.burial = -1;
        break;
    end;
    
    % remove negative values
    burtestpos = burtest.*(burtest>=0) + (burtest<0).*max(abs(burtest));
    
    % find index for min burtest > 0
    [m,idx] = min(burtestpos);
    
    % set new min and max e
    maxe = erov(idx);
    mine = erov(min(idx+1,eN));
end;
if min(abs(burtest))<=burlim;
    % pick out erosion and burial
    out.erosion = interp1(burtest,erov,0);
    out.burial = interp1(burtest,bur10,0);
end;

% if calculating external uncertainties
if extunc == 1;
    % if solution exist
    if out.erosion >= 0;
        % decay constants
        l10mce = repmat(sample.l10mc,eN,1);
        l26mce = repmat(sample.l26mc,eN,1);
        % fix prod rates
        Psp10mce = repmat(sample.Psp10mc,size(sample.z)) .* repmat(exp(-sample.z./sample.Lsp),1,mc);
        Psp26mce = repmat(sample.Psp26mc,size(sample.z)) .* repmat(exp(-sample.z./sample.Lsp),1,mc);
        Pmu10mce = repmat(Pmu10,1,mc);
        Pmu26mce = repmat(Pmu26,1,mc);
        P10mce = Psp10mce + Pmu10mce;
        P26mce = Psp26mce + Pmu26mce;
        
        % calculations in mc subfunction
        [out.erosion_uncext out.burial_uncext] = erobur_mc(sample.z,N10mc,N26mc,l10mce,l26mce,...
            P10mce,P26mce,minemc,maxemc,mc,eN,burlim,out.erosion,out.burial);
        % convert erosion to mm/ka
        out.erosion_uncext = out.erosion_uncext./sample.rho.*1E4;
    else;
        out.erosion_uncext = -1;
        out.burial_uncext = -1;
    end;
end;

% if calculating internal uncertainties
if intunc == 1;
    % if solution exist
    if out.erosion >= 0;
        % decay constants
        l10mci = repmat(l10,eN,mc);
        l26mci = repmat(l26,eN,mc);
        % fix prod rates
        P10mci = repmat(P10,1,mc);
        P26mci = repmat(P26,1,mc);
        
        % calculations in mc subfunction
        [out.erosion_uncint out.burial_uncint] = erobur_mc(sample.z,N10mc,N26mc,l10mci,l26mci,...
            P10mci,P26mci,minemc,maxemc,mc,eN,burlim,out.erosion,out.burial);
        % convert erosion to mm/ka
        out.erosion_uncint = out.erosion_uncint./sample.rho.*1E4;
    else;
        out.erosion_uncint = -1;
        out.burial_uncint = -1;
    end;
end;
% convert erosion to mm/ka
if out.erosion >= 0;
    out.erosion = out.erosion./sample.rho.*1E4;
end;
% end subfunction get_erosion_burial ===============================================================


% subfunction erobur_mc ============================================================================
function [erosion_unc burial_unc] = erobur_mc(z,N10mc,N26mc,l10mc,l26mc,P10mc,P26mc,minemc,maxemc,mc,eN,burlim,ero_mid,bur_mid);
% this function calculates the mc uncertainties for erosion+burial
burtestmc = burlim+1; % parameter for while loop
loopN = 1; % for first round in while loop
while mean(min(abs(burtestmc)))>burlim;
    % fix erosion matrix
    if min(minemc) > 0;
        %erovmc = 1E1.^linspace(log10(maxemc),log10(minemc),eN)'; % does not work in Matlab...
        erovmc = 1E1.^(repmat((0:eN-1)',size(minemc)).*...
            repmat((log10(minemc)-log10(maxemc))./(eN-1),eN,1)+repmat(log10(maxemc),eN,1));
    else;
        mc1idx = find(minemc>0);
        mc0idx = find(minemc==0);
        %erovmc(:,mc1idx) = 1E1.^linspace(log10(maxemc(mc1idx)),log10(minemc(mc1idx)),eN)';
        erovmc(:,mc1idx) = 1E1.^(repmat((0:eN-1)',size(minemc(mc1idx))).*repmat((log10(minemc...
            (mc1idx))-log10(maxemc(mc1idx)))./(eN-1),eN,1)+repmat(log10(maxemc(mc1idx)),eN,1));
        %erovmc(1:end-1,mc0idx) = 1E1.^linspace(log10(maxemc(mc0idx)),...
            %log10(maxemc(mc0idx).*1E-4),eN-1)'; % does not work in Matlab...
        erovmc(1:end-1,mc0idx) = 1E1.^(repmat((0:eN-2)',size(maxemc(mc0idx))).*...
            repmat((log10(maxemc(mc0idx).*1E-4)-log10(maxemc(mc0idx)))./(eN-2),eN-1,1)+...
            repmat(log10(maxemc(mc0idx)),eN-1,1));
        erovmc(end,mc0idx) = 0;
    end;
    
    % fix for zero erosion
    erovmcfix = (erovmc>0).*erovmc + (erovmc==0);
    
    % calculate N with continuous erosion rate erov
    for j = 1:eN;
        tvmc = repmat(z,1,mc)./repmat(erovmcfix(j,:),size(z));
        idx0 = find(erovmc(j,:)==0);
        tvmc(:,idx0) = repmat(linspace(0,1E7,size(z,1))',size(idx0));
        P10mcfix = P10mc;
        P26mcfix = P26mc;
        P10mcfix(:,idx0) = repmat(P10mc(1,idx0),size(z));
        P26mcfix(:,idx0) = repmat(P26mc(1,idx0),size(z));
        dcf10mc = exp(-tvmc.*repmat(l10mc(j,:),size(z)));
        dcf26mc = exp(-tvmc.*repmat(l26mc(j,:),size(z)));
        N10eromc(j,:) = trapz_m(tvmc,P10mcfix.*dcf10mc);
        N26eromc(j,:) = trapz_m(tvmc,P26mcfix.*dcf26mc);
    end;
    
    % calculate burial needed to reach N10ero/N26ero
    bur10mc = log(N10mc./N10eromc)./-l10mc;
    bur26mc = log(N26mc./N26eromc)./-l26mc;
    
    % calculate difference in burial duration for N10 and N26
    burtestmc = bur26mc-bur10mc;
    
    % only in first round
    if loopN == 1;
        % remove mc points without solution
        mcrm = find((max(burtestmc)<0)|(min(burtestmc)>0)==1);
        mc = mc-numel(mcrm);
        burtestmc(:,mcrm) = [];
        bur10mc(:,mcrm) = []; bur26mc(:,mcrm) = [];
        N10mc(:,mcrm) = []; N26mc(:,mcrm) = [];
        clear N10eromc; clear N26eromc;
        P10mc(:,mcrm) = []; P26mc(:,mcrm) = [];
        erovmc(:,mcrm) = [];
        loopN = 2;
    end;
    
    % remove negative values
    burtestmcpos = burtestmc.*(burtestmc>=0) + (burtestmc<0).*max(max(abs(burtestmc)));
    
    % find index for min burtest > 0
    [mmc,idxmc] = min(burtestmcpos);
    idxmc1 = sub2ind(size(burtestmc),idxmc,(1:1:mc));
    idxmc2 = sub2ind(size(burtestmc),min(idxmc+1,eN),(1:1:mc));
    
    % set new min and max e
    maxemc = erovmc(idxmc1);
    minemc = erovmc(idxmc2);
end;

% pick out mc erosion and burial
erosionmc = -burtestmc(idxmc2)./(-burtestmc(idxmc2)+burtestmc(idxmc1)).*...
    (erovmc(idxmc1)-erovmc(idxmc2))+erovmc(idxmc2);
burialmc = -burtestmc(idxmc2)./(-burtestmc(idxmc2)+burtestmc(idxmc1)).*...
    (bur10mc(idxmc1)-bur10mc(idxmc2))+bur10mc(idxmc2);
% uncertainty based on mc deviation from erosion and burial
erosion_unc = sqrt(sum((erosionmc-ero_mid).^2)./mc);
burial_unc = sqrt(sum((burialmc-bur_mid).^2)./mc);
% end subfunction erobur_mc ========================================================================


% plotting =========================================================================================
function plotting(pl,consts);
% plots banana
set(gcf,'visible','off'); % hide plots
hold on;
box on;

% simple exposure line,
if pl.expline > 0;
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
    
    if pl.exposure_burial == 1;
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
        semilogx(bu_be50,bu_al50./bu_be50,'linestyle','--','color','black');
        semilogx(bu_be100,bu_al100./bu_be100,'linestyle','--','color','black');
        semilogx(bu_be150,bu_al150./bu_be150,'linestyle','--','color','black');
        semilogx(bu_be200,bu_al200./bu_be200,'linestyle','--','color','black');
        semilogx(bu_be250,bu_al250./bu_be250,'linestyle','--','color','black');
        semilogx(bu_be300,bu_al300./bu_be300,'linestyle','--','color','black');
        
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
        semilogx(path_be10,path_al10./path_be10,'linestyle','--','color','black');
        semilogx(path_be100,path_al100./path_be100,'linestyle','--','color','black');
        semilogx(path_be1000,path_al1000./path_be1000,'linestyle','--','color','black');
        semilogx(path_be10000,path_al10000./path_be10000,'linestyle','--','color','black');
    end;
end;

% erosion line
if pl.eline > 0;
    fprintf(1,'calculating erosion end-point line...');
    % Precompute P_mu(z) to ~200,000 g/cm2 
    % start at the average sample mid-depth.
    z = [0 logspace(0,5.3,100)]';
    z_mu = z+(mean(pl.thick_rho)./2);
    P_mu_z = P_mu_expage(z_mu,mean(pl.atm),mean(pl.RcEst),consts.SPhiInf,1,1,consts,'no');
    P_mu_z10 = P_mu_z.mu10'.*mean(pl.othercorr);
    P_mu_z26 = P_mu_z.mu26'.*mean(pl.othercorr);
    P10sp_z = mean(pl.P10sp).*exp(-z./mean(pl.Lsp));
    P26sp_z = mean(pl.P26sp).*exp(-z./mean(pl.Lsp));
    P10z = (P_mu_z10+P10sp_z);
    P26z = (P_mu_z26+P26sp_z);
    fprintf(1,' done!\n');
    
    tempe = logspace(-5,1,100);
    zm = repmat(z,size(tempe));
    tvm = zm./repmat(tempe,size(z));
    dcf10 = exp(-tvm.*consts.l10);
    dcf26 = exp(-tvm.*consts.l26);
    N10tv = trapz_m(tvm,repmat(P10z,size(tempe)).*dcf10);
    N26tv = trapz_m(tvm,repmat(P26z,size(tempe)).*dcf26);
    bee = N10tv./P10z(1);
    ale = N26tv./P26z(1);
    bee(1,2:end+1) = bee;
    ale(1,2:end+1) = ale;
    bee(1,1) = (1/consts.l10)*(1-exp(-consts.l10*1E7));
    ale(1,1) = (1/consts.l26)*(1-exp(-consts.l26*1E7));
    ale_bee = ale./bee;
    
    semilogx(bee,ale_bee,'color','blue');
    
    if pl.erosion_burial == 1;
        % make data for simple burial lines
        bu_bee50 = bee.*exp(-consts.l10*500000);
        bu_ale50 = ale.*exp(-consts.l26*500000);
        bu_bee100 = bee.*exp(-consts.l10*1000000);
        bu_ale100 = ale.*exp(-consts.l26*1000000);
        bu_bee150 = bee.*exp(-consts.l10*1500000);
        bu_ale150 = ale.*exp(-consts.l26*1500000);
        bu_bee200 = bee.*exp(-consts.l10*2000000);
        bu_ale200 = ale.*exp(-consts.l26*2000000);
        bu_bee250 = bee.*exp(-consts.l10*2500000);
        bu_ale250 = ale.*exp(-consts.l26*2500000);
        bu_bee300 = bee.*exp(-consts.l10*3000000);
        bu_ale300 = ale.*exp(-consts.l26*3000000);
        
        % plot simple burial lines
        semilogx(bu_bee50,bu_ale50./bu_bee50,'linestyle','--','color','blue');
        semilogx(bu_bee100,bu_ale100./bu_bee100,'linestyle','--','color','blue');
        semilogx(bu_bee150,bu_ale150./bu_bee150,'linestyle','--','color','blue');
        semilogx(bu_bee200,bu_ale200./bu_bee200,'linestyle','--','color','blue');
        semilogx(bu_bee250,bu_ale250./bu_bee250,'linestyle','--','color','blue');
        semilogx(bu_bee300,bu_ale300./bu_bee300,'linestyle','--','color','blue');
        
        % make data for burial pathways
        burv = (0:1E5:1E7);
        pathv = [1E-4 1E-3 1E-2 1E-1];
        zm = repmat(z,size(pathv));
        tvm = zm./repmat(pathv,size(z));
        dcf10 = exp(-tvm.*consts.l10);
        dcf26 = exp(-tvm.*consts.l26);
        N10p = trapz_m(tvm,repmat(P10z,size(pathv)).*dcf10);
        N26p = trapz_m(tvm,repmat(P26z,size(pathv)).*dcf26);
        beep = N10p./P10z(1);
        alep = N26p./P26z(1);
        path_bee1 = (1/consts.l10)*(1-exp(-consts.l10*1E7)).*exp(-consts.l10.*burv);
        path_bee2 = beep(1).*exp(-consts.l10.*burv);
        path_bee3 = beep(2).*exp(-consts.l10.*burv);
        path_bee4 = beep(3).*exp(-consts.l10.*burv);
        path_bee5 = beep(4).*exp(-consts.l10.*burv);
        path_ale1 = (1/consts.l26)*(1-exp(-consts.l26*1E7)).*exp(-consts.l26.*burv);
        path_ale2 = alep(1).*exp(-consts.l26.*burv);
        path_ale3 = alep(2).*exp(-consts.l26.*burv);
        path_ale4 = alep(3).*exp(-consts.l26.*burv);
        path_ale5 = alep(4).*exp(-consts.l26.*burv);
        
        % plot burial pathways
        semilogx(path_bee1,path_ale1./path_bee1,'linestyle','--','color','blue');
        semilogx(path_bee2,path_ale2./path_bee2,'linestyle','--','color','blue');
        semilogx(path_bee3,path_ale3./path_bee3,'linestyle','--','color','blue');
        semilogx(path_bee4,path_ale4./path_bee4,'linestyle','--','color','blue');
        semilogx(path_bee5,path_ale5./path_bee5,'linestyle','--','color','blue');
    end;
end;

% probability contours
if pl.probcontours > 0;
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
if pl.sigma1line+pl.sigma2line > 0;
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
        
        if pl.sigma1line > 0;
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
        
        if pl.sigma2line > 0;
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
if pl.points > 0;
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


% subfunction trapz_m ==============================================================================
function z = trapz_m(x,y);
    dim = 1;
    nd = ndims(y);
    sz = size(y);
    n = sz(dim);
    idx1 = repmat ({':'}, [nd, 1]);
    idx2 = idx1;
    idx1{dim} = 2 : n;
    idx2{dim} = 1 : (n - 1);
    z = 0.5 * sum (diff (x, 1, dim) .* (y(idx1{:}) + y(idx2{:})), dim);
% end subfunction trapz_m ==========================================================================
