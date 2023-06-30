function burial();

% expage 10Be and 26Al exposure/burial duration calculator.
% Read and fix input, use get_exposure_burial to calculate exposure and burial durations and/or
% get_erosion_burial to calculate erosion rate and burial duration, and fix and write output. If
% plotting is on, display 26/10 banana plot.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2023 (jakob.heyman@gu.se)

clear; close all;
tic();

% What version is this?
ver = '202306';

% make choices here ================================================================================
% exposure + burial and/or erosion + burial? (1 = yes)
par.exposure_burial = 1; % calculate burial assuming prior surface exposure
par.erosion_burial = 0; % calculate burial assuming prior constant erosion

% calculate external and/or internal uncertainties? (1 = yes)
par.extunc = 1;
par.intunc = 1;
% Monte Carlo size for uncertainty estimate
mc = 1E5;

% plotting? (1 = yes)
pl.expline = 1;              % plot simple exposure line (and burial paths/lines)
pl.eline = 1;                % plot erosion end point line (and burial paths/lines)
pl.points = 1;               % plot sample points
pl.sigma1line = 1;           % plot sample 1 sigma line
pl.sigma2line = 0;           % plot sample 2 sigma line
pl.contourarea = 0;          % plot probability contour area
pl.contourline = 0;          % plot probability contour line
pl.extunc = 1;               % plot external uncertainty (including prod rate unc) for sample conc
pl.intunc = 0;               % plot internal uncertainty (excluding prod rate unc) for sample conc
pl.uncarea = 1;              % plot uncertainties as semi-transparent polygons
pl.sampleclr = 'red';        % color for sample conc
pl.explineclr = 'black';     % color for simple exposure line (and burial paths/lines)
pl.elineclr = [0.7 0.7 0.7]; % color for erosion end point line (and burial paths/lines)
pl.pointmarker = '.';        % marker for conc points
%~ pl.pointmarker = 'o';        % marker for conc points
pl.normalize = 1;            % normalize 26/10-P on Y-axis?

% write normalized conc data? (for later plotting)
par.normNout = 1;
% ==================================================================================================

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','lat','long','elv','Pflag','thick','dens','shield','erosion','N10','N10unc',...
    'std10','N26','N26unc','std26','samplingyr','pressure'};
vartypes = {'%s','%n','%n','%n','%s','%n','%n','%n','%n','%n','%n','%s','%n','%n','%s','%n','%n'};
% read input file
fid = fopen('input.txt');
varsin = strsplit(fgetl(fid)); % read first line
if(ismember(varsin,varnames)); % if first line contain only variable names
    [testi,vari] = ismember(varsin,varnames); % find index of varnames
    typestr = vartypes{vari(1)}; % fix type string
    for i = 2:numel(vari); % fix type string
        typestr = [typestr ' ' vartypes{vari(i)}];
    end;
elseif numel(varsin) == 15; % if no variable names in first line
    frewind(fid); % read from first line
    varsin = {'sample','lat','long','elv','Pflag','thick','dens','shield','N10','N10unc','std10',...
        'N26','N26unc','std26','samplingyr'};
    typestr = '%s %n %n %n %s %n %n %n %n %n %s %n %n %s %n';
else;
    fprintf(1,'ERROR! Something is wrong with the input\n');
    return;
end;
indata = textscan(fid,typestr,'CommentStyle','%','TreatAsEmpty','-'); % scan data
for i = 1:numel(varsin); % fix variables
    samplein.(varsin{i}) = indata{i};
end;
fclose(fid);
% ==================================================================================================

% run and load expage constants
make_consts_expage;
load consts_expage;

% Prefs
Pref10 = consts.Pref10; Pref10unc = consts.Pref10unc;
Pref26 = consts.Pref26; Pref26unc = consts.Pref26unc;
% decay constant
l10 = consts.l10; dell10 = consts.l10unc;
l26 = consts.l26; dell26 = consts.l26unc;

% display external/internal uncertainty info
if par.extunc == 1 && par.intunc == 1;
    fprintf(1,'external and (internal) uncertainties\n');
elseif par.extunc == 1;
    fprintf(1,'external uncertainties\n');
elseif par.intunc == 1;
    fprintf(1,'internal uncertainties\n');
end;

% if there is no N10 or N26 in input: fill with 0
if isfield(samplein,'N10') == 0; samplein.N10(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N10unc') == 0; samplein.N10unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std10') == 0; samplein.std10(1:numel(samplein.sample),1) = {'0'}; end;
if isfield(samplein,'N26') == 0; samplein.N26(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N26unc') == 0; samplein.N26unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std26') == 0; samplein.std26(1:numel(samplein.sample),1) = {'0'}; end;

% if there is NaN in N10 and N26: replace with 0
samplein.N10(isnan(samplein.N10)) = 0;
samplein.N10unc(isnan(samplein.N10unc)) = 0;
samplein.N26(isnan(samplein.N26)) = 0;
samplein.N26unc(isnan(samplein.N26unc)) = 0;

% convert 10Be concentrations according to standards
[testi,stdi] = ismember(samplein.std10,consts.std10); % find index of standard conversion factors
mult10 = consts.std10_cf(stdi); % pick out conversion factor
samplein.N10 = samplein.N10 .* mult10;
samplein.N10unc = samplein.N10unc .* mult10;

% convert 26Al concentrations according to standards
[testi,stdi] = ismember(samplein.std26,consts.std26); % find index of standard conversion factors
mult26 = consts.std26_cf(stdi); % pick out conversion factor
samplein.N26 = samplein.N26 .* mult26;
samplein.N26unc = samplein.N26unc .* mult26;

% fix longitude values
samplein.long(find(samplein.long < 0)) = samplein.long(find(samplein.long < 0)) + 360;

% fix sample pressure
if isfield(samplein,'pressure') == 0;
    % if there is no pressure flag: use std
    if isfield(samplein,'Pflag') == 0; samplein.Pflag(1:numel(samplein.sample),1) = {'std'}; end;
    stdv = strcmp(samplein.Pflag,'std');
    antv = strcmp(samplein.Pflag,'ant');
    prev = strcmp(samplein.Pflag,'pre');
    samplein.pressure(stdv) = ERA40atm(samplein.lat(stdv),samplein.long(stdv),samplein.elv(stdv));
    samplein.pressure(antv) = antatm(samplein.elv(antv));
    samplein.pressure(prev) = samplein.elv(prev);
end;

% declare pl parameters
pl.N10n = []; pl.N10uncintn = []; pl.N10uncextn = [];
pl.N26n = []; pl.N26uncintn = []; pl.N26uncextn = [];
pl.P10sp = []; pl.P26sp = []; pl.P_mu10 = []; pl.P_mu26 = []; pl.P10z = []; pl.P26z = [];
pl.Lsp = []; pl.RcEst = []; pl.SPhiAv = []; pl.thick_dens = []; pl.atm = []; pl.shield = [];

% fix for output =========================================================================
output(1,1) = {'sample'};
outcol.sample = 1;
outputstr = {'exposure(yr)','exp-uncext(yr)','exp-uncint(yr)';...
    'burial(yr)','bur-uncext(yr)','bur-uncint(yr)';...
    'erosion(mm/ka)','ero-uncext(mm/ka)','ero-uncint(mm/ka)';...
    'burial(yr)','bur-uncext(yr)','bur-uncint(yr)';...
    'N10norm','N10uncext','N10uncint';...
    'N26norm','N26uncext','N26uncint'};
outcolstr = {'exp','expextu','expintu';...
    'expbur','expburextu','expburintu';...
    'ero','eroextu','erointu';...
    'erobur','eroburextu','eroburintu';...
    'N10norm','N10uncext','N10uncint';...
    'N26norm','N26uncext','N26uncint'};
chstr = {'exposure_burial','erosion_burial','normNout'};
if sum(samplein.N10+samplein.N26) > 0;
    for i = 1:3;
        if par.(chstr{i}) == 1;
            output(1,end+1) = outputstr(i*2-1,1);
            outcol.(outcolstr{i*2-1,1}) = max(cell2mat(struct2cell(outcol)))+1;
            if par.extunc == 1;
                output(1,end+1) = outputstr(i*2-1,2);
                outcol.(outcolstr{i*2-1,2}) = max(cell2mat(struct2cell(outcol)))+1;
            end;
            if par.intunc == 1;
                output(1,end+1) = outputstr(i*2-1,3);
                outcol.(outcolstr{i*2-1,3}) = max(cell2mat(struct2cell(outcol)))+1;
            end;
            output(1,end+1) = outputstr(i*2,1);
            outcol.(outcolstr{i*2,1}) = max(cell2mat(struct2cell(outcol)))+1;
            if par.extunc == 1;
                output(1,end+1) = outputstr(i*2,2);
                outcol.(outcolstr{i*2,2}) = max(cell2mat(struct2cell(outcol)))+1;
            end;
            if par.intunc == 1;
                output(1,end+1) = outputstr(i*2,3);
                outcol.(outcolstr{i*2,3}) = max(cell2mat(struct2cell(outcol)))+1;
            end;
        end;
    end;
end;
% ========================================================================================

% pick out samples one by one
for i = 1:numel(samplein.lat);
    % pick out sample data
    samplefields = fieldnames(samplein);
    for j = 1:numel(samplefields);
        sample.(samplefields{j}) = samplein.(samplefields{j})(i);
        % remove fields with input '-' or numeric fields with NaN input
        if strcmp(sample.(samplefields{j}),'-') || ...
                (isnumeric(sample.(samplefields{j})) && isnan(sample.(samplefields{j})));
            sample = rmfield(sample,samplefields{j});
        end;
    end;
    
    % write sample name to output
    output(i+1,1) = sample.sample;
    
    % check that there is 10Be and 26Al data
    if sample.N10 .* sample.N26 <= 0;
        continue;
    end;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample{1});
    
    % Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
    LSDfix = LSD_fix(sample.lat,sample.long,1E7,-1,sample.samplingyr,consts);
    
    % interpolate Lsp (Sato, 2008; Marrero et al., 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,LSDfix.RcEst);
    
    % thickness scaling factor.
    if sample.thick > 0;
        thickSF = (sample.Lsp./(sample.dens.*sample.thick)).* ...
            (1-exp(((-1.*sample.dens.*sample.thick)./sample.Lsp)));
    else 
        thickSF = 1;
    end;
    
    % if calculating uncertainties
    if par.extunc+par.intunc >= 1;
        % randomized 10Be and 26Al concentrations for mc simulation
        sample.N10mc = normrnd(sample.N10,sample.N10unc,[1 mc]);
        sample.N26mc = normrnd(sample.N26,sample.N26unc,[1 mc]);
        % too low N not allowed!
        sample.N10mc(sample.N10mc<1) = 1;
        sample.N26mc(sample.N26mc<1) = 1;
        % randomized decay constants for mc simulation
        sample.l10mc = normrnd(l10,dell10,[1 mc]);
        sample.l26mc = normrnd(l26,dell26,[1 mc]);
    end;
    
    % spallation production scaling
    Psp = P_sp_expage(sample.pressure,LSDfix.RcEst,consts.SPhiInf,LSDfix.w,consts,1,1,0);
    
    % randomized ref prod rates for Monte Carlo simulation
    Pref10mc = normrnd(Pref10,Pref10unc,[1 mc]);
    Pref26mc = normrnd(Pref26,Pref26unc,[1 mc]);
    
    % spallation prod rates
    sample.Psp10 = Psp.sp10.*Pref10.*thickSF.*sample.shield;
    sample.Psp26 = Psp.sp26.*Pref26.*thickSF.*sample.shield;
    sample.Psp10mc = Psp.sp10.*Pref10mc.*thickSF.*sample.shield;
    sample.Psp26mc = Psp.sp26.*Pref26mc.*thickSF.*sample.shield;
    
    % fix shielding depth vector z (g/cm2/yr)
    if par.erosion_burial == 1;
        %~ sample.z = [0 logspace(0,5.3,100)]';
        sample.z = [(0:3:12) logspace(1.18,5.3,50)]';
    else;
        sample.z = 0;
    end;
    zmu = sample.z + sample.thick.*sample.dens./2; % add half sample depth
    
    % muon production
    P_mu = P_mu_expage(zmu,sample.pressure,LSDfix.RcEst,consts.SPhiInf,1,1,0,consts,'no');
    sample.Pmu10 = P_mu.mu10'.*sample.shield;
    sample.Pmu26 = P_mu.mu26'.*sample.shield;
    
    % normalized N (sent out in results for plotting)
    N10norm = sample.N10./(sample.Psp10+sample.Pmu10(1));
    N10norm_uncint = sample.N10unc./(sample.Psp10+sample.Pmu10(1));
    N10norm_uncext = sqrt(N10norm_uncint.^2 + (Pref10unc./Pref10.*N10norm).^2);
    N26norm = sample.N26./(sample.Psp26+sample.Pmu26(1));
    N26norm_uncint = sample.N26unc./(sample.Psp26+sample.Pmu26(1));
    N26norm_uncext = sqrt(N26norm_uncint.^2 + (Pref26unc./Pref26.*N26norm).^2);
    
    % calculate exposure + burial
    if par.exposure_burial == 1;
        expbur = get_exposure_burial(sample,consts,mc,par.extunc,par.intunc);
        
        % write to output
        output(i+1,outcol.exp) = {num2str(expbur.exposure,'%.0f')};
        output(i+1,outcol.expbur) = {num2str(expbur.burial,'%.0f')};
        if par.extunc == 1;
            output(i+1,outcol.expextu) = {num2str(expbur.exposure_uncext,'%.0f')};
            output(i+1,outcol.expburextu) = {num2str(expbur.burial_uncext,'%.0f')};
        end;
        if par.intunc == 1;
            output(i+1,outcol.expintu) = {num2str(expbur.exposure_uncint,'%.0f')};
            output(i+1,outcol.expburintu) = {num2str(expbur.burial_uncint,'%.0f')};
        end;
        
        % display exposure/burial and fix output if no solution
        if par.erosion_burial == 1; fprintf(1,'\n'); end; % separate lines if both exposure/erosion
        if expbur.exposure >= 0; % if initial exposure within bounds
            if par.extunc == 1 && par.intunc == 1;
                fprintf(1,' \texposure = %s ± %s (%s) yr',output{i+1,outcol.exp},...
                    output{i+1,outcol.expextu},output{i+1,outcol.expintu});
                fprintf(1,' \tburial = %s ± %s (%s) yr',output{i+1,outcol.expbur},output{i+1,outcol.expburextu},...
                    output{i+1,outcol.expburintu});
            elseif par.extunc == 1;
                fprintf(1,' \texposure = %s ± %s yr',output{i+1,outcol.exp},output{i+1,outcol.expextu});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outcol.expbur},output{i+1,outcol.expburextu});
            elseif par.intunc == 1;
                fprintf(1,' \texposure = %s ± %s yr',output{i+1,outcol.exp},output{i+1,outcol.expintu});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outcol.expbur},output{i+1,outcol.expburintu});
            else;
                fprintf(1,' \texposure = %s yr',output{i+1,outcol.exp});
                fprintf(1,' \tburial = %s yr',output{i+1,outcol.expbur});
            end;
        else;
            fprintf(1,' \tno solution!');
            output(i+1,outcol.exp) = {'-'}; output(i+1,outcol.expbur) = {'-'};
            if par.extunc == 1; output(i+1,outcol.expextu) = {'-'}; output(i+1,outcol.expburextu) = {'-'}; end;
            if par.intunc == 1; output(i+1,outcol.expintu) = {'-'}; output(i+1,outcol.expburintu) = {'-'}; end;
        end;
    end;
    
    % calculate erosion + burial
    if par.erosion_burial == 1;
        erobur = get_erosion_burial(sample,consts,mc,par.extunc,par.intunc);
        
        % write to output
        if erobur.erosion>0.1;
            output(i+1,outcol.ero) = {num2str(erobur.erosion,'%.2f')};
            if par.extunc == 1; output(i+1,outcol.eroextu) = {num2str(erobur.erosion_uncext,'%.2f')}; end;
            if par.intunc == 1; output(i+1,outcol.erointu) = {num2str(erobur.erosion_uncint,'%.2f')}; end;
        else;
            output(i+1,outcol.ero) = {num2str(erobur.erosion,'%.3f')};
            if par.extunc == 1; output(i+1,outcol.eroextu) = {num2str(erobur.erosion_uncext,'%.3f')}; end;
            if par.intunc == 1; output(i+1,outcol.erointu) = {num2str(erobur.erosion_uncint,'%.3f')}; end;
        end;
        output(i+1,outcol.erobur) = {num2str(erobur.burial,'%.0f')};
        if par.extunc == 1; output(i+1,outcol.eroburextu) = {num2str(erobur.burial_uncext,'%.0f')}; end;
        if par.intunc == 1; output(i+1,outcol.eroburintu) = {num2str(erobur.burial_uncint,'%.0f')}; end;
        
        % display erosion/burial and fix output if no solution
        if par.exposure_burial == 1; fprintf(1,'\n'); end; % separate lines if both exposure/erosion
        if erobur.erosion >= 0; % if initial erosion within bounds
            if par.extunc == 1 && par.intunc == 1;
                fprintf(1,' \terosion = %s ± %s (%s) mm/ka',output{i+1,outcol.ero},...
                    output{i+1,outcol.eroextu},output{i+1,outcol.erointu});
                fprintf(1,' \tburial = %s ± %s (%s) yr',output{i+1,outcol.erobur},...
                    output{i+1,outcol.eroburextu},output{i+1,outcol.eroburintu});
            elseif par.extunc == 1;
                fprintf(1,' \terosion = %s ± %s mm/ka',output{i+1,outcol.ero},output{i+1,outcol.eroextu});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outcol.erobur},output{i+1,outcol.eroburextu});
            elseif par.intunc == 1;
                fprintf(1,' \terosion = %s ± %s mm/ka',output{i+1,outcol.ero},output{i+1,outcol.erointu});
                fprintf(1,' \tburial = %s ± %s yr',output{i+1,outcol.erobur},output{i+1,outcol.eroburintu});
            else;
                fprintf(1,' \terosion = %s mm/ka',output{i+1,outcol.ero});
                fprintf(1,' \tburial = %s yr',output{i+1,outcol.erobur});
            end;
        else;
            fprintf(1,' \tno solution!');
            output(i+1,outcol.ero) = {'-'}; output(i+1,outcol.erobur) = {'-'};
            if par.extunc == 1; output(i+1,outcol.eroextu) = {'-'}; output(i+1,outcol.eroburextu) = {'-'}; end;
            if par.intunc == 1; output(i+1,outcol.erointu) = {'-'}; output(i+1,outcol.eroburintu) = {'-'}; end;
        end;
    end;
    
    % if writing normalized conc data
    if par.normNout > 0;
        output(i+1,outcol.N10norm) = {num2str(N10norm,'%.6f')};
        if par.extunc == 1; output(i+1,outcol.N10uncext) = {num2str(N10norm_uncext,'%.6f')}; end;
        if par.intunc == 1; output(i+1,outcol.N10uncint) = {num2str(N10norm_uncint,'%.6f')}; end;
        output(i+1,outcol.N26norm) = {num2str(N26norm,'%.6f')};
        if par.extunc == 1; output(i+1,outcol.N26uncext) = {num2str(N26norm_uncext,'%.6f')}; end;
        if par.intunc == 1; output(i+1,outcol.N26uncint) = {num2str(N26norm_uncint,'%.6f')}; end;
    end;
    
    % new line
    fprintf(1,'\n');
    
    % for plotting
    pl.N10n(end+1,1) = N10norm;
    pl.N10uncintn(end+1,1) = N10norm_uncint;
    pl.N10uncextn(end+1,1) = N10norm_uncext;
    pl.N26n(end+1,1) = N26norm;
    pl.N26uncintn(end+1,1) = N26norm_uncint;
    pl.N26uncextn(end+1,1) = N26norm_uncext;
    pl.P10sp(end+1,1) = sample.Psp10;
    pl.P26sp(end+1,1) = sample.Psp26;
    pl.P_mu10(end+1,1) = sample.Pmu10(1);
    pl.P_mu26(end+1,1) = sample.Pmu26(1);
    pl.Lsp(end+1,1) = sample.Lsp;
    pl.RcEst(end+1,1) = LSDfix.RcEst;
    pl.thick_dens(end+1,1) = sample.thick * sample.dens;
    pl.atm(end+1,1) = sample.pressure;
    pl.shield(end+1,1) = sample.shield;
    
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
if (pl.expline + pl.eline + pl.points + pl.sigma1line + pl.contourarea + pl.contourline) > 0;
    pl.exposure_burial = par.exposure_burial;
    pl.erosion_burial = par.erosion_burial;
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
        out.erosion_uncext = out.erosion_uncext./sample.dens.*1E4;
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
        out.erosion_uncint = out.erosion_uncint./sample.dens.*1E4;
    else;
        out.erosion_uncint = -1;
        out.burial_uncint = -1;
    end;
end;
% convert erosion to mm/ka
if out.erosion >= 0;
    out.erosion = out.erosion./sample.dens.*1E4;
end;
% end subfunction get_erosion_burial ===============================================================


% subfunction erobur_mc ============================================================================
function [erosion_unc burial_unc] = erobur_mc(z,N10mc,N26mc,l10mc,l26mc,P10mc,P26mc,minemc,...
    maxemc,mc,eN,burlim,ero_mid,bur_mid);
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
banana_plot = figure('name','Banana-plot','NumberTitle','off','visible','off');
%set(gcf,'visible','off'); % hide plots
hold on;
box on;

% fix for normalization
if pl.normalize == 1;
    pl.norm2610 = 1;
else; % use mean 26/10 ratio for all samples
    norm2610samples = (pl.P26sp+pl.P_mu26)./(pl.P10sp+pl.P_mu10);
    pl.norm2610 = mean(norm2610samples);
    pl.N26n = pl.N26n .* norm2610samples;
    pl.N26uncintn = pl.N26uncintn .* norm2610samples;
    pl.N26uncextn = pl.N26uncextn .* norm2610samples;
end;

% probability contours
if pl.contourarea == 1 || pl.contourline == 1;
    if pl.intunc == 1; pl_contours(pl,pl.N10uncintn,pl.N26uncintn); end;
    if pl.extunc == 1; pl_contours(pl,pl.N10uncextn,pl.N26uncextn); end;
end;

% simple exposure line,
if pl.expline > 0;
    pl_expline(pl,consts);
end;

% erosion line
if pl.eline == 1;
    pl_eline(pl,consts);
end;

% sigma lines
if pl.sigma1line == 1 || pl.sigma2line == 1;
    if pl.intunc == 1; pl_sigmalines(pl,pl.N10uncintn,pl.N26uncintn); end;
    if pl.extunc == 1; pl_sigmalines(pl,pl.N10uncextn,pl.N26uncextn); end;
end;

% sample points
if pl.points == 1 && strcmp(pl.pointmarker,'.');
    plot(pl.N10n,pl.N26n./pl.N10n,'.','color',pl.sampleclr,'markersize',15);
end;
if pl.points == 1 && strcmp(pl.pointmarker,'o');
    plot(pl.N10n,pl.N26n./pl.N10n,'o','color',pl.sampleclr,'markersize',3);
end;

% fix and display plot
% fix axis
if pl.normalize == 1;
    axis([1000 3000000 0.2 1.2]);
    ylabel('[^{26}Al]*/[^{10}Be]*');
else;
    axis([1000 3000000 1.5 8]);
    ylabel('[^{26}Al]/[^{10}Be]');
end;
xlabel('[^{10}Be]*');

set(gca,'layer','top'); % plot axis on top
set(gca,'XScale','log'); % fix for matlab
set(gcf,'visible','on'); % display plot
hold off;
figure(banana_plot);
% end subfunction plotting =========================================================================


% subfunction pl_expline ===========================================================================
function pl_expline(pl,consts);
    % create data for the simple-exposure line including ratio uncertainties
    tempt = [1 (100:100:900) logspace(3,7,60)];
    be = (1/consts.l10)*(1-exp(-consts.l10*tempt));
    al = (1/consts.l26)*(1-exp(-consts.l26*tempt));
    al_be = al./be;

    % color
    clr = pl.explineclr;
    
    % plot simple exposure line
    semilogx(be,al_be.*pl.norm2610,'color',clr);
    
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
        semilogx(bu_be50,bu_al50./bu_be50.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(bu_be100,bu_al100./bu_be100.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(bu_be150,bu_al150./bu_be150.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(bu_be200,bu_al200./bu_be200.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(bu_be250,bu_al250./bu_be250.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(bu_be300,bu_al300./bu_be300.*pl.norm2610,'linestyle','--','color',clr);
        
        % make data for burial pathways
        burv = (0:1E5:1E7);
        path_be10k = (1/consts.l10)*(1-exp(-consts.l10*1E4)).*exp(-consts.l10.*burv);
        path_be100k = (1/consts.l10)*(1-exp(-consts.l10*1E5)).*exp(-consts.l10.*burv);
        path_be1M = (1/consts.l10)*(1-exp(-consts.l10*1E6)).*exp(-consts.l10.*burv);
        path_be10M = (1/consts.l10)*(1-exp(-consts.l10*1E7)).*exp(-consts.l10.*burv);
        path_al10k = (1/consts.l26)*(1-exp(-consts.l26*1E4)).*exp(-consts.l26.*burv);
        path_al100k = (1/consts.l26)*(1-exp(-consts.l26*1E5)).*exp(-consts.l26.*burv);
        path_al1M = (1/consts.l26)*(1-exp(-consts.l26*1E6)).*exp(-consts.l26.*burv);
        path_al10M = (1/consts.l26)*(1-exp(-consts.l26*1E7)).*exp(-consts.l26.*burv);
        
        % plot burial pathways
        semilogx(path_be10k,path_al10k./path_be10k.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(path_be100k,path_al100k./path_be100k.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(path_be1M,path_al1M./path_be1M.*pl.norm2610,'linestyle','--','color',clr);
        semilogx(path_be10M,path_al10M./path_be10M.*pl.norm2610,'linestyle','--','color',clr);
    end;
% end subfunction pl_expline =======================================================================


% subfunction pl_eline =============================================================================
function pl_eline(pl,consts);
    fprintf(1,'calculating erosion end-point line...');
    % Precompute P_mu(z) to ~200,000 g/cm2 
    % start at the average sample mid-depth.
    z = [(0:3:12) logspace(1.18,5.3,50)]';
    z_mu = z+(mean(pl.thick_dens)./2);
    P_mu_z = P_mu_expage(z_mu,mean(pl.atm),mean(pl.RcEst),consts.SPhiInf,1,1,0,consts,'no');
    P_mu_z10 = P_mu_z.mu10'.*mean(pl.shield);
    P_mu_z26 = P_mu_z.mu26'.*mean(pl.shield);
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
    
    semilogx(bee,ale_bee.*pl.norm2610,'color',pl.elineclr);
    
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
        semilogx(bu_bee50,bu_ale50./bu_bee50.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(bu_bee100,bu_ale100./bu_bee100.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(bu_bee150,bu_ale150./bu_bee150.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(bu_bee200,bu_ale200./bu_bee200.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(bu_bee250,bu_ale250./bu_bee250.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(bu_bee300,bu_ale300./bu_bee300.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        
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
        semilogx(path_bee1,path_ale1./path_bee1.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(path_bee2,path_ale2./path_bee2.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(path_bee3,path_ale3./path_bee3.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(path_bee4,path_ale4./path_bee4.*pl.norm2610,'linestyle','--','color',pl.elineclr);
        semilogx(path_bee5,path_ale5./path_bee5.*pl.norm2610,'linestyle','--','color',pl.elineclr);
    end;
% end subfunction pl_eline =========================================================================


% subfunction pl_contours ==========================================================================
function pl_contours(pl,N10unc,N26unc);
    % Estimate range and create mesh
    Rv = pl.N26n./pl.N10n;
    delRv = sqrt((N26unc./pl.N26n).^2 + (N10unc./pl.N10n).^2);
    
    xmin = min(pl.N10n - 4.*N10unc);
    xmax = max(pl.N10n + 4.*N10unc);
    Rmin = min(Rv.*(1 - 4.*delRv));
    Rmax = max(Rv.*(1 + 4.*delRv));
    
    [xa,ya] = meshgrid(xmin:0.01 * mean(N10unc):xmax,Rmin:0.01.*mean(Rv).*mean(delRv):Rmax);
    
    Proba = zeros(size(xa));
    
    for j = 1:numel(pl.N10n(:,1));
        % calculate PDF
        Proba1 = exp(-0.5.*(((ya.*xa-pl.N26n(j))./N26unc(j)).^2+((xa-pl.N10n(j))./...
            N10unc(j)).^2));
        Proba = Proba + Proba1;
    end;
    minProba = min(min(Proba));
    maxProba = max(max(Proba));
    minlevel = minProba+(maxProba-minProba)/11;
    if pl.contourarea == 1;
        contourf(xa,ya,Proba,[linspace(minlevel,maxProba,100)],'linestyle','none');
        colormap jet;
    end;
    if pl.contourline == 1;
        contour(xa,ya,Proba,[linspace(minlevel,maxProba,10)]);
        colormap jet;
    end;
% end subfunction pl_contours ======================================================================


% subfunction pl_sigmalines ========================================================================
function pl_sigmalines(pl,N10unc,N26unc);
    % loop for each sample with al/be data
    for j = 1:numel(pl.N10n);
        % Estimate range and create mesh
        R = (pl.N26n(j)/pl.N10n(j));
        delR = sqrt((N26unc(j) / pl.N26n(j))^2 + (N10unc(j) / pl.N10n(j))^2);
        
        [x,y] = meshgrid((pl.N10n(j) - 4*N10unc(j)):(0.1*N10unc(j)):(pl.N10n(j) +...
            4*N10unc(j)),(R*(1 - 4*delR)):(0.1*R*delR):(R*(1 + 4*delR)));
        
        % calculate PDF
        Prob = x.*exp(-0.5.*((((y.*x) - pl.N26n(j))./N26unc(j)).^2 + ...
            ((x - pl.N10n(j))./N10unc(j)).^2));
        
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
            if pl.uncarea == 1;
                sd1 = patch(x1,y1,pl.sampleclr,'EdgeColor','none','FaceAlpha',...
                    max(min([0.8./numel(N10unc),0.4]),0.2));
            else;
                sd1 = plot(x1,y1,'color',pl.sampleclr);
            end;
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
            if pl.uncarea == 1;
                sd2 = patch(x2,y2,pl.sampleclr,'EdgeColor','none','FaceAlpha',...
                    max(min([0.8./numel(N10unc),0.4]),0.2));
            else;
                sd2 = plot(x2,y2,'color',pl.sampleclr);
            end;
        end;
    end;
% end subfunction pl_sigmalines ====================================================================


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
