function expage_submerge()

% expage 10Be/26Al/14C exposure age calculator.
% Read and fix input, calculate exposure ages, expected exposure ages given deglaciation age and
% subsequesnt sealevel history, and fix and write output.
% This code calculates exposure ages like expage.m plus expected exposure ages given deglaciation
% age and subsequesnt sealevel history. The output includes the standard simple exposure age plus
% the difference between that exposure age and the expected exposure age. A positive difference
% implies cosmogenic nuclide inheritance.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2023 (jakob.heyman@gu.se)

clear;
close all;
tic();

% What version is this?
ver = '202306';

% plotting? (1 = yes) ==============================================================================
plots.pointages = 0; % plots exposure ages as points
plots.probdens = 0;  % plots exposure ages as probability density curves (for single groups of ages)
% plots sealevels and exposure ages against elevation =====================================
plots.seaageelv = 1; % plot sealevels and ages against elevation
plots.submarea = 0;  % plot submerged time-elevation domain as area
plots.submline = 0;  % plot submerged time-elevation domain as line (polygon)
plots.deglac = 0;    % plot vertical lines at sample input delgaciation times
plots.ymax = 80;    % max elevation for sealevel plot (m a.s.l.) =========================
% plot elevation change over time - only valid if input has isostP ========================
plots.elv = 0;       % plot sample elevation against time (all samples in one plot)
plots.delv = 1;      % plot sample elevation change against time (all samples in one plot)
plots.elv1 = 0;      % plot single sample elevation against time (NOTE! one plot per sample)
plots.delv1 = 0;     % plot single sample elevation change against time (NOTE! one plot per sample)
plots.maxT = 14000;  % max time for elevation plots =======================================
plots.color10 = 'red'; % color for 10Be
plots.color26 = 'blue'; % color for 26Al
plots.color14 = 'green'; % color for 14C
% ==================================================================================================

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','isostsubm','isostP','lat','long','elv','thick',...
    'dens','shield','erosion','N10','N10unc','N26','N26unc','N14','N14unc','samplingyr','deglac',...
    'pressure','isostsubmmod','isostPmod'};
vartypes = {'%s','%s','%s','%s','%s','%s','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n','%n','%n','%n','%n','%n','%n'};
% read input file
fid = fopen('input.txt');
varsin = strsplit(fgetl(fid)); % read first line
if(ismember(varsin,varnames)); % if first line contain only variable names
    [testi,vari] = ismember(varsin,varnames); % find index of varnames
    typestr = vartypes{vari(1)}; % fix type string
    for i = 2:numel(vari); % fix type string
        typestr = [typestr ' ' vartypes{vari(i)}];
    end;
else; % if no or wrong variable names in first line
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

% if there is no erosion in input: assume zero erosion
if isfield(samplein,'erosion') == 0; samplein.erosion(1:numel(samplein.sample),1) = 0; end;

% if there is no N10/N26/N14 in input: fill with 0
if isfield(samplein,'N10') == 0; samplein.N10(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N10unc') == 0; samplein.N10unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std10') == 0; samplein.std10(1:numel(samplein.sample),1) = {'0'}; end;
if isfield(samplein,'N26') == 0; samplein.N26(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N26unc') == 0; samplein.N26unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std26') == 0; samplein.std26(1:numel(samplein.sample),1) = {'0'}; end;
if isfield(samplein,'N14') == 0; samplein.N14(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N14unc') == 0; samplein.N14unc(1:numel(samplein.sample),1) = 0; end;

% if there is NaN in N10/N26/N14: replace with 0
samplein.N10(isnan(samplein.N10)) = 0;
samplein.N10unc(isnan(samplein.N10unc)) = 0;
samplein.N26(isnan(samplein.N26)) = 0;
samplein.N26unc(isnan(samplein.N26unc)) = 0;
samplein.N14(isnan(samplein.N14)) = 0;
samplein.N14unc(isnan(samplein.N14unc)) = 0;

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

% fix erosion rate unit (mm/ka -> cm/yr)
samplein.erosion = samplein.erosion .* 1E-4;

% fix isostatic adjustment input
if isfield(samplein,'isostP');
    samplein.isostP = isost_data_load(samplein.isostP);
end;
if isfield(samplein,'isostsubm');
    samplein.isostsubm = isost_data_load(samplein.isostsubm);
end;

% initial Lsp for simple age calculation
Lsp1 = consts.Lsp;

% fix simple thickness scaling factor
thick_v = find(samplein.thick > 0);
samplein.thickSF1(1:numel(samplein.thick),1) = 1;
samplein.thickSF1(thick_v) = (Lsp1./(samplein.dens(thick_v).*samplein.thick(thick_v))).*...
    (1 - exp(((-1.*samplein.dens(thick_v).*samplein.thick(thick_v))./Lsp1)));

% sealevel matrices
sealevel.tv = [];
sealevel.shorel = [];
sealevel.nanm = [];

% fix for output
output(1,1) = {'sample'};
if sum(samplein.N10+samplein.N10unc) > 0;
    output(1,end+1:end+4) = {'10age(yr)','10uncext(yr)','10uncint(yr)','10diff(yr)'};
    outidx.n10 = (size(output,2)-3:size(output,2));
end;
if sum(samplein.N26+samplein.N26unc) > 0;
    output(1,end+1:end+4) = {'26age(yr)','26uncext(yr)','26uncint(yr)','26diff(yr)'};
    outidx.n26 = (size(output,2)-3:size(output,2));
end;
if sum(samplein.N14+samplein.N14unc) > 0;
    output(1,end+1:end+4) = {'14age(yr)','14uncext(yr)','14uncint(yr)','14diff(yr)'};
    outidx.n14 = (size(output,2)-3:size(output,2));
end;

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
    % remove potential isostatic adjustment if having atmospheric pressure as input
    if isfield(sample,'isostP') && strcmp(sample.Pflag,'pre');
        sample = rmfield(sample,'isostP');
    end;
    
    % write sample name to output
    output(i+1,1) = sample.sample;
    
    % Set nucl and mt to 0 for both 10/26 and check if there is N10/N26
    nucl10 = 0; nucl26 = 0; nucl14 = 0; mt10 = 0; mt26 = 0; mt14 = 0;
    if (sample.N10 + sample.N10unc) > 0; nucl10 = 1; end;
    if (sample.N26 + sample.N26unc) > 0; nucl26 = 1; end;
    if (sample.N14 + sample.N14unc) > 0; nucl14 = 1; end;
    
    % if no 10Be/26Al/14C: move on
    if nucl10 + nucl26 + nucl14 == 0;
        continue;
    end;
    
    % Prefs
    sample.Pref10 = consts.Pref10; sample.Pref10unc = consts.Pref10unc;
    sample.Pref26 = consts.Pref26; sample.Pref26unc = consts.Pref26unc;
    sample.Pref14 = consts.Pref14; sample.Pref14unc = consts.Pref14unc;
    % decay constants
    sample.l10 = consts.l10; sample.l10unc = consts.l10unc;
    sample.l26 = consts.l26; sample.l26unc = consts.l26unc;
    sample.l14 = consts.l14; sample.l14unc = consts.l14unc;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample{1});
    
    % Find P scaling factor according to Stone/Lal
    P_St_SF = stone2000(sample.lat,sample.pressure,1) * sample.thickSF1 * sample.shield;
    
    % if 10Be measured: calculate max time
    if nucl10 == 1;
        [tsimple10,mt10] = get_mt(sample,sample.Pref10,P_St_SF,sample.l10,Lsp1,sample.N10);
    end;
    % if 26Al measured: calculate max time
    if nucl26 == 1;
        [tsimple26,mt26] = get_mt(sample,sample.Pref26,P_St_SF,sample.l26,Lsp1,sample.N26);
    end;
    % if 14C measured: calculate max time
    if nucl14 == 1;
        P_St_SF = P_St_SF/0.8;
        [tsimple14,mt14] = get_mt(sample,sample.Pref14,P_St_SF,sample.l14,Lsp1,sample.N14);
    end;
    
    % pick largest of mt10 and mt26 as max time
    mt = max([mt10 mt26 mt14]);
    mt = mt + 11000; % skb-fix
    
    % Age Relative to t0=2010 - LSD tv from LSDfix
    % tv = [0:10:50 60:100:2960 3060:200:74860 75060:1000:799060 800060:2000:2000060 1E7];
    
    % Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
    LSDfix = LSD_fix(sample.lat,sample.long,mt,-1,sample.samplingyr,consts);
    
    % include tv in sample
    sample.tv = LSDfix.tv;

    % add points in time for tv, Rc, SPhi, delv for emergence/submergence events
    if isfield(sample,'isostsubm');
        % fix delv
        sample.delv = isost_elv(sample.isostsubm{1},sample);
        if numel(sample.delv) == 1; sample.delv = repmat(sample.delv,size(sample.tv)); end;
        if isfield(sample,'isostsubmmod');
            sample.delv = (sample.delv - sample.elv) .* sample.isostsubmmod + sample.elv;
        end;
        % pick out boundary points for water/air change
        bp1 = find(diff(sample.delv < 0) ~= 0); % boundary points 1
        bp2 = bp1 + 1; % boundary points 2
        % check if any elevation boundary point is exactly 0 and remove
        rmv = find(sample.delv(bp1) == 0); bp1(rmv) = []; bp2(rmv) = [];
        rmv = find(sample.delv(bp2) == 0); bp1(rmv) = []; bp2(rmv) = [];
        % calculate ratio of step between sample.delv(bp1) and sample.delv(bp2)
        stepratio = -sample.delv(bp1)./(sample.delv(bp2)-sample.delv(bp1));
        % interpolate tv, Rc, SPhi
        sample.tv = [sample.tv,sample.tv(bp1)+(sample.tv(bp2)-sample.tv(bp1)).*stepratio];
        LSDfix.Rc = [LSDfix.Rc,LSDfix.Rc(bp1)+(LSDfix.Rc(bp2)-LSDfix.Rc(bp1)).*stepratio];
        LSDfix.SPhi = [LSDfix.SPhi,LSDfix.SPhi(bp1)+(LSDfix.SPhi(bp2)-LSDfix.SPhi(bp1)).*stepratio];
        % add 0s to delv
        sample.delv = [sample.delv,repmat(0,1,numel(bp1))];
        % sort tv, Rc, SPhi, delv
        [sample.tv,sortidx] = sort(sample.tv);
        LSDfix.Rc = LSDfix.Rc(sortidx);
        LSDfix.SPhi = LSDfix.SPhi(sortidx);
        sample.delv = sample.delv(sortidx);
    end;

    % add a point to deglac for tv, Rc, SPhi, delv
    if ~any(sample.tv == sample.deglac);
        % pick out boundary points for deglac
        bp1 = find(diff(sample.tv < sample.deglac) ~= 0); % boundary point 1
        bp2 = bp1 + 1; % boundary point 2
        % calculate ratio of step between sample.tv(bp1) and sample.tv(bp2)
        stepratio = (sample.deglac-sample.tv(bp1))./(sample.tv(bp2)-sample.tv(bp1));
        % interpolate Rc, SPhi
        LSDfix.Rc = [LSDfix.Rc,LSDfix.Rc(bp1)+(LSDfix.Rc(bp2)-LSDfix.Rc(bp1)).*stepratio];
        LSDfix.SPhi = [LSDfix.SPhi,LSDfix.SPhi(bp1)+(LSDfix.SPhi(bp2)-LSDfix.SPhi(bp1)).*stepratio];
        % add deglac to tv
        sample.tv = [sample.tv,sample.deglac];
        % sort tv, Rc, SPhi, delv
        [sample.tv,sortidx] = sort(sample.tv);
        LSDfix.Rc = LSDfix.Rc(sortidx);
        LSDfix.SPhi = LSDfix.SPhi(sortidx);
        % interpolate delv
        if isfield(sample,'delv');
            sample.delv = [sample.delv,sample.delv(bp1)+(sample.delv(bp2)-...
                sample.delv(bp1)).*stepratio];
            sample.delv = sample.delv(sortidx);
        end;
    end;
    
    % Production from muons
    if sample.erosion <= 0;
        Pmu = P_mu_expage(sample.thick.*sample.dens./2,sample.pressure,LSDfix.RcEst,...
            consts.SPhiInf,nucl10,nucl26,nucl14,consts,'no');
        if nucl10 == 1; sample.mu10 = Pmu.mu10 .* sample.shield; end;
        if nucl26 == 1; sample.mu26 = Pmu.mu26 .* sample.shield; end;
        if nucl14 == 1; sample.mu14 = Pmu.mu14 .* sample.shield; end;
    else;
        tv_z = (sample.tv.*sample.erosion + sample.thick./2) .* sample.dens; % T-d vect (g/cm^2)
        if nucl10 == 1;
            sample.mu10 = get_PmuE(sample,tv_z,tsimple10,LSDfix.RcEst,consts,1,0,0);
        end;
        if nucl26 == 1;
            sample.mu26 = get_PmuE(sample,tv_z,tsimple26,LSDfix.RcEst,consts,0,1,0);
        end;
        if nucl14 == 1;
            sample.mu14 = get_PmuE(sample,tv_z,tsimple14,LSDfix.RcEst,consts,0,0,1);
        end;
    end;
    
    % fix atmospheric pressure if using isostatic adjustment
    if isfield(sample,'isostP');
        % calculate elevation development
        sample.elvv = isost_elv(sample.isostP{1},sample);
        if isfield(sample,'isostPmod');
            sample.elvv = (sample.elvv - sample.elv) .* sample.isostPmod + sample.elv;
        end;
        % calculate atmospheric pressure
        if strcmp(sample.Pflag,'std');
            sample.pressure = ERA40atm(sample.lat,sample.long,sample.elvv);
        elseif strcmp(sample.Pflag,'ant');
            sample.pressure = antatm(sample.elvv);
        end;
        % change Pref to isostatic calibration values
        sample.Pref10 = consts.Pref10iso; sample.Pref10unc = consts.Pref10isounc;
        sample.Pref26 = consts.Pref26iso; sample.Pref26unc = consts.Pref26isounc;
        sample.Pref14 = consts.Pref14iso; sample.Pref14unc = consts.Pref14isounc;
        % fill plot matrices
        plots.tv(1:size(sample.tv,2),i) = sample.tv(:);
        plots.elvm(1:size(sample.elvv,2),i) = sample.elvv;
        plots.delvm(1:size(sample.elvv,2),i) = sample.elvv - sample.elv;
        plots.one(1:size(sample.tv,2),i) = 1;
    end;
    
    % spallation production scaling
    sample.Psp = P_sp_expage(sample.pressure,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,nucl10,nucl26,...
        nucl14);
    
    % interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,LSDfix.Rc);
    
    % Thickness scaling factor.
    if sample.thick > 0;
        sample.thickSF = (sample.Lsp./(sample.dens.*sample.thick)).*...
            (1 - exp(((-1.*sample.dens.*sample.thick)./sample.Lsp)));
    else;
        sample.thickSF = 1;
    end;
    
    % spallation depth dependence
    sample.dpfs = exp(-sample.tv.*sample.erosion.*sample.dens./sample.Lsp);
    
    % submergence fix: calculate post-glacial uplift through water
    if isfield(sample,'isostsubm');
        % check for submergence only for postglacial time
        seaidx = find(sample.tv <= sample.deglac,1,'last');
        sample.delv2 = sample.delv(1:seaidx);
        underwaterv = (sample.delv2 < 0);
        overwaterv = (sample.delv2 >= 0);
        waterdv = -sample.delv2 .* 1E2 .* underwaterv; % (cm)
        
        % submergence fix: fix matrices for plotting sealevel
        seatv = [sample.tv(sample.tv<sample.deglac) sample.deglac];
        sealevel.tv(1:numel(seatv),end+1) = seatv';
        sealevel.shorel(1:numel(seatv),end+1) = interp1(sample.tv,-sample.delv,seatv)' + sample.elv;
        sealevel.nanm(1:numel(seatv),end+1) = 1;
        
        % submergence fix: fix erosion vector and then depth vector zv
        Ev = repmat(sample.erosion,1,numel(sample.tv));
        Ev(1:seaidx) = Ev(1:seaidx) .* overwaterv;
        zv = cumtrapz(sample.tv,Ev) .* sample.dens; % erosion rate zv (g/cm2)
        zv(1:seaidx) = zv(1:seaidx) + waterdv; % full depth vector (g/cm2)
        
        % submergence fix: production from muons - no atmospheric pressure isostatic adjustment!
        z_mu1 = [(0:3:12) logspace(1.18,5.3,50)];
        muidx = min(find(z_mu1 >= max(zv)+sample.thick.*sample.dens./2));
        z_mu = z_mu1(1:muidx);
        P_mu_z = P_mu_expage(z_mu+(sample.thick.*sample.dens./2),sample.pressure(1),LSDfix.RcEst,...
            consts.SPhiInf,nucl10,nucl26,nucl14,consts,'no');
        if nucl10 == 1;
            sample.mu10subm = interp1(z_mu,P_mu_z.mu10,zv,'pchip') .* sample.shield; %P_mu over time
        end;
        if nucl26 == 1;
            sample.mu26subm = interp1(z_mu,P_mu_z.mu26,zv,'pchip') .* sample.shield; %P_mu over time
        end;
        if nucl14 == 1;
            sample.mu14subm = interp1(z_mu,P_mu_z.mu14,zv,'pchip') .* sample.shield; %P_mu over time
        end;
        
        % include in sample
        sample.dpfs_subm = exp(-zv./sample.Lsp); % submergence fix: spallation depth dependence
    end;

    % fix and calculate exposure ages
    if nucl10 == 1;
        [output,plots] = fix_and_calculate_inh(output,plots,sample,'10',1E7,outidx,i,' \t10Be');
    end;
    if nucl26 == 1;
        [output,plots] = fix_and_calculate_inh(output,plots,sample,'26',6E6,outidx,i,' \t26Al');
    end;
    if nucl14 == 1;
        [output,plots] = fix_and_calculate_inh(output,plots,sample,'14',6E4,outidx,i,' \t14C');
    end;
    
    fprintf(1,'\n')
    clear sample;
end;

% fix and save output ==============================================================================
if sum(samplein.N10+samplein.N10unc+samplein.N26+samplein.N26unc+samplein.N14+samplein.N14unc) > 0;
    % fix output string
    outstr = '%s';
    for j = 1:size(output,2)-1;
        outstr = strcat(outstr,'\t%s');
    end;
    outstr = strcat(outstr,'\n');

    % fill empty cells with '-'
    nullidx = cellfun(@isempty,output);
    output(nullidx) = {'-'};
    
    % write out-glacialE.txt
    out = fopen('out-expage.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    fclose(out);
end;
% ==================================================================================================

% plotting ========================================
if plots.pointages == 1 || plots.probdens == 1 || plots.seaageelv == 1;
    plot_fixing(plots,samplein,sealevel);
end;
if (plots.elv == 1 || plots.delv == 1 || plots.elv1 == 1 || plots.delv1) && isfield(plots,'tv');
    plot_Telv(plots,samplein.sample);
end;
% =================================================

toc()
clear;
% end expage function ==============================================================================


% subfunction get_mt ===============================================================================
function [tsimple,mt] = get_mt(sample,Pref,P_St_SF,l,Lsp1,N);
    % Find P according to Stone/Lal - no muon production here!
    % Pref used as ref prod rate
    P_St = Pref * P_St_SF;
    A = l + sample.dens * sample.erosion / Lsp1;
    if (N < (P_St./A)); % if not saturated: do calculation
        tsimple = (-1/A)*log(1-(N * A / P_St));
        mt = tsimple .* 2;
        if mt < 12060;
            mt = 12060; % don't allow unreasonably short times
        elseif mt > 1e7;
            mt = 1e7;
        end;
    else; % if saturated: use full tv
        tsimple = 5e6;
        mt = 1e7;
    end;
% end subfunction get_mt ===========================================================================


% subfunction get_PmuE =============================================================================
function out = get_PmuE(sample,tv_z,tsimple,RcEst,consts,nucl10,nucl26,nucl14);
    agez = sample.erosion .* tsimple .* sample.dens; % shielding depth at t simple
    agezv = [0 25 50 75 125 250]; % shielding depth break values for number of mu_z points (5-10)
    idx = find(agez>agezv,1,'last');
    mu_z = [(linspace(0,agez,idx+3) + sample.thick.*sample.dens./2) max(tv_z)];
    Pmu_d = P_mu_expage(mu_z,sample.pressure,RcEst,consts.SPhiInf,nucl10,nucl26,nucl14,consts,'no');
    if nucl10 == 1; Pmud = Pmu_d.mu10;
    elseif nucl26 == 1; Pmud = Pmu_d.mu26;
    elseif nucl14 == 1; Pmud = Pmu_d.mu14; end;
    out = interp1(mu_z,Pmud,tv_z,'pchip') .* sample.shield; % P_mu
% end subfunction get_PmuE =========================================================================


% subfunction fix_and_calculate_inh ================================================================
function [output,plots] = fix_and_calculate_inh(output,plots,sample,nucl,maxage,outidx,i,dispstr);
    % sample surface spallation production rate over time
    sample.sp = sample.Psp.(['sp' nucl]).*sample.(['Pref' nucl]).*sample.thickSF.*sample.shield;
    
    % sample muon P
    sample.mu = sample.(['mu' nucl]);
    if isfield(sample,'isostsubm');
        sample.musubm = sample.(['mu' nucl 'subm']);
    end;
    
    % various parameters
    sample.N = sample.(['N' nucl]); sample.Nunc = sample.(['N' nucl 'unc']);
    sample.Pref = sample.(['Pref' nucl]); sample.Prefunc = sample.(['Pref' nucl 'unc']);
    sample.l = sample.(['l' nucl]);
    
    % get age
    results = get_1026_age_inh(sample,maxage);
    
    % fill output
    output(i+1,outidx.(['n' nucl])) = results.outstr;
    
    % fill plot matrix
    plots.(['i' nucl])(i,1) = i;
    plots.(['agem' nucl])(i,1:3) = results.num(1:3);
    plots.(['Pref' nucl])(i,1) = sample.(['Pref' nucl]);
    plots.(['Pref' nucl 'unc'])(i,1) = sample.(['Pref' nucl 'unc']);
    plots.(['elv' nucl])(i,1) = sample.elv;
    
    % display age
    fprintf(1,strcat(dispstr,results.outdisp{1}));
% end subfunction fix_and_calculate_inh ============================================================


% subfunction get_1026_age_inh =====================================================================
function results = get_1026_age_inh(sample,maxage);
    % Calculate N(t) including decay and erosion
    dcf = exp(-sample.tv.*sample.l); % decay factor;
    Nv = cumtrapz(sample.tv,(sample.sp.*dcf.*sample.dpfs + sample.mu.*dcf)); % pot N back in time
    
    % Look for saturation
    if sample.N <= max(Nv); % if not saturated
        age = min(interp1(Nv,sample.tv,sample.N),maxage); % get age by reverse-interpolation

        % uncertainty propagation based on CRONUS v. 2 ====================================
        if age > 0;
            % A with integrated average Lsp
            Lsp_avg = interp1(sample.tv,cumtrapz(sample.tv,sample.Lsp.*...
                exp(-sample.l.*sample.tv)),min(age,max(sample.tv)))/interp1(sample.tv,...
                cumtrapz(sample.tv,exp(-sample.l.*sample.tv)),min(age,max(sample.tv)));
            A = sample.l + sample.dens * sample.erosion ./Lsp_avg;
            % do most of computation
            FP = (sample.N.*A)./(1 - exp(-A.*age));
            delFP = (sample.Prefunc/sample.Pref) * FP;
            dtdN = 1./(FP - sample.N.*A);
            dtdP = -sample.N./(FP.*FP - sample.N.*A.*FP);
            % make respective delt's
            extunc = sqrt(dtdN.^2 * sample.Nunc.^2 + dtdP.^2 * delFP.^2);
            intunc = sqrt(dtdN.^2 * sample.Nunc.^2);
        else; % t = 0, estimate uncertainty based on conc + unc
            intunc = interp1(Nv,sample.tv,sample.N+sample.Nunc);
            extunc = intunc * (1 + sample.Prefunc/sample.Pref);
        end; % end uncertainty block ======================================================
    else;
        age = maxage;
        intunc = maxage;
        extunc = maxage;
        age_diff = maxage;
    end;

    if isfield(sample,'isostsubm');
        % calculate expected N from deglaciation age
        Nv_deglac = cumtrapz(sample.tv,(sample.sp.*dcf.*sample.dpfs_subm + sample.musubm.*dcf));
        N_deglac = interp1(sample.tv,Nv_deglac,sample.deglac);
        
        % calculate expected simple exposure age and mismatch with actual simple exposure age
        age_expected = interp1(Nv,sample.tv,N_deglac); % get age by reverse-interpolation
        age_diff = age - age_expected;
    else;
        age_diff = 0;
    end;
    
    % fix results
    results.num(1) = round(age);
    results.num(2) = round(extunc);
    results.num(3) = round(intunc);
    results.num(4) = round(age_diff);
    results.outstr(1) = {num2str(age,'%.0f')};
    results.outstr(2) = {num2str(extunc,'%.0f')};
    results.outstr(3) = {num2str(intunc,'%.0f')};
    results.outstr(4) = {num2str(age_diff,'%.0f')};
    results.outdisp = strcat({' = '},results.outstr{1},{' ± '},results.outstr{2},{' yr'});
    if sample.N+sample.Nunc > max(Nv); % notification for saturated samples
        results.outdisp = strcat(results.outdisp{1},{' (saturated!)'});
    end;
    if round(age_diff) > 0;
        results.outdisp = strcat(results.outdisp{1},{' (+'},results.outstr{4},{' yr)'});
    else;
        results.outdisp = strcat(results.outdisp{1},{' ('},results.outstr{4},{' yr)'});
    end;
% end subfunction get_1026_age_inh =================================================================


% subfunction plot_fixing ==========================================================================
function plot_fixing(plots,samplein,sealevel);
    % This subfunction fixes the data for plotting 10Be and/or 26Al ages as:
    % points with external uncertainties
    % probability density plots
    % nuclide vector and arrays
    nucla = {'10','26','14'};
    maxagev = [1E7 6E6 6E4];
    lega = {'^{10}Be','^{26}Al','^{14}C'};
    % fix for axis limits
    for j = 1:numel(nucla);
        if isfield(plots,['i' nucla{j}]) == 0; plots.(['i' nucla{j}]) = []; end;
    end;
    for j = 1:numel(nucla);
        % if no agem: move on
        if isfield(plots,['agem' nucla{j}]) == 0; continue; end;
        % find nuclide index and fix nn
        nuclidx = find(strcmp(nucla,nucla{j})==1);
        nn = nucla{nuclidx};
        % fix axis limits for point plot
        axislim = [0 max([plots.i10;plots.i26;plots.i14])+1 0 ...
            min(max(plots.(['agem' nn])(:,1)+plots.(['agem' nn])(:,2)),maxagev(nuclidx))];
        % remove samples without N
        rmidx = find(plots.(['agem' nn])(:,2) == 0);
        plots.(['i' nn])(rmidx) = [];
        plots.(['agem' nn])(rmidx,:) = [];
        plots.(['elv' nn])(rmidx,:) = [];
        plots.(['Pref' nn])(rmidx) = [];
        plots.(['Pref' nn 'unc'])(rmidx) = [];
        % fix legend and x/y-max
        nucl_legend = [lega{nuclidx} ' exposure age (yr)'];
        % plot points
        if plots.pointages == 1;
            plot_points(plots.(['i' nn]),plots.(['agem' nn]),nucl_legend,axislim);
        end;
        % plot probability density
        if plots.probdens == 1;
            plot_probdens(plots.(['agem' nn]),plots.(['Pref' nn]),...
                plots.(['Pref' nn 'unc']),nucl_legend,maxagev(nuclidx));
        end;
    end;
    % plot exposure ages, elevation, and sealevel history
    if plots.seaageelv == 1;
        plot_seaageelv(samplein,plots,sealevel);
    end;
% end subfunction plot_fixing ======================================================================


% subfunction plot_points ==========================================================================
function plot_points(ploti,agem,nucl_legend,axislim);
    % This subfunction plots exposure ages as points with external uncertainties
    figure('name','Exposure ages','NumberTitle','off');
    hold on; box on;
        
    % plot uncertainty lines
    ymax = [];
    plot([ploti';ploti'],[agem(:,1)'+agem(:,2)';agem(:,1)'-agem(:,2)'],'color','black');

    % plot points
    plot(ploti,agem(:,1),'.','markersize',15,'color','black');

    % fix plot
    axis(axislim);
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    xlabel('Samples');
    ylabel(nucl_legend);
    hold off;
% end subfunction plot_points ======================================================================


% subfunction plot_probdens ========================================================================
function plot_probdens(agem,Pref,Prefunc,nucl_legend,xmaxmax);
    % This subfunction plots probability density curves
    % red curves: single age probability density curves using internal uncertainty
    % black curve: summed probability density curve
    % black vertical line: weighted mean age
    % grey area: weighted uncertainty with propagated production rate uncertainty added
    figure('name','Exposure age probability','NumberTitle','off');
    hold on; box on;

    % find number of saturated samples and remove them
    rmidx = find(agem(:,2) == xmaxmax);
    agem(rmidx,:) = [];
    nsat = numel(rmidx);

    % vectors and matrices for probability estimation
    ages = agem(:,1);
    uncs = agem(:,3);
    agemin = min(ages - 4.*uncs);
    agemax = max(ages + 4.*uncs);
    timev = linspace(agemin,agemax,500);
    timem = repmat(timev,numel(ages),1);
    agesm = repmat(ages,1,numel(timev));
    uncsm = repmat(uncs,1,numel(timev));

    % estimate probability distr and sum
    probdensmatr = normpdf(timem,agesm,uncsm);
    probsum = sum(probdensmatr);

    if numel(ages)>1;
        % calculate weigthed mean age and uncertainty
        [wage,wunc] = w_mean_unc(ages,uncs);
        wuncp = sqrt(wunc^2 + (wage.*mean(Prefunc./Pref))^2); % add prodrate uncertainty
        
        % plot grey age uncertainty regions
        uncage = linspace(max(wage-wunc,agemin),min(wage+wunc,agemax),100);
        uncagep = linspace(max(wage-wuncp,agemin),min(wage+wuncp,agemax),100);
        uncprob = interp1(timev,probsum,uncage,'pchip');
        uncprobp = interp1(timev,probsum,uncagep,'pchip');
        uncage = [uncage,wage+wunc,wage-wunc];
        uncagep = [uncagep,wage+wuncp,wage-wuncp];
        uncprob = [uncprob 0 0];
        uncprobp = [uncprobp 0 0];
        patch(uncagep,uncprobp,[0.9 0.9 0.9],'EdgeColor','none');
        patch(uncage,uncprob,[0.8 0.8 0.8],'EdgeColor','none');
        
        % plot Pref line
        wageprob = interp1(timev,probsum,wage,'pchip');
        plot([wage wage],[0 wageprob],'color','black');
    end;

    % plot individual sample prob dens curves
    plot(timev',probdensmatr','color','red');

    if numel(ages)>1;
        % plot summed prob dens curve
        plot(timev,sum(probdensmatr),'color','black');
    end;

    % fix plot
    xlim([max(agemin,0) min(agemax,xmaxmax)]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    xlabel(nucl_legend);
    ylabel('Relative probability');
    set(gca,'layer','top'); % plot axis on top
    if nsat > 0; % if saturated samples: display number of excluded samples
        if nsat == 1; satstr = 'sample'; else satstr = 'samples'; end;
        xlims = xlim;
        ylims = ylim;
        xmin = xlims(1) + 0.75*(xlims(2)-xlims(1));
        text(xmin,0.9*ylims(2),{[num2str(nsat,'%.0f'),' saturated'],[satstr,' excluded']});
    end;
    hold off;
% end subfunction plot_probdens ====================================================================


% subfunction plot_seaageelv =======================================================================
function plot_seaageelv(samplein,plots,sealevel);
    % This subfunction plots 10Be/26Al/14C exposure ages against sample elevation plus sample
    % sea-level curves back to the input deglaciation age.
    % Optional plotting of the submerged domain as line or area (patch) of maximum submergence.
    % Optional plotting of deglaciation lines

    % return if nothing to plot
    if numel(sealevel.tv) == 0; return; end;
    
    figure('name','Exposure age - Sealevel','NumberTitle','off'); hold on; box on;
    
    legh = []; legin = {};
    % fix and plot sealevel(s)
    sealevel.tv(sealevel.nanm~=1) = nan;
    sealevel.shorel(sealevel.nanm~=1) = nan;
    % fix and plot submerged domain
    tvsea = max(sealevel.tv');
    shoreline = max(sealevel.shorel');
    if plots.submarea == 1;
        shoreline(shoreline>plots.ymax) = plots.ymax;
        legh(end+1) = patch('XData',[tvsea tvsea(end) 0]./1000,'YData',[shoreline 0 0],...
            'FaceColor',[0.8 0.8 1],'EdgeColor','none');
        legin(end+1) = {'Submerged'};
    end;
    if plots.submline == 1;
        shoreline(shoreline>plots.ymax) = plots.ymax;
        legh(end+1) = plot([tvsea tvsea(end) 0]./1E3,[shoreline 0 0],'color',[0.8 0.8 1]);
        legin(end+1) = {'Submerged-line'};
    end
    % plot deglaciation line(s)
    if plots.deglac == 1;
        samplein.deglac(samplein.N10+samplein.N26+samplein.N14==0) = [];
        deglx = [samplein.deglac';samplein.deglac']./1E3;
        degly = [repmat(plots.ymax,1,numel(samplein.deglac));repmat(0,1,numel(samplein.deglac))];
        legtemp = plot(deglx,degly,'color',[0.5 0.5 0.5]);
        legh(end+1) = legtemp(1);
        legin(end+1) = {'Deglaciation'};
    end;
    % plot sealevel line(s)
    legtemp = plot(sealevel.tv./1E3,sealevel.shorel,'color',[0.5 0.5 1]);
    legh(end+1) = legtemp(1);
    legin(end+1) = {'Sealevel'};

    % check which nuclides exist and fix arrays
    nucla = {}; lega = {};
    if sum(samplein.N14) > 0;
        nucla(end+1) = '14';
        lega(end+1) = '^{14}C';
    end;
    if sum(samplein.N26) > 0;
        nucla(end+1) = '26';
        lega(end+1) = '^{26}Al';
    end;
    if sum(samplein.N26) > 0;
        nucla(end+1) = '10';
        lega(end+1) = '^{10}Be';
    end;
    % plot 14C/26Al/10Be samples
    for j = 1:numel(nucla);
        nn = nucla{j};
        % fix legend
        legh(2:end+1) = legh(1:end);
        legin(2:end+1) = legin(1:end);
        plot([plots.(['agem' nn])(:,1)+plots.(['agem' nn])(:,2),plots.(['agem' nn])(:,1)-...
            plots.(['agem' nn])(:,2)]'./1E3,[plots.(['elv' nn]),plots.(['elv' nn])]','color',...
            plots.(['color' nn]));
        legh(1) = plot(plots.(['agem' nn])(:,1)./1E3,plots.(['elv' nn]),'marker','o',...
            'markersize',4,'linestyle','none','color',plots.(['color' nn]));
        legin(1) = lega(j);
    end;
    % fix plot
    ylim([0 plots.ymax]);
    xlim([0 plots.maxT/1000]);
    set(gca(),'xdir','reverse');
    xlabel('Age (ka)');
    ylabel('Elevation (m a.s.l.)');
    legend(legh,legin,'location','northwest');
    set(gca,'layer','top'); % plot axis on top
    hold off;
% end subfunction plot_seaageelv ===================================================================


% subfunction plot_Telv ============================================================================
function plot_Telv(plots,sample);
    % This subfunction fixes the data for plotting time-elevation graphs:
    % fix for missing agem10/agem26/agem14
    if isfield(plots,'agem10') == 0; plots.agem10 = zeros(size(plots.tv,2),2); end;
    if isfield(plots,'agem26') == 0; plots.agem26 = zeros(size(plots.tv,2),2); end;
    if isfield(plots,'agem14') == 0; plots.agem14 = zeros(size(plots.tv,2),2); end;
    % fix agev10/agev26/agev14
    agev10 = plots.agem10(:,1)';
    agev26 = plots.agem26(:,1)';
    agev14 = plots.agem14(:,1)';
    % fix for agev10/agev26/agev14 with missing numbers
    if numel(agev10) < size(plots.tv,2); agev10(end+1:size(plots.tv,2)) = 0; end;
    if numel(agev26) < size(plots.tv,2); agev26(end+1:size(plots.tv,2)) = 0; end;
    if numel(agev14) < size(plots.tv,2); agev14(end+1:size(plots.tv,2)) = 0; end;
    % remove empty columns
    nullv = (sum(plots.one) == 0);
    plots.tv(:,nullv) = [];
    plots.one(:,nullv) = [];
    plots.elvm(:,nullv) = [];
    plots.delvm(:,nullv) = [];
    agev10(nullv) = [];
    agev26(nullv) = [];
    agev14(nullv) = [];
    sample(nullv) = [];
    % remove ending zeros (division by zero yields Inf)
    plots.tv = plots.tv ./ plots.one;
    plots.elvm = plots.elvm ./ plots.one;
    plots.delvm = plots.delvm ./ plots.one;
    
    % plot sample elevation against time (all samples in one plot)
    if plots.elv == 1;
        figure('name','T-elv','NumberTitle','off'); hold on; box on;
        plot_Telv_plot(plots.tv,plots.elvm,agev10,agev26,agev14,plots.maxT,'Elevation (m)',...
            plots.one);
    end;
    if plots.delv == 1;
        figure('name','T-delv','NumberTitle','off'); hold on; box on;
        plot_Telv_plot(plots.tv,plots.delvm,agev10,agev26,agev14,plots.maxT,...
            'Elevation change (m)',plots.one);
    end;
    if plots.elv1 == 1;
        for j = 1:size(plots.tv,2);
            figure('name',['T-elv ' sample{j}],'NumberTitle','off'); hold on; box on;
            title(sample{j});
            plot_Telv_plot(plots.tv(:,j),plots.elvm(:,j),agev10(j),agev26(j),agev14(j),...
                plots.maxT,'Elevation (m)',plots.one(:,j));
        end;
    end;
    if plots.delv1 == 1;
        for j = 1:size(plots.tv,2);
            figure('name',['T-delv ' sample{j}],'NumberTitle','off'); hold on; box on;
            title(sample{j});
            plot_Telv_plot(plots.tv(:,j),plots.delvm(:,j),agev10(j),agev26(j),agev14(j),...
                plots.maxT,'Elevation change (m)',plots.one(:,j));
        end;
    end;
% end subfunction plot_Telv ========================================================================


% subfunction plot_Telv_plot =======================================================================
function plot_Telv_plot(tvm,elvm,agev10,agev26,agev14,maxT,yleg,onem);
    % this subfunction plots sample elevation or elevation change against time
    % plot elevation curve
    plot(tvm,elvm,'color','black');
    % plot 10Be/26Al/14C ages on elevation curve
    leg = []; legin = {};
    if sum(agev10) > 0;
        nullm10 = (agev10==0);
        agev10(nullm10) = [];
        elvm10 = elvm;
        elvm10(:,nullm10) = [];
        if size(tvm,2) > 1;
            ageelv10 = diag(interp1(max(tvm')',elvm10,agev10))';
        else;
            ageelv10 = interp1(tvm(onem==1),elvm10(onem==1),agev10);
        end;
        leg(end+1) = plot(agev10,ageelv10,'o','markersize',2,'LineWidth',3,'color','red');
        legin(end+1) = {'^{10}Be age'};
    end;
    if sum(agev26) > 0;
        nullm26 = (agev26==0);
        agev26(nullm26) = [];
        elvm26 = elvm;
        elvm26(:,nullm26) = [];
        if size(tvm,2) > 1;
            ageelv26 = diag(interp1(max(tvm')',elvm26,agev26))';
        else;
            ageelv26 = interp1(tvm(onem==1),elvm26(onem==1),agev26);
        end;
        leg(end+1) = plot(agev26,ageelv26,'o','markersize',2,'LineWidth',3,'color','blue');
        legin(end+1) = {'^{26}Al age'};
    end;
    if sum(agev14) > 0;
        nullm14 = (agev14==0);
        agev14(nullm14) = [];
        elvm14 = elvm;
        elvm14(:,nullm14) = [];
        if size(tvm,2) > 1;
            ageelv14 = diag(interp1(max(tvm')',elvm14,agev14))';
        else;
            ageelv14 = interp1(tvm(onem==1),elvm14(onem==1),agev14);
        end;
        leg(end+1) = plot(agev14,ageelv14,'o','markersize',2,'LineWidth',3,'color','green');
        legin(end+1) = {'^{14}C age'};
    end;
    xlabel('Age (yr)');
    ylabel(yleg);
    xlim([0 maxT]);
    set(gca(),'xdir','reverse');
    legend(leg,legin,'location','northwest');
    set(gca,'layer','top'); % plot axis on top
    hold off;
% end subfunction plot_Telv_plot ===================================================================

% subfunction evm ==================================================================================
function [wage,wunc] = evm(agev,uncv);
    % this subfunction calculates weighted average and uncertainty using the expected value method
    % (Birch and Singh 2014)
    % fix data (make matrices)
    agem = repmat(agev,numel(agev),1);
    uncm = repmat(uncv,numel(uncv),1);
    xm = repmat(agev',1,numel(agev));
    
    % calculate probability for each sample age
    Mui = sum(sqrt(2./(pi.*(2.*uncm).^2)).*exp(-(xm-agem).^2./(2.*uncm.^2)))./numel(agev);
    
    % calculate summed probability for all samples
    Muj = sum(Mui);
    
    % calculate sample weights
    wi = Mui./Muj;
    
    % calculate weighted mean age
    wage = sum(wi.*agev);
    
    % uncertainty estimation
    uncint = sqrt(sum(wi.^2.*uncv.^2));
    uncext = sqrt(sum(wi.*(agev-wage).^2));
    wunc = max(uncint,uncext);
% end subfunction evm ==============================================================================
