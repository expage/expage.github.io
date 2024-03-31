function prodrate()

% Function for 10Be/26Al/14C reference production rate calibration.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2024 (jakob.heyman@gu.se)

clear;
close all;
tic();

% What version is this?
ver = '202403';

% do choices here ==================================================================================
% plotting?
plsites = 0; % plot Pref probability curves (1 = yes)
plfull.plot = 1; % plot Pref for samples + sites + total
plfull.xmax = 8; % max x value in plot (number of sites) - only used if >0
plfull.printsvg = 1; % print svg file?
plfull.printscale = '-S1500,400'; % scale printing - only used if printsvg = 1
plfull.printscale = '-S1000,500';

% cluster test?
cl.Pcluster = 1; % exclude outliers to try to achieve a well-clustered group Pref (1 = yes)
% parameters for cluster test / outlier rejection - only used when cl.Pcluster = 1
cl.chiprob = 0.05; % lower p-value limit for chi-square probability test
cl.mingroupn = 3; % minimum number of samples in well-clustered group
cl.maxoutratio = 1/3; % maximum outlier ratio
% minimum uncertainty ratio used in cluster analysis only (not weighted mean production rate
% calculation) to not punish samples with low uncertainties
cl.minuncratio = 0.05;
% ==================================================================================================

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','isostP','isostsubm','site','lat','long','elv',...
    'thick','dens','shield','erosion','N10','N10unc','N26','N26unc','N14','N14unc','samplingyr',...
    'pressure','calage','calageunc','isostPmod','isostsubmmod'};
vartypes = {'%s','%s','%s','%s','%s','%s','%s','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n','%n','%n','%n','%n','%n','%n','%n','%n'};
% read input file
fid = fopen('input.txt');
varsin = strsplit(fgetl(fid)); % read first line
if(ismember(varsin,varnames)); % if first line contain only variable names
    [testi,vari] = ismember(varsin,varnames); % find index of varnames
    typestr = vartypes{vari(1)}; % fix type string
    for i = 2:numel(vari); % fix type string
        typestr = [typestr ' ' vartypes{vari(i)}];
    end;
elseif numel(varsin) == 18; % if no variable names in first line
    frewind(fid); % read from first line
    varsin = {'sample','lat','long','elv','Pflag','thick','dens','shield','erosion','N10',...
        'N10unc','std10','N26','N26unc','std26','samplingyr','calage','calageunc'};
    typestr = '%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n %n %n';
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

% decay constant
vars.l10 = consts.l10; l10unc = consts.l10unc;
vars.l26 = consts.l26; l26unc = consts.l26unc;
vars.l14 = consts.l14; l14unc = consts.l14unc;

% if there is no site in input: assume all samples are from one site
if isfield(samplein,'site') == 0; samplein.site(1:numel(samplein.sample),1) = {'noname'}; end;

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

% full number of samples
fulln = numel(samplein.sample);

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

% fix row numbers
samplein.rown = (1:1:numel(samplein.sample))';

% initial Lsp for simple age calculation
vars.Lsp1 = consts.Lsp;

% fix simple thickness scaling factor
thick_v = find(samplein.thick > 0);
samplein.thickSF1(1:fulln,1) = 1;
samplein.thickSF1(thick_v) = (vars.Lsp1./(samplein.dens(thick_v).*samplein.thick(thick_v))).*...
    (1 - exp(((-1.*samplein.dens(thick_v).*samplein.thick(thick_v))./vars.Lsp1)));

% declare Pfull
Pfull10.Psite = []; Pfull10.Psiteunc = []; Pfull10.PsiteOK = [];
Pfull26.Psite = []; Pfull26.Psiteunc = []; Pfull26.PsiteOK = [];
Pfull14.Psite = []; Pfull14.Psiteunc = []; Pfull14.PsiteOK = [];
Pfull10.P = []; Pfull10.Punc = []; Pfull10.PidxOK = []; Pfull10.site = {};
Pfull26.P = []; Pfull26.Punc = []; Pfull26.PidxOK = []; Pfull26.site = {};
Pfull14.P = []; Pfull14.Punc = []; Pfull14.PidxOK = []; Pfull14.site = {};

% fix for output
output(1,1) = {'sample'};
col.sample = 1;
if sum(samplein.N10+samplein.N10unc) > 0;
    output(1,2:3) = {'P10(at/g/yr)','P10unc(at/g/yr)'};
    col.P10 = size(output,2)-1; col.P10unc = size(output,2);
    if cl.Pcluster == 1;
        output(1,end+1) = {'P10-included?'};
        col.P10incl = size(output,2);
    end;
end;
if sum(samplein.N26+samplein.N26unc) > 0;
    output(1,end+1:end+2) = {'P26(at/g/yr)','P26unc(at/g/yr)'};
    col.P26 = size(output,2)-1; col.P26unc = size(output,2);
    if cl.Pcluster == 1;
        output(1,end+1) = {'P26-included?'};
        col.P26incl = size(output,2);
    end;
end;
if sum(samplein.N14+samplein.N14unc) > 0;
    output(1,end+1:end+2) = {'P14(at/g/yr)','P14unc(at/g/yr)'};
    col.P14 = size(output,2)-1; col.P14unc = size(output,2);
    if cl.Pcluster == 1;
        output(1,end+1) = {'P14-included?'};
        col.P14incl = size(output,2);
    end;
end;

% fix for site output
outputsite(1,1) = {'site'};
scol.site = 1;
if sum(samplein.N10+samplein.N10unc) > 0;
    outputsite(1,2:5) = {'P10(at/g/yr)','P10unc(at/g/yr)','P10-χ2R','P10-p-value'};
    scol.P10 = size(outputsite,2)-3; scol.P10unc = size(outputsite,2)-2;
    scol.chi10 = size(outputsite,2)-1; scol.pvalue10 = size(outputsite,2);
    if cl.Pcluster == 1;
        outputsite(1,end+1) = {'P10-included?'};
        scol.P10incl = size(outputsite,2);
    end;
end;
if sum(samplein.N26+samplein.N26unc) > 0;
    outputsite(1,end+1:end+4) = {'P26(at/g/yr)','P26unc(at/g/yr)','P26-χ2R','P26-p-value'};
    scol.P26 = size(outputsite,2)-3; scol.P26unc = size(outputsite,2)-2;
    scol.chi26 = size(outputsite,2)-1; scol.pvalue26 = size(outputsite,2);
    if cl.Pcluster == 1;
        outputsite(1,end+1) = {'P26-included?'};
        scol.P26incl = size(outputsite,2);
    end;
end;
if sum(samplein.N14+samplein.N14unc) > 0;
    outputsite(1,end+1:end+4) = {'P14(at/g/yr)','P14unc(at/g/yr)','P14-χ2R','P14-p-value'};
    scol.P14 = size(outputsite,2)-3; scol.P14unc = size(outputsite,2)-2;
    scol.chi14 = size(outputsite,2)-1; scol.pvalue14 = size(outputsite,2);
    if cl.Pcluster == 1;
        outputsite(1,end+1) = {'P14-included?'};
        scol.P14incl = size(outputsite,2);
    end;
end;

% fix for total P output
outputfull = {'P-total','P(at/g/yr)','Punc(at/g/yr)'};

% loop for single sites
while numel(samplein.site) > 0;
    % display site name
    fprintf(1,'\nSite: %s\n',samplein.site{1});
    
    % pick out all samples from the same site
    siteidx = find(strcmp(samplein.site,samplein.site(1)));

    % pick out sample data
    samplefields = fieldnames(samplein);
    for j = 1:numel(samplefields);
        sitein.(samplefields{j}) = samplein.(samplefields{j})(siteidx);
    end;
    
    % calculate sample Prefs
    [output,P10,P26,P14] = get_sample_Prefs(sitein,vars,output,col,consts);

    % fix for site calculations
    if numel(P10.Prefv)+numel(P26.Prefv)+numel(P14.Prefv) > 0;
        siterow = size(outputsite,1)+1;
        outputsite(siterow,scol.site) = sitein.site(1);
    end;

    if numel(P10.Prefv) >= 1;
        [output,outputsite,Pfull10] = fix_and_calc_site_Pref(output,outputsite,siterow,...
            P10,Pfull10,cl,col,scol,plsites,sitein,'10');
    end;
    if numel(P26.Prefv) >= 1;
        [output,outputsite,Pfull26] = fix_and_calc_site_Pref(output,outputsite,siterow,...
            P26,Pfull26,cl,col,scol,plsites,sitein,'26');
    end;
    if numel(P14.Prefv) >= 1;
        [output,outputsite,Pfull14] = fix_and_calc_site_Pref(output,outputsite,siterow,...
            P14,Pfull14,cl,col,scol,plsites,sitein,'14');
    end;

    % remove calculated samples
    for j = 1:numel(samplefields);
        samplein.(samplefields{j})(siteidx) = [];
    end;
    
    % clear variables
    clear sitein; clear P10; clear P26; clear P14; clear siteP10; clear siteP26; clear siteP14;
end;

% calculate total average
if numel(Pfull10.Psite) > 1 || numel(Pfull26.Psite) > 1 || numel(Pfull14.Psite) > 1;
    fprintf(1,'\nTotal Pref\n');
end;
if numel(Pfull10.Psite) > 1;
    [outputfull,Pfull10] = get_total_P(outputfull,Pfull10,'10');
end;
if numel(Pfull26.Psite) > 1;
    [outputfull,Pfull26] = get_total_P(outputfull,Pfull26,'26');
end;
if numel(Pfull14.Psite) > 1;
    [outputfull,Pfull14] = get_total_P(outputfull,Pfull14,'14');
end;

% add line break
fprintf(1,'\n');

% display plots
plotv = findobj('type','figure');
for i = 1:numel(plotv); figure(plotv(i)); end;

% plot sample + site + global P
if numel(Pfull10.Psite) > 1 && plfull.plot == 1; plot_fullP(Pfull10,plfull,10); end;
if numel(Pfull26.Psite) > 1 && plfull.plot == 1; plot_fullP(Pfull26,plfull,26); end;
if numel(Pfull14.Psite) > 1 && plfull.plot == 1; plot_fullP(Pfull14,plfull,14); end;

% fix and save output ============================
if size(output,1) > 1;
    % fix output string
    outstr = '%s';
    for j = 1:size(output,2)-1;
        outstr = strcat(outstr,'\t%s');
    end;
    outstr = strcat(outstr,'\n');

    % fill empty cells with '-'
    nullidx = cellfun(@isempty,output);
    output(nullidx) = {'-'};
    
    % fix outputsite
    if size(outputsite,1) > 1;
        % fix output string
        outstr2 = '%s';
        for j = 1:size(outputsite,2)-1;
            outstr2 = strcat(outstr2,'\t%s');
        end;
        outstr2 = strcat(outstr2,'\n');
        
        % fill empty cells with '-'
        nullidx = cellfun(@isempty,outputsite);
        outputsite(nullidx) = {'-'};
    end;
    
    % write out-prodrate.txt
    out = fopen('out-prodrate.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    if size(outputsite,1) > 1;
        fprintf(out,'\n');
        for i = 1:size(outputsite,1);
            fprintf(out,outstr2,outputsite{i,:});
        end;
    end;
    if size(outputfull,1) > 1;
        fprintf(out,'\n');
        for i = 1:size(outputfull,1);
            fprintf(out,'%s\t%s\t%s\n',outputfull{i,:});
        end;
    end;
    fclose(out);
end;
% ================================================

toc()
clear;
% end prodrate function ============================================================================


% subfunction get_sample_Prefs =====================================================================
function [output,P10,P26,P14] = get_sample_Prefs(sitein,vars,output,col,consts);
    % number of 10Be/26Al/14C samples
    num10 = 0; num26 = 0; num14 = 0;
    
    % declare Pref and Prow matrices
    P10.Prefv = []; P10.Puncv = []; P10.Prow = [];
    P26.Prefv = []; P26.Puncv = []; P26.Prow = [];
    P14.Prefv = []; P14.Puncv = []; P14.Prow = [];
    
    for i = 1:numel(sitein.sample);
        % pick out sample data
        samplefields = fieldnames(sitein);
        for j = 1:numel(samplefields);
            sample.(samplefields{j}) = sitein.(samplefields{j})(i);
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
        output(sample.rown+1,1) = sample.sample;
        
        % Set nucl and mt to 0 for 10/26/14 and check if there is N10/N26/N14
        nucl10 = 0; nucl26 = 0; nucl14 = 0; mt10 = 0; mt26 = 0; mt14 = 0;
        if (sample.N10 + sample.N10unc) > 0; nucl10 = 1; num10 = num10 + 1; end;
        if (sample.N26 + sample.N26unc) > 0; nucl26 = 1; num26 = num26 + 1; end;
        if (sample.N14 + sample.N14unc) > 0; nucl14 = 1; num14 = num14 + 1; end;
        
        if nucl10 + nucl26 + nucl14 == 0;
            continue;
        end;
        
        % Pref for simple age estimates
        Pref10 = consts.Pref10;
        Pref26 = consts.Pref26;
        Pref14 = consts.Pref14;
        
        % display sample name
        fprintf(1,'%.0f. %s',i,sample.sample{1});
        
        % Find P scaling factor according to Stone/Lal
        P_St_SF = stone2000(sample.lat,sample.pressure,1) * sample.thickSF1 * sample.shield;
        
        % if 10Be measured: calculate tsimple for case with erosion
        if nucl10 == 1;
            [tsimple10,mt10] = get_mt(sample,Pref10,P_St_SF,vars.l10,vars.Lsp1,sample.N10);
        end;
        % if 26Al measured: calculate tsimple for case with erosion
        if nucl26 == 1;
            [tsimple26,mt26] = get_mt(sample,Pref26,P_St_SF,vars.l26,vars.Lsp1,sample.N26);
        end;
        % if 14C measured: calculate tsimple for case with erosion
        if nucl14 == 1;
            [tsimple14,mt14] = get_mt(sample,Pref14,P_St_SF,vars.l14,vars.Lsp1,sample.N14);
        end;
        
        % use calibration age for scaling
        mt = sample.calage;
        
        % Age Relative to t0=2010 - LSD tv from LSDfix
        % tv = [0:10:50 60:100:2960 3060:200:74860 75060:1000:799060 800060:2000:2000060 1E7];
        
        % Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
        LSDfix = LSD_fix(sample.lat,sample.long,mt,-1,sample.samplingyr,consts);
        
        % include tv in sample
        sample.tv = LSDfix.tv;
        
        % fix submergence vectors if using isostatic adjustment
        if isfield(sample,'isostsubm');
            sample.elvv = isost_elv(sample.isostsubm{1},sample);
            if isfield(sample,'isostsubmmod');
                sample.elvv = (sample.elvv - sample.elv) .* sample.isostsubmmod + sample.elv;
            end;
            sample.overwater = (sample.elvv >= 0);
            sample.underwater = (sample.elvv < 0);
            sample.waterdepth = -sample.elvv .* sample.underwater .* 1E2; % cm (assumed density: 1)
            if sum(sample.underwater) == 0;
                sample = rmfield(sample,'isostsubm');
            end;
        end;
        
        % Production from muons
        if sample.erosion <= 0 && isfield(sample,'isostsubm') == 0;
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
            if isfield(sample,'isostPmod') && sample.isostPmod>0;
                sample.elvv = (sample.elvv - sample.elv) .* sample.isostPmod + sample.elv;
            end;
            % calculate atmospheric pressure
            if strcmp(sample.Pflag,'std');
                sample.pressure = ERA40atm(sample.lat,sample.long,sample.elvv);
            elseif strcmp(sample.Pflag,'ant');
                sample.pressure = antatm(sample.elvv);
            end;
            % estimate muogenic production
            if numel(sample.pressure) > 1 && min(sample.pressure)<max(sample.pressure);
                sample = get_Pmu_isost(sample,LSDfix.RcEst,consts,nucl10,nucl26,nucl14);
            end;
        end;
        
        % spallation production scaling
        sample.Psp = P_sp_expage(sample.pressure,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,nucl10,...
            nucl26,nucl14);
        
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
        
        if isfield(sample,'isostsubm');
            sample.zdepth = cumtrapz(sample.tv,sample.erosion.*sample.overwater).*sample.dens;
            sample.zdepth = sample.zdepth + sample.waterdepth;
            sample.dpfs = exp(-sample.zdepth./sample.Lsp);
        end;

        % fix and calculate Pref
        if nucl10 == 1;
            [output,P10] = fix_and_calc_Pref(output,sample,vars,col,P10,num10,'10');
        end;
        if nucl26 == 1;
            [output,P26] = fix_and_calc_Pref(output,sample,vars,col,P26,num26,'26');
        end;
        if nucl14 == 1;
            [output,P14] = fix_and_calc_Pref(output,sample,vars,col,P14,num14,'14');
        end;
        
        fprintf(1,'\n');
        clear sample;
    end;
% end subfunction get_sample_Prefs =================================================================


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
    if sample.erosion > 0;
        agez = sample.erosion .* tsimple .* sample.dens; % shielding depth at t simple
        agezv = [0 25 50 75 125 250]; % shielding depth break values for mu_z points (5-10)
        idx = find(agez>agezv,1,'last');
        mu_z = [(linspace(0,agez,idx+3) + sample.thick.*sample.dens./2) max(tv_z)];
    end;
    if isfield(sample,'isostsubm');
        mu_z = [sample.thick.*sample.dens./2];
        if sample.erosion > 0;
            idx0m = find(diff(sample.underwater)==1,1,'first');
            tsimple = sample.tv(idx0m);
            agez = sample.erosion .* tsimple .* sample.dens; % shielding depth at t simple
            agezv = [0 25 50 75 125 250]; % shielding depth break values for mu_z points (5-10)
            idx = find(agez>agezv,1,'last');
            mu_z = [(linspace(0,agez,idx+3) + sample.thick.*sample.dens./2)];
            tv_z = (cumtrapz(sample.tv,sample.erosion.*sample.overwater) + sample.thick./2) .* ...
                sample.dens;
        end;
        mu_z(end+1:end+5) = logspace(log10(10+mu_z(end)),log10(max(sample.waterdepth+1)+...
            mu_z(end)),5);
        tv_z = tv_z + sample.waterdepth;
    end;
    Pmu_d = P_mu_expage(mu_z,sample.pressure,RcEst,consts.SPhiInf,nucl10,nucl26,nucl14,consts,'no');
    if nucl10 == 1; Pmud = Pmu_d.mu10;
    elseif nucl26 == 1; Pmud = Pmu_d.mu26;
    elseif nucl14 == 1; Pmud = Pmu_d.mu14; end;
    out = interp1(mu_z,Pmud,tv_z,'pchip') .* sample.shield; % P_mu
    out(isnan(out)) = out(end); % fix for nan issue...
% end subfunction get_PmuE =========================================================================


% subfunction get_Pmu_isost ========================================================================
function sample = get_Pmu_isost(sample,RcEst,consts,nucl10,nucl26,nucl14);
    % calculate Pmu at 4 points from min to max atmospheric pressure
    atmv = linspace(min(sample.pressure),max(sample.pressure),4);
    for j = 1:numel(atmv);
        Pmu_isost = P_mu_expage(sample.thick.*sample.dens./2,atmv(j),RcEst,consts.SPhiInf,nucl10,...
            nucl26,nucl14,consts,'no');
        if nucl10 == 1; Pmu_isost10(j) = Pmu_isost.mu10; end;
        if nucl26 == 1; Pmu_isost26(j) = Pmu_isost.mu26; end;
        if nucl14 == 1; Pmu_isost14(j) = Pmu_isost.mu14; end;
    end;
    % interpret Pmu based on sample pressure and multiply
    if nucl10 == 1;
        Pmu_atm10 = interp1(atmv,Pmu_isost10,sample.pressure,'pchip') .* sample.shield;
        sample.mu10 = Pmu_atm10 ./ sample.mu10(1) .* sample.mu10;
    end;
    if nucl26 == 1;
        Pmu_atm26 = interp1(atmv,Pmu_isost26,sample.pressure,'pchip') .* sample.shield;
        sample.mu26 = Pmu_atm26 ./ sample.mu26(1) .* sample.mu26;
    end;
    if nucl14 == 1;
        Pmu_atm14 = interp1(atmv,Pmu_isost14,sample.pressure,'pchip') .* sample.shield;
        sample.mu14 = Pmu_atm14 ./ sample.mu14(1) .* sample.mu14;
    end;
% end subfunction get_Pmu_isost ====================================================================


% subfunction fix_and_calc_Pref ====================================================================
function [output,Pnucl] = fix_and_calc_Pref(output,sample,vars,col,Pnucl,numnucl,nucl);
    % various parameters
    sample.dcf = exp(-sample.tv.*vars.(['l' nucl]));
    sample.N = sample.(['N' nucl]); sample.Nunc = sample.(['N' nucl 'unc']);
    sample.l = vars.(['l' nucl]);
    sample.sp = sample.Psp.(['sp' nucl]);
    sample.mu = sample.(['mu' nucl]);
    
    % get sample Pref
    [Prefi,Prefiunc,Prefiuncint] = Pref_calc(sample,nucl);
    
    % fill output
    output(sample.rown+1,col.(['P' nucl]):col.(['P' nucl 'unc'])) = ...
        {num2str(Prefi,'%.3f'),num2str(Prefiunc,'%.3f')};
    
    % display Pref
    fprintf(1,[' \tP' nucl ' = %s ± %s at/g/yr'],output{sample.rown+1,col.(['P' nucl])},...
        output{sample.rown+1,col.(['P' nucl 'unc'])});
    
    % fill Pref vector for uncertainty estimation, plotting, and cluster analysis
    Pnucl.Prefv(numnucl) = Prefi;
    Pnucl.Puncv(numnucl) = Prefiuncint;
    Pnucl.Prow(numnucl) = sample.rown; % input row number
    Pnucl.calage(numnucl) = sample.calage;
    Pnucl.calageunc(numnucl) = sample.calageunc;
    Pnucl.site(numnucl) = sample.site;
% end subfunction fix_and_calc_Pref ================================================================


% subfunction Pref_calc ============================================================================
function [Pref,Prefunc,Prefuncint] = Pref_calc(sample,nuclstr);
    % sample spallation production rate scaling over time including decay and erosion
    Psp = sample.sp.*sample.dcf.*sample.dpfs.*sample.thickSF.*sample.shield;
    
    % sample muon P
    Pmu = sample.mu.*sample.dcf;
    
    % calculate N produced by muons and spallation
    Nmu = trapz(sample.tv,Pmu);
    Nsp = sample.N-Nmu;
    
    % calculate N without Pref produced by spallation
    Nsp_scaling = trapz(sample.tv,Psp);
    
    % calculate Pref
    Pref = Nsp./Nsp_scaling;

    % error propagation based on Balco et al. (2008) method
    % decay-weighted average Lsp used in thicksf and A
    Lsp_avg = trapz(sample.tv,sample.Lsp.*exp(-sample.l.*sample.tv))/...
        trapz(sample.tv,exp(-sample.l.*sample.tv));
    thicksf = (Lsp_avg./(sample.dens.*sample.thick)).*...
        (1 - exp(((-1.*sample.dens.*sample.thick)./Lsp_avg)));
    A = sample.l + sample.dens.*sample.erosion./Lsp_avg;
    SF = (Nsp.*A)./(Pref.*sample.shield.*thicksf.*(1-exp(-A.*sample.calage)));
    dPdN = A./(SF.*sample.shield.*thicksf.*(1-exp(-A.*sample.calage)));
    dPdt = Nsp.*A./(SF.*thicksf.*sample.shield).*A.*exp(-A.*sample.calage).*...
        ((1-exp(-A.*sample.calage)).^-2);
    % calculate internal Pref uncertainty
    Prefuncint = sample.Nunc.*dPdN;
    % add calibration age uncertainty to get external Pref uncertainty
    Prefunc = sqrt(Prefuncint.^2 + (sample.calageunc.*dPdt).^2);
% end subfunction Pref_calc ========================================================================


% subfunction fix_and_calc_site_Pref ===============================================================
function [output,outputsite,Pfulln] = fix_and_calc_site_Pref(output,outputsite,siterow,Pn,Pfulln,...
        cl,col,scol,plsites,sitein,nucl);
    % if too few samples for site calculation: fill Pfulln with sample Pref
    if (numel(Pn.Prefv)<2) || (cl.Pcluster==1 && numel(Pn.Prefv)<cl.mingroupn);
        Pfulln.P(Pn.Prow) = Pn.Prefv;
        Pfulln.Punc(Pn.Prow) = Pn.Puncv;
        Pfulln.site(Pn.Prow) = Pn.site;
        Pfulln.PsiteOK(end+1) = 0;
        Pfulln.Psite(end+1) = 0;
        Pfulln.Psiteunc(end+1) = 0;
        return;
    end;
    
    % calculate site Pref
    sitePn = get_site_Pref(Pn);
    
    % if doing cluster analysis
    if cl.Pcluster == 1;
        % find clustered Pref by removing outliers
        sitePn = get_cluster(Pn,sitePn,cl);
        
        % mark OK samples in output
        output(sitePn.OKrow+1,col.(['P' nucl 'incl'])) = {'X'};
    end;

    % fill outputsite
    outputsite(siterow,scol.(['P' nucl])) = {num2str(sitePn.Pref,'%.3f')};
    outputsite(siterow,scol.(['P' nucl 'unc'])) = {num2str(sitePn.Prefunc,'%.3f')};
    outputsite(siterow,scol.(['chi' nucl])) = {num2str(sitePn.rchisq,'%.3f')};
    outputsite(siterow,scol.(['pvalue' nucl])) = {num2str(sitePn.Pvalue,'%.3f')};
    % OK site Pref for output and full Pref calculation
    Pfulln.PsiteOK(end+1) = 1;
    if cl.Pcluster==1;
        if sitePn.Pvalue>=cl.chiprob && numel(sitePn.OKrow)>=cl.mingroupn;
            outputsite(siterow,scol.(['P' nucl 'incl'])) = {'X'};
        else;
            Pfulln.PsiteOK(end) = 0;
        end;
    end;
    % display site Pref
    fprintf(1,['Site Pref' nucl ' = %s ± %s at/g/yr'],outputsite{siterow,scol.(['P' nucl])},...
        outputsite{siterow,scol.(['P' nucl 'unc'])});
    if cl.Pcluster==1 && numel(sitePn.OKrow)<numel(Pn.Prow);
        fprintf(1,'  (%.0f outlier',numel(Pn.Prow)-numel(sitePn.OKrow));
        if numel(Pn.Prow)-numel(sitePn.OKrow) > 1; fprintf(1,'s)'); else; fprintf(1,')'); end;
    end;
    fprintf(1,'\nR-chi2 = %s    P-value = %s\n',outputsite{siterow,scol.(['chi' nucl])},...
        outputsite{siterow,scol.(['pvalue' nucl])});
    % fill Pfull for total average calculation and plotting
    Pfulln.Psite(end+1) = sitePn.Pref;
    Pfulln.Psiteunc(end+1) = sitePn.Prefunc;
    Pfulln.P(Pn.Prow) = Pn.Prefv;
    Pfulln.Punc(Pn.Prow) = Pn.Puncv;
    Pfulln.site(Pn.Prow) = Pn.site;
    if isfield(sitePn,'OKrow');
        Pfulln.PidxOK(end+1:end+numel(sitePn.OKrow)) = sitePn.OKrow;
    else;
        Pfulln.PidxOK(Pn.Prow) = Pn.Prow;
    end;
    % do plotting here...
    if plsites == 1;
        if isfield(sitePn,'include') == 0; sitePn.include = (1:1:numel(Pn.Prefv)); end;
        plot_Pref(Pn,sitePn,str2num(nucl),sitein.site{1});
    end;
% end subfunction fix_and_calc_site_Pref ===========================================================

% subfunction get_site_Pref ========================================================================
function siteP = get_site_Pref(Pn);
    % calculate weighted Pref, internal Pref uncertainty, reduced chi-square, and P-value
    [siteP.Pref,Puncint,siteP.rchisq,siteP.Pvalue] = w_mean_unc(Pn.Prefv,Pn.Puncv,1);
    
    % calculate calibration age uncertainty part
    calage_unc = mean(Pn.calageunc./Pn.calage);
    
    % calculate total prod rate uncertainty
    siteP.Prefunc = sqrt(Puncint^2 + (calage_unc.*siteP.Pref)^2);
% end subfunction get_site_Pref ====================================================================


% subfunction get_cluster ==========================================================================
function siteP = get_cluster(Pn,siteP,cl);
    % in variables
    sitePin = siteP;
    
    % define maximum number of outliers
    remove = floor(numel(Pn.Prefv)*cl.maxoutratio);
    r = 1; % outlier counter
    
    % well-clustered samples for plotting
    Pn.include = (1:1:numel(Pn.Prow));

    % calculate modified Pref uncerntainty based on minimum uncertainty ratio cl.minuncratio
    Puncvmod = Pn.Puncv;
    Puncvmod(Puncvmod<(siteP.Pref*cl.minuncratio)) = siteP.Pref*cl.minuncratio;
    % calculate modified rchi2 and Pvalue
    siteP.rchisq = 1/(numel(Pn.Prefv)-1) .* sum(((Pn.Prefv-siteP.Pref)./Puncvmod).^2);
    siteP.Pvalue = 1 - chi2cdf(siteP.rchisq.*(numel(Pn.Prefv)-1),numel(Pn.Prefv)-1);
    
    % loop for outlier removal
    while siteP.Pvalue<cl.chiprob && remove>=r && numel(Pn.Prefv)>=cl.mingroupn;
        % calculate deviation for all individual samples
        chidev = ((Pn.Prefv-siteP.Pref)./Puncvmod).^2;
        
        % find index of sample with largest dev
        [value,rmv_idx] = max(chidev);
        
        % remove outlier
        Pn.calage(rmv_idx) = []; Pn.calageunc(rmv_idx) = [];
        Pn.Prefv(rmv_idx) = []; Pn.Puncv(rmv_idx) = [];
        Pn.Prow(rmv_idx) = []; Pn.include(rmv_idx) = [];
        Puncvmod(rmv_idx) = [];
        
        % calculate Pref
        siteP = get_site_Pref(Pn);

        % calculate modified rchi2 and Pvalue
        siteP.rchisq = 1/(numel(Pn.Prefv)-1) .* sum(((Pn.Prefv-siteP.Pref)./Puncvmod).^2);
        siteP.Pvalue = 1 - chi2cdf(siteP.rchisq.*(numel(Pn.Prefv)-1),numel(Pn.Prefv)-1);
        
        % add outlier number
        r = r+1;
    end;
    
    % if not well-clustered: use input variables and remove row numbers
    if siteP.Pvalue < cl.chiprob || numel(Pn.Prow) < cl.mingroupn;
        % use input variables
        siteP = sitePin;
        Pn.Prow = []; Pn.include = [];
    end;
    
    % fix output
    siteP.OKrow = Pn.Prow;
    siteP.include = Pn.include;
% end subfunction get_cluster ======================================================================


% subfunction get_total_P ==========================================================================
function [outputfull,Pfull] = get_total_P(outputfull,Pfull,nucl);
    % calculate aritmetic mean of clustered site P
    Ptot = mean(Pfull.Psite(Pfull.PsiteOK==1));
    % calculate unbiased deviation from Ptot for all clustered single samples
    P = Pfull.P(Pfull.PidxOK);
    Punc = Pfull.Punc(Pfull.PidxOK);
    Ptotunc = sqrt(sum(1./Punc.^2.*(P-Ptot).^2) / ...
        (sum(1./Punc.^2)-(sum(1./Punc.^4)/sum(1./Punc.^2))));
    % display total P
    fprintf(1,['P' nucl ' = %.2f ± %.2f at/g/yr\n'],Ptot,Ptotunc);
    % fill outputfull and Pfull
    outputfull(size(outputfull,1)+1,1) = {['P' nucl]};
    outputfull(size(outputfull,1),2) = {num2str(Ptot,'%.2f')};
    outputfull(size(outputfull,1),3) = {num2str(Ptotunc,'%.2f')};
    Pfull.Ptot = Ptot;
    Pfull.Ptotunc = Ptotunc;
% end subfunction get_total_P ======================================================================


% subfunction plot_Pref ============================================================================
function plot_Pref(Pn,siteP,nucl,pltitle);
    % fix for 10Be/26Al/14C
    if nucl == 10;
        step = 0.01;
        nuclstr = '10Be';
        lablestr = '^{10}Be';
        clrin = 'red';
    elseif nucl == 26;
        step = 0.02;
        nuclstr = '26Al';
        lablestr = '^{26}Al';
        clrin = 'blue';
    elseif nucl == 14;
        step = 0.02;
        nuclstr = '14C';
        lablestr = '^{14}C';
        clrin = 'green';
    end;
    
    % split Pref vectors in clustered and outliers
    cPrefv = Pn.Prefv(siteP.include)';
    cPuncv = Pn.Puncv(siteP.include)';
    oPrefv = Pn.Prefv'; oPrefv(siteP.include) = [];
    oPuncv = Pn.Puncv'; oPuncv(siteP.include) = [];
    
    % define min and max Pref for plot
    plmin = floor(min(Pn.Prefv-4.*Pn.Puncv));
    plmax = ceil(max(Pn.Prefv+4.*Pn.Puncv));
    
    % vectors and matrices for probability estimation
    xv = linspace(plmin,plmax,(plmax-plmin)/step+1);
    cxm = repmat(xv,numel(cPrefv),1); if numel(siteP.include) == 0; cxm = []; end;
    cPrefm = repmat(cPrefv,1,numel(xv));
    cPuncm = repmat(cPuncv,1,numel(xv));
    oxm = repmat(xv,numel(oPrefv),1);
    oPrefm = repmat(oPrefv,1,numel(xv));
    oPuncm = repmat(oPuncv,1,numel(xv));
    
    % estimate probability distribution
    cprobdensmatr = normpdf(cxm,cPrefm,cPuncm);
    oprobdensmatr = normpdf(oxm,oPrefm,oPuncm);
    
    % pick out summed probability distribution
    if numel(siteP.include) > 0;
        probsum = sum(cprobdensmatr);
    else;
        probsum = sum(oprobdensmatr);
    end;
    
    % start new figure
    figure('name',[pltitle ' – ' nuclstr],'NumberTitle','off','visible','off');
    hold on; box on;
    
    % plot grey Pref uncertainty region
    uncP = linspace(siteP.Pref-siteP.Prefunc,siteP.Pref+siteP.Prefunc,100);
    uncprob = interp1(xv,probsum,uncP,'pchip');
    uncP = [uncP,siteP.Pref+siteP.Prefunc,siteP.Pref-siteP.Prefunc];
    uncprob = [uncprob 0 0];
    patch(uncP,uncprob,[0.85 0.85 0.85],'EdgeColor','none');
    
    % plot Pref line
    Prefprob = interp1(xv,probsum,siteP.Pref,'pchip');
    plot([siteP.Pref siteP.Pref],[0 Prefprob],'color','black');
    
    % plot individual sample prob dens curves (outliers)
    for j = 1:size(oprobdensmatr,1);
        plot(xv',oprobdensmatr(j,:)','color',[0.7 0.7 0.7]);
    end
    
    % plot individual sample prob dens curves (clustered samples)
    for j = 1:size(cprobdensmatr,1);
        plot(xv',cprobdensmatr(j,:)','color',clrin);
    end
    
    % plot summed prob dens curve
    plot(xv',probsum','color','black');
    
    xlabel(strcat({'Ref '},lablestr,{' prodrate (atoms/g/yr)'}));
    set(gca,'ytick',[]);
    ylabel('Relative probability');
    set(gca,'layer','top'); % plot axis on top
    title(pltitle);
    hold off;
% end subfunction plot_Pref ========================================================================


% subfunction plot_fullP ===========================================================================
function plot_fullP(Pfull,plfull,nucl);
    % fix nuclide specific parameters
    if nucl == 10; nuclstr = '10Be'; yv = [1.9 6];
    elseif nucl == 26; nuclstr = '26Al'; yv = [20 42];
    elseif nucl == 14; nuclstr = '14C'; yv = [5 23]; end;
    
    % fix sites and number of sites
    sites = Pfull.site;
    sites(cellfun(@isempty,sites)) = [];
    sites = unique(sites,'stable');
    numsites = numel(sites);

    % plot global Pref and uncertainty
    figure('name',['Pref – ' nuclstr],'NumberTitle','off');
    hold on; box on;
    Pmaxtot = Pfull.Ptot + Pfull.Ptotunc;
    Pmintot = Pfull.Ptot - Pfull.Ptotunc;
    plot([0 numsites numsites 0 0],[Pmaxtot Pmaxtot Pmintot Pmintot Pmaxtot],'color','black');
    plot([0 numsites],[Pfull.Ptot Pfull.Ptot],'color','black');

    % fix and plot site Punc
    idxin = find(Pfull.PsiteOK == 1);
    sitenr = (1:1:numsites)(idxin);
    siteP = Pfull.Psite(idxin);
    sitePunc = Pfull.Psiteunc(idxin);
    siteunc_x = [sitenr-1;sitenr;sitenr;sitenr-1;sitenr-1];
    Pmaxsite = siteP + sitePunc;
    Pminsite = siteP - sitePunc;
    siteunc_y = [Pmaxsite;Pmaxsite;Pminsite;Pminsite;Pmaxsite];
    plot(siteunc_x,siteunc_y,'color','blue');
    % fix and plot site Pref
    siteP_x = [sitenr-1;sitenr];
    siteP_y = [siteP;siteP];
    plot(siteP_x,siteP_y,'color','blue');

    % fix Px, P, Punc, Pin
    Px = []; P = []; Punc = []; Pin = [];
    Pincl = zeros(1,numel(Pfull.P));
    Pincl(Pfull.PidxOK) = 1;
    for i = 1:numsites;
        idx = find(strcmp(Pfull.site,sites(i)));
        n = numel(idx);
        Px(end+1:end+n) = i - 1 + linspace(1/(n+1),1-1/(n+1),n);
        P(end+1:end+n) = Pfull.P(idx);
        Punc(end+1:end+n) = Pfull.Punc(idx);
        Pin(end+1:end+n) = Pincl(idx);
    end;
    
    % fix included and excluded sample P
    idxin = find(Pin == 1);
    idxex = find(Pin == 0);
    P_in = P(idxin);
    Punc_in = Punc(idxin);
    P_ex = P(idxex);
    Punc_ex = Punc(idxex);
    
    % fix Pin and Puncin
    Pin_x = Px(idxin);
    Puncin_x = [Pin_x;Pin_x];
    Puncin_y = [P_in-Punc_in;P_in+Punc_in];
    % fix Pex and Puncex
    Pex_x = Px(idxex);
    Puncex_x = [Pex_x;Pex_x];
    Puncex_y = [P_ex-Punc_ex;P_ex+Punc_ex];
    % plot sample Puncex and Pex
    plot(Puncex_x,Puncex_y,'color','red');
    plot(Pex_x,P_ex,'marker','o','markersize',2,'linestyle','none','color','red');
    % plot sample Puncin and Pin
    plot(Puncin_x,Puncin_y,'color','black');
    plot(Pin_x,P_in,'marker','o','markersize',2,'linestyle','none','color','black');

    % fix axes
    if plfull.xmax <= 0; plfull.xmax = numsites; end;
    xlim([0 plfull.xmax]);
    ylim(yv);

    % fix x-ticks
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);

    % fix axis labels and tics
    if plfull.printsvg ~= 1;
        ylabel([nuclstr ' Pref (atoms/g/yr)']);
    else;
        set(gca,'yticklabel',[]);
        print(['Pref-' nuclstr '.svg'],'-dsvg',plfull.printscale);
    end;
% end subfunction plot_fullP =======================================================================
