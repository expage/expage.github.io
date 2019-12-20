function prodrate()

% Function for 10Be and 26Al reference production rate calibration.
% Calibrate site reference 10Be and/or 26Al production rate using chi-square minimization.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

clear;
close all;
tic();

% What version is this?
ver = '201912';

% do choices here ==================================================================================
% plotting / cluster test?
plotting = 1; % plot Pref probability curves (1 = yes)
Pcluster = 0; % exclude outliers to try to achieve a well-clustered group Pref (1 = yes)

% parameters for cluster test / outlier rejection
chiprob = 0.05; % lower p-value limit for chi-square probability test
mingroupn = 3; % minimum number of samples in well-clustered group
maxoutratio = 1/3; % maximum outlier ratio
% ==================================================================================================

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','isostP','isostsubm','site','lat','long','elv',...
    'thick','dens','shield','erosion','N10','N10unc','N26','N26unc','samplingyr','pressure',...
    'calage','calageunc'};
vartypes = {'%s','%s','%s','%s','%s','%s','%s','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n','%n','%n','%n','%n'};
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
indata = textscan(fid,typestr,'CommentStyle','%'); % scan data
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

% if there is no site in input: assume all samples are from one site
if isfield(samplein,'site') == 0; samplein.site(1:numel(samplein.sample),1) = {'noname'}; end;

% if there is no erosion in input: assume zero erosion
if isfield(samplein,'erosion') == 0; samplein.erosion(1:numel(samplein.sample),1) = 0; end;

% if there is no N10 or N26 in input: fill with 0
if isfield(samplein,'N10') == 0; samplein.N10(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N10unc') == 0; samplein.N10unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std10') == 0; samplein.std10(1:numel(samplein.sample),1) = {'0'}; end;
if isfield(samplein,'N26') == 0; samplein.N26(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N26unc') == 0; samplein.N26unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std26') == 0; samplein.std26(1:numel(samplein.sample),1) = {'0'}; end;

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

% fix for output
output(1,1) = {'sample'};
col.sample = 1;
if sum(samplein.N10+samplein.N10unc) > 0;
    output(1,2:3) = {'P10(at/g/yr)','P10unc(at/g/yr)'};
    col.P10 = col.sample+1; col.P10unc = col.P10+1;
    if Pcluster == 1;
        output(1,end+1) = {'P10-included?'};
        col.P10incl = col.P10unc+1;
    end;
end;
if sum(samplein.N26+samplein.N26unc) > 0;
    output(1,max(cell2mat(struct2cell(col)))+1:max(cell2mat(struct2cell(col)))+2) = ...
        {'P26(at/g/yr)','P26unc(at/g/yr)'};
    col.P26 = max(cell2mat(struct2cell(col)))+1; col.P26unc = col.P26+1;
    if Pcluster == 1;
        output(1,end+1) = {'P26-included?'};
        col.P26incl = col.P26unc+1;
    end;
end;

% fix for site output
outputsite(1,1) = {'site'};
scol.site = 1;
if sum(samplein.N10+samplein.N10unc) > 0;
    outputsite(1,2:5) = {'P10(at/g/yr)','P10unc(at/g/yr)','P10-χ2R','P10-p-value'};
    scol.P10 = scol.site+1; scol.P10unc = scol.P10+1;
    scol.chi10 = scol.P10unc+1; scol.pvalue10 = scol.chi10+1;
    if Pcluster == 1;
        outputsite(1,scol.pvalue10+1) = {'P10-included?'};
        scol.P10incl = scol.pvalue10+1;
    end;
end;
if sum(samplein.N26+samplein.N26unc) > 0;
    scolmax = max(cell2mat(struct2cell(scol)));
    outputsite(1,scolmax+1:scolmax+4) = {'P26(at/g/yr)','P26unc(at/g/yr)','P26-χ2R','P26-p-value'};
    scol.P26 = scolmax+1; scol.P26unc = scol.P26+1;
    scol.chi26 = scol.P26unc+1; scol.pvalue26 = scol.chi26+1;
    if Pcluster == 1;
        outputsite(1,scol.pvalue26+1) = {'P26-included?'};
        scol.P26incl = scol.pvalue26+1;
    end;
end;

% loop for single sites
while numel(samplein.site) > 0;
    % display site name
    fprintf(1,'\nSite: %s\n',samplein.site{1});
    
    % pick out all samples with same sampleID
    siteidx = find(strcmp(samplein.site,samplein.site(1)));
    
    % pick out variables
    sitein.sample = samplein.sample(siteidx);
    sitein.lat = samplein.lat(siteidx);
    sitein.long = samplein.long(siteidx);
    sitein.thick = samplein.thick(siteidx);
    sitein.dens = samplein.dens(siteidx);
    sitein.shield = samplein.shield(siteidx);
    sitein.erosion = samplein.erosion(siteidx);
    sitein.N10 = samplein.N10(siteidx);
    sitein.N10unc = samplein.N10unc(siteidx);
    sitein.N26 = samplein.N26(siteidx);
    sitein.N26unc = samplein.N26unc(siteidx);
    sitein.samplingyr = samplein.samplingyr(siteidx);
    sitein.calage = samplein.calage(siteidx);
    sitein.calageunc = samplein.calageunc(siteidx);
    sitein.pressure = samplein.pressure(siteidx);
    sitein.thickSF1 = samplein.thickSF1(siteidx);
    sitein.rown = samplein.rown(siteidx);
    if isfield(samplein,'isostP');
        sitein.isostP = samplein.isostP(siteidx);
        sitein.elv = samplein.elv(siteidx);
        sitein.Pflag = samplein.Pflag(siteidx);
    end;
    if isfield(samplein,'isostsubm');
        sitein.isostsubm = samplein.isostsubm(siteidx);
        sitein.elv = samplein.elv(siteidx);
    end;
    
    % calculate sample Prefs
    [output,P10,P26] = get_sample_Prefs(sitein,vars,output,col,consts);
    
    % calculate site Pref for 10Be
    if numel(P10.Prefv) >= mingroupn || (Pcluster==0 && numel(P10.Prefv)>1);
        % calculate site Pref
        siteP10 = get_site_Pref(P10);
        
        % if doing cluster analysis
        if Pcluster == 1;
            % fix variable cl for function get_cluster
            cl.maxoutratio = maxoutratio; cl.mingroupn = mingroupn; cl.chiprob = chiprob;
            
            % find clustered Pref by removing outliers
            siteP10 = get_cluster(P10,siteP10,cl);
            
            % mark OK samples in output
            output(siteP10.OKrow+1,col.P10incl) = {'X'};
        end;
    end;
    
    % calculate site Pref for 26Al
    if numel(P26.Prefv) >= mingroupn || (Pcluster==0 && numel(P26.Prefv)>1);
        % calculate site Pref
        siteP26 = get_site_Pref(P26);
        
        % if doing cluster analysis
        if Pcluster == 1;
            % fix variable cl for function get_cluster
            cl.maxoutratio = maxoutratio; cl.mingroupn = mingroupn; cl.chiprob = chiprob;
            
            % find clustered Pref by removing outliers
            siteP26 = get_cluster(P26,siteP26,cl);
            
            % mark OK samples in output
            output(siteP26.OKrow+1,col.P26incl) = {'X'};
        end;
    end;
    
    % fill outputsite and plot
    if numel(P10.Prefv)+numel(P26.Prefv) > 0;
        siterow = size(outputsite,1)+1;
        outputsite(siterow,scol.site) = samplein.site(1);
    end;
    if exist('siteP10');
        outputsite(siterow,scol.P10) = {num2str(siteP10.Pref,'%.3f')};
        outputsite(siterow,scol.P10unc) = {num2str(siteP10.Prefunc,'%.3f')};
        outputsite(siterow,scol.chi10) = {num2str(siteP10.rchisq,'%.3f')};
        outputsite(siterow,scol.pvalue10) = {num2str(siteP10.Pvalue,'%.3f')};
        % mark OK Pref in output
        if Pcluster==1 && siteP10.Pvalue>=chiprob && numel(siteP10.OKrow)>=mingroupn;
            outputsite(siterow,scol.P10incl) = {'X'};
        end;
        % display site Pref
        fprintf(1,'Site Pref10 = %s ± %s at/g/yr',outputsite{siterow,scol.P10},...
            outputsite{siterow,scol.P10unc});
        if Pcluster==1 && numel(siteP10.OKrow)<numel(P10.Prow);
            fprintf(1,'  (%.0f outlier',numel(P10.Prow)-numel(siteP10.OKrow));
            if numel(P10.Prow)-numel(siteP10.OKrow) > 1; fprintf(1,'s)'); else; fprintf(1,')'); end;
        end;
        fprintf(1,'\nR-chi2 = %s    P-value = %s\n',outputsite{siterow,scol.chi10},...
            outputsite{siterow,scol.pvalue10});
        % do plotting here...
        if plotting == 1;
            if isfield(siteP10,'include') == 0; siteP10.include = (1:1:numel(P10.Prefv)); end;
            plot_Pref(P10,siteP10,10,samplein.site{1});
        end;
    end;
    if exist('siteP26');
        outputsite(siterow,scol.P26) = {num2str(siteP26.Pref,'%.3f')};
        outputsite(siterow,scol.P26unc) = {num2str(siteP26.Prefunc,'%.3f')};
        outputsite(siterow,scol.chi26) = {num2str(siteP26.rchisq,'%.3f')};
        outputsite(siterow,scol.pvalue26) = {num2str(siteP26.Pvalue,'%.3f')};
        % mark OK Pref in output
        if Pcluster==1 && siteP26.Pvalue>=chiprob && numel(siteP26.OKrow)>=mingroupn;
            outputsite(siterow,scol.P26incl) = {'X'};
        end;
        % display site Pref
        fprintf(1,'Site Pref26 = %s ± %s at/g/yr',outputsite{siterow,scol.P26},...
            outputsite{siterow,scol.P26unc});
        if Pcluster==1 && numel(siteP26.OKrow)<numel(P26.Prow);
            fprintf(1,'  (%.0f outlier',numel(P26.Prow)-numel(siteP26.OKrow));
            if numel(P26.Prow)-numel(siteP26.OKrow) > 1; fprintf(1,'s)'); else; fprintf(1,')'); end;
        end;
        fprintf(1,'\nR-chi2 = %s    P-value = %s\n',outputsite{siterow,scol.chi26},...
            outputsite{siterow,scol.pvalue26});
        % do plotting here...
        if plotting == 1;
            if isfield(siteP26,'include') == 0; siteP26.include = (1:1:numel(P26.Prefv)); end;
            plot_Pref(P26,siteP26,26,samplein.site{1});
        end;
    end;
    
    % remove calculated samples
    samplein.site(siteidx) = [];
    samplein.sample(siteidx) = [];
    samplein.lat(siteidx) = [];
    samplein.long(siteidx) = [];
    samplein.thick(siteidx) = [];
    samplein.dens(siteidx) = [];
    samplein.shield(siteidx) = [];
    samplein.erosion(siteidx) = [];
    samplein.N10(siteidx) = [];
    samplein.N10unc(siteidx) = [];
    samplein.N26(siteidx) = [];
    samplein.N26unc(siteidx) = [];
    samplein.samplingyr(siteidx) = [];
    samplein.calage(siteidx) = [];
    samplein.calageunc(siteidx) = [];
    samplein.pressure(siteidx) = [];
    samplein.thickSF1(siteidx) = [];
    samplein.rown(siteidx) = [];
    samplein.elv(siteidx) = [];
    samplein.Pflag(siteidx) = [];
    if isfield(samplein,'isostP');
        samplein.isostP(siteidx) = [];
    end;
    if isfield(samplein,'isostsubm');
        samplein.isostsubm(siteidx) = [];
    end;
    
    % clear variables
    clear sitein; clear P10; clear P26; clear siteP10; clear siteP26;
end;

% display plots
plotv = findobj('type','figure');
for i = 1:numel(plotv); figure(plotv(i)); end;

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
        for i = 1:size(outputsite,1);
            fprintf(out,outstr2,outputsite{i,:});
        end;
    end;
    fclose(out);
end;
% ================================================

toc()
clear;
% end prodrate function ============================================================================


% subfunction get_sample_Prefs =====================================================================
function [output,P10,P26] = get_sample_Prefs(sitein,vars,output,col,consts);
    % number of 10Be and 26Al samples
    num10 = 0;
    num26 = 0;
    
    % declare Pref and Prow matrices
    P10.Prefv = []; P10.Puncv = []; P10.Prow = [];
    P26.Prefv = []; P26.Puncv = []; P26.Prow = [];
    
    for i = 1:numel(sitein.sample);
        % pick out single sample data
        sample.sample = sitein.sample{i};
        sample.lat = sitein.lat(i);
        sample.long = sitein.long(i);
        sample.thick = sitein.thick(i);
        sample.dens = sitein.dens(i);
        sample.shield = sitein.shield(i);
        sample.erosion = sitein.erosion(i);
        sample.N10 = sitein.N10(i);
        sample.N10unc = sitein.N10unc(i);
        sample.N26 = sitein.N26(i);
        sample.N26unc = sitein.N26unc(i);
        sample.samplingyr = sitein.samplingyr(i);
        sample.calage = sitein.calage(i);
        sample.calageunc = sitein.calageunc(i);
        sample.pressure = sitein.pressure(i);
        sample.thickSF1 = sitein.thickSF1(i);
        sample.rown = sitein.rown(i);
        if isfield(sitein,'isostP');
            sample.isostP = sitein.isostP{i};
            sample.elv = sitein.elv(i);
            sample.Pflag = sitein.Pflag{i};
            if strcmp(sample.isostP,'-') || strcmp(sample.Pflag,'pre');
                sample = rmfield(sample,'isostP');
            end;
        end;
        if isfield(sitein,'isostsubm');
            sample.isostsubm = sitein.isostsubm{i};
            sample.elv = sitein.elv(i);
            if strcmp(sample.isostsubm,'-'); sample = rmfield(sample,'isostsubm'); end;
        end;
        
        % write sample name to output
        output(sample.rown+1,1) = sample.sample;
        
        % Set nucl and mt to 0 for both 10/26 and check if there is N10/N26
        nucl10 = 0; nucl26 = 0; mt10 = 0; mt26 = 0;
        if (sample.N10 + sample.N10unc) > 0; nucl10 = 1; num10 = num10 + 1; end;
        if (sample.N26 + sample.N26unc) > 0; nucl26 = 1; num26 = num26 + 1; end;
        
        if nucl10 + nucl26 == 0;
            continue;
        end;
        
        % Pref for simple age estimates
        Pref10 = consts.Pref10;
        Pref26 = consts.Pref26;
        
        % change Pref to isostatic calibration values if using isostatic adjustment
        if isfield(sample,'isostP');
            Pref10 = consts.Pref10iso;
            Pref26 = consts.Pref26iso;
        end;
        
        % display sample name
        fprintf(1,'%.0f. %s',i,sample.sample);
        
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
        
        % use calibration age for scaling
        mt = sample.calage;
        
        % Age Relative to t0=2010 - LSD tv from LSDfix
        % tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
        
        % Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
        LSDfix = LSD_fix(sample.lat,sample.long,mt,-1,sample.samplingyr,consts);
        
        % include tv in sample
        sample.tv = LSDfix.tv;
        
        % fix submergence vectors if using isostatic adjustment
        if isfield(sample,'isostsubm');
            sample.elv = isost_elv(sample.isostsubm,sample);
            sample.overwater = (sample.elv >= 0);
            sample.underwater = (sample.elv < 0);
            sample.waterdepth = -sample.elv .* sample.underwater .* 1E2; % cm (dens assumed to be 1)
            if sum(sample.underwater) == 0;
                sample = rmfield(sample,'isostsubm');
            end;
        end;
        
        % Production from muons
        if sample.erosion <= 0 && isfield(sample,'isostsubm') == 0;
            Pmu = P_mu_expage(sample.thick.*sample.dens./2,sample.pressure,LSDfix.RcEst,...
                consts.SPhiInf,nucl10,nucl26,consts,'no');
            if nucl10 == 1; sample.mu10 = Pmu.mu10 .* sample.shield; end;
            if nucl26 == 1; sample.mu26 = Pmu.mu26 .* sample.shield; end;
        else;
            tv_z = (sample.tv.*sample.erosion + sample.thick./2) .* sample.dens; % T-d vect (g/cm^2)
            if nucl10 == 1;
                sample.mu10 = get_PmuE(sample,tv_z,tsimple10,LSDfix.RcEst,consts,1,0);
            end;
            if nucl26 == 1;
                sample.mu26 = get_PmuE(sample,tv_z,tsimple26,LSDfix.RcEst,consts,0,1);
            end;
        end;
        
        % fix atmospheric pressure if using isostatic adjustment
        if isfield(sample,'isostP');
            % calculate elevation development
            sample.elvv = isost_elv(sample.isostP,sample);
            % calculate atmospheric pressure
            if strcmp(sample.Pflag,'std');
                sample.pressure = ERA40atm(sample.lat,sample.long,sample.elvv);
            elseif strcmp(sample.Pflag,'ant');
                sample.pressure = antatm(sample.elvv);
            end;
        end;
        
        % spallation production scaling
        Psp = P_sp_expage(sample.pressure,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,nucl10,nucl26);
        
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
        
        if nucl10 == 1;
            % various parameters
            sample.dcf = exp(-sample.tv.*vars.l10);
            sample.N = sample.N10; sample.Nunc = sample.N10unc;
            sample.l = vars.l10;
            sample.sp = Psp.sp10;
            sample.mu = sample.mu10;
            
            % get sample Pref
            [Prefi,Prefiunc,Prefiuncint] = Pref_calc(sample,'10');
            
            % fill output
            output(sample.rown+1,col.P10:col.P10unc) = {num2str(Prefi,'%.3f'),...
                num2str(Prefiunc,'%.3f')};
            
            % display Pref
            fprintf(1,' \tP10 = %s ± %s at/g/yr',output{sample.rown+1,col.P10},...
                output{sample.rown+1,col.P10unc});
            
            % fill Pref vector for uncertainty estimation, plotting, and cluster analysis
            P10.Prefv(num10) = Prefi;
            P10.Puncv(num10) = Prefiuncint;
            P10.Prow(num10) = sample.rown; % input row number
            P10.calage(num10) = sample.calage;
            P10.calageunc(num10) = sample.calageunc;
        end;
        
        if nucl26 == 1;
            % various parameters
            sample.dcf = exp(-sample.tv.*vars.l26);
            sample.N = sample.N26; sample.Nunc = sample.N26unc;
            sample.l = vars.l26;
            sample.sp = Psp.sp26;
            sample.mu = sample.mu26;
            
            % get sample Pref
            [Prefi,Prefiunc,Prefiuncint] = Pref_calc(sample,'26');
            
            % fill output
            output(sample.rown+1,col.P26:col.P26unc) = {num2str(Prefi,'%.3f'),...
                num2str(Prefiunc,'%.3f')};
            
            % display Pref
            fprintf(1,' \tP26 = %s ± %s at/g/yr',output{sample.rown+1,col.P26},...
                output{sample.rown+1,col.P26unc});
            
            % fill Pref vector for uncertainty estimation, plotting, and cluster analysis
            P26.Prefv(num26) = Prefi;
            P26.Puncv(num26) = Prefiuncint;
            P26.Prow(num26) = sample.rown; % input row number
            P26.calage(num26) = sample.calage;
            P26.calageunc(num26) = sample.calageunc;
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
function out = get_PmuE(sample,tv_z,tsimple,RcEst,consts,nucl10,nucl26);
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
    Pmu_d = P_mu_expage(mu_z,sample.pressure,RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
    if nucl10 == 1; Pmud = Pmu_d.mu10; elseif nucl26 == 1; Pmud = Pmu_d.mu26; end;
    out = interp1(mu_z,Pmud,tv_z,'pchip') .* sample.shield; % P_mu
    out(isnan(out)) = out(end); % fix for nan issue...
% end subfunction get_PmuE =========================================================================


% subfunction Pref_calc ============================================================================
function [Pref,Prefunc,Prefuncint] = Pref_calc(sample,nuclstr);
    % sample spallation production rate scaling over time including decay and erosion
    Psp = sample.sp.*sample.dcf.*sample.dpfs.*sample.thickSF.*sample.shield;
    
    % sample muon P
    Pmu = sample.mu.*sample.dcf;
    
    % calculate N produced by muons
    Nmu = trapz(sample.tv,Pmu);
    
    % calculate N without Pref produced by spallation
    Nsp_scaling = trapz(sample.tv,Psp);
    
    % calculate Pref
    Pref = (sample.N-Nmu)./Nsp_scaling;
    
    % age error propagation from Balco et al. (2008) CRONUS calculator
    % A with decay-weighted average Lsp
    Lsp_avg = trapz(sample.tv,sample.Lsp.*exp(-sample.l.*sample.tv))/...
        trapz(sample.tv,exp(-sample.l.*sample.tv));
    A = sample.l + sample.dens.*sample.erosion./Lsp_avg;
    FP = (sample.N.*A)./(1-exp(-A.*sample.calage));
    dtdN = 1./(FP-sample.N.*A);
    ageunc = sqrt(dtdN.^2.*sample.Nunc.^2);
    
    % calculate internal and external (including calibration age uncertainty) Pref uncertainty
    Prefuncint = ageunc./sample.calage.*Pref;
    Prefunc = sqrt(Prefuncint.^2 + (sample.calageunc./sample.calage.*Pref).^2);
% end subfunction Pref_calc ========================================================================


% subfunction get_site_Pref ========================================================================
function siteP = get_site_Pref(Pn);
    % calculate weighted Pref and internal uncertainty
    [siteP.Pref,Puncint] = evm(Pn.Prefv,Pn.Puncv);
    
    % number of samples
    nn = numel(Pn.Prefv);
    
    % calculate reduced chi square
    siteP.rchisq = 1/(nn-1) .* sum(((Pn.Prefv-siteP.Pref)./Pn.Puncv).^2);
    
    % calculate P value
    siteP.Pvalue = 1 - chi2cdf(siteP.rchisq.*(nn-1),nn-1);
    
    % calculate calibration age uncertainty part
    calage_unc = mean(Pn.calageunc./Pn.calage);
    
    % calculate total prod rate uncertainty
    siteP.Prefunc = sqrt(Puncint^2 + (calage_unc.*siteP.Pref)^2);
% end subfunction get_site_Pref ====================================================================


% subfunction evm ==================================================================================
function [Pref,Punc] = evm(Prefv,Puncv);
    % fix data (make matrices)
    Prefm = repmat(Prefv,numel(Prefv),1);
    Puncm = repmat(Puncv,numel(Puncv),1);
    xm = repmat(Prefv',1,numel(Prefv));
    
    % calculate probability for each sample Pref
    Mui = sum(sqrt(2./(pi.*(2.*Puncm).^2)).*exp(-(xm-Prefm).^2./(2.*Puncm.^2)))./numel(Prefv);
    
    % calculate summed probability for all samples
    Muj = sum(Mui);
    
    % calculate sample weights
    wi = Mui./Muj;
    
    % calculate weighted mean Pref
    Pref = sum(wi.*Prefv);
    
    % uncertainty estimation
    uncint = sqrt(sum(wi.^2.*Puncv.^2));
    uncext = sqrt(sum(wi.*(Prefv-Pref).^2));
    Punc = max(uncint,uncext);
% end subfunction evm ==============================================================================



% subfunction get_cluster ==========================================================================
function siteP = get_cluster(Pn,siteP,cl);
    % in variables
    sitePin = siteP;
    
    % define maximum number of outliers
    remove = floor(numel(Pn.Prefv)*cl.maxoutratio);
    r = 1; % outlier counter
    
    % well-clustered samples for plotting
    Pn.include = (1:1:numel(Pn.Prow));
    
    % loop for outlier removal
    while siteP.Pvalue<cl.chiprob && remove>=r && numel(Pn.Prefv)>=cl.mingroupn;
        % calculate deviation for all individual samples
        chidev = ((Pn.Prefv-siteP.Pref)./Pn.Puncv).^2;
        
        % find index of sample with largest dev
        [value,rmv_idx] = max(chidev);
        
        % remove outlier
        Pn.calage(rmv_idx) = []; Pn.calageunc(rmv_idx) = [];
        Pn.Prefv(rmv_idx) = []; Pn.Puncv(rmv_idx) = [];
        Pn.Prow(rmv_idx) = []; Pn.include(rmv_idx) = [];
        
        % calculate Pref
        siteP = get_site_Pref(Pn);
        
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


% subfunction plot_Pref ============================================================================
function plot_Pref(Pn,siteP,nucl,pltitle);
    % fix for 10Be/26Al    
    if nucl == 10;
        step = 0.01;
        nuclstr = '^{10}Be';
    elseif nucl == 26;
        step = 0.02;
        nuclstr = '^{26}Al';
    end;
    
    % split Pref vectors in clustered and outliers
    cPrefv = Pn.Prefv(siteP.include)';
    cPuncv = Pn.Puncv(siteP.include)';
    oPrefv = Pn.Prefv'; oPrefv(siteP.include) = [];
    oPuncv = Pn.Puncv'; oPuncv(siteP.include) = [];
    
    % define min and max Pref for plot
    plmin = floor(min(Pn.Prefv-4.*Pn.Puncv));
    plmax = ceil(max(Pn.Prefv+4.*Pn.Puncv));
    %~ plmin = floor(siteP.Pref-4.*siteP.Prefunc);
    %~ plmax = ceil(siteP.Pref+4.*siteP.Prefunc);
    
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
    figure,set(gcf,'visible','off'); hold on; box on;
    
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
        plot(xv',oprobdensmatr(j,:)','color','blue');
    end
    
    % plot individual sample prob dens curves (clustered samples)
    for j = 1:size(cprobdensmatr,1);
        plot(xv',cprobdensmatr(j,:)','color','red');
    end
    
    % plot summed prob dens curve
    plot(xv',probsum','color','black');
    
    xlabel(strcat({'Ref '},nuclstr,{' prodrate (atoms/g/yr)'}));
    set(gca,'ytick',[]);
    ylabel('Relative probability');
    set(gca,'layer','top'); % plot axis on top
    title(pltitle);
    hold off;
% end subfunction plot_Pref ========================================================================
