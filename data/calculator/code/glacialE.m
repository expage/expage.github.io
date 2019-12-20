function glacialE();

% 10Be and 26Al glacial erosion calculator.
% Function for quantification/investigation of glacial erosion based on 10Be and/or 26Al conc.
% Interpolates the best fitting glacial erosion to yield the 10Be and/or 26Al conc given as input
% based on the time-dependent cosmogenic nuclide production rate.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman (jakob.heyman@gu.se) 2018-2019

clear; close all; tic();

% What version is this?
ver = '201912';

% =============== MAKE CHOICES HERE ================================================================
% max time (determines the start of the simulation - max 1E7)
mt = 1E7;

% calculate glacial erosion rate and/or incremental erosion depth steps? 1 = yes, 0 = no
glacErate = 1;
glacEstep = 0;

% Monte Carlo iterations for uncertainty estimation
mc = 1E3;

% calculate combined 10Be-26Al erosion? 1 = yes, 0 = no
nucl1026 = 1;

% plotting choices
plch.depth = 1;         % plot sample depth
plch.P = 0;             % plot sample production - NOT AVAILABLE YET
plch.N = 0;             % plot nuclide build-up - NOT AVAILABLE YET
plch.combined = 0;      % plot 10Be and 26Al data together?
plch.uncline = 0;       % plot uncertainties as lines instead of areas?
plch.banana = 0;        % plot 26/10 banana - NOT AVAILABLE YET
plch.banana_path = 10;  % number of OK sample 26Al/10Be ratio paths to plot in banana
plch.maxt = 1E5;        % max time for plotting (yr)
if plch.maxt > mt; plch.maxt = mt; end; % don't allow too large maxt
plch.maxt = plch.maxt.*1E-3; % change to ka
plch.maxd = 10;         % max depth for plotting (m) - comment out to not use
plch.clr10 = [1 0 0];   % color for 10Be: red
plch.clr26 = [0 0 1];   % color for 26Al: blue
plch.clr1026 = [0 0 0]; % color for combined 10Be+26Al: black

% time vector (yr) for writing sample depth (m) in output
tdv = [1E5 1E6];
tdv(tdv>mt) = []; % remove time points beyond max time

% fixed parameters for specific time period(s) =======================
% use multiple rows for multiple periods - comment out if not using
%fixE.nonglacE_tv = [0 200000]; % yr  [min max]
%fixE.nonglacE = 4; % mm/ka - non-glacial erosion rate
%fixE.glacEr_tv = [0 30000;100000 300000]; % yr  [min max]
%fixE.glacEr = [0;10]; % mm/ka - glacial erosion rate
%fixE.glacEi_tv = [10000 11000]; % yr  [min max]
%fixE.glacEi = [80]; % cm/glac - incremental glacial erosion steps
%fixE.fix0yr = 1950; % start year for fixE vectors

% fixed ice-cover and icefree for specific time period(s) ============
% use multiple rows for multiple periods - comment out if not using
%fixice.ice = [40000 70000]; % yr  [min max]
%fixice.noice = [10000 23000]; % yr  [min max]
%fixice.fix0yr = 1950; % start year for fixice vectors
% ====================================================================

% ice cover proxy file - default is Lisiecki and Raymo (2005) d18O ice volume record
% column 1: year   column 2: proxy value
iceproxy = 'LR04.txt';
iceproxy_startyr = 1950; % start year (= 0) in iceproxy file

% glacial erosion vector
Etestv = [(0:0.01:0.99) logspace(0,4,900)]; % mm/ka or cm/glac
% ==================================================================================================

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','isostsubm','isostP','lat','long','elv','thick',...
    'dens','densunc','shield','erosion','erosionunc','N10','N10unc','N26','N26unc','samplingyr',...
    'pressure','deglac','deglacunc','icevalue','icevalueunc','glacE','glacEunc','isostsubmunc',...
    'burialdepth'};
vartypes = {'%s','%s','%s','%s','%s','%s','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n'};

% read input file
fid = fopen('input.txt');
varsin = strsplit(fgetl(fid)); % read first line
if(ismember(varsin,varnames)); % if first line contain only variable names
    [testi,vari] = ismember(varsin,varnames); % find index of varnames
    typestr = vartypes{vari(1)}; % fix type string
    for i = 2:numel(vari); % fix type string
        typestr = [typestr ' ' vartypes{vari(i)}];
    end;
else; % if no variable names in first line
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

% Decay constant
l10 = consts.l10; l10unc = consts.l10unc;
l26 = consts.l26; l26unc = consts.l26unc;

% if there is no N10 or N26 in input: fill with 0
if isfield(samplein,'N10') == 0; samplein.N10(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N10unc') == 0; samplein.N10unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std10') == 0; samplein.std10(1:numel(samplein.sample),1) = {'0'}; end;
if isfield(samplein,'N26') == 0; samplein.N26(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N26unc') == 0; samplein.N26unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std26') == 0; samplein.std26(1:numel(samplein.sample),1) = {'0'}; end;

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

% read iceproxy file
fid = fopen(iceproxy);
inproxy = textscan(fid,'%n %n','CommentStyle','%');
iceproxy_tv = inproxy{1};
iceproxyin = inproxy{2};
fclose(fid);
% extrapolate end values to 10 Ma and youngest sampling year
if iceproxy_tv(end)+2010-iceproxy_startyr < 1E7;
    iceproxy_tv = [iceproxy_tv; 1E7+iceproxy_startyr-2010];
    iceproxyin = [iceproxyin; iceproxyin(end)];
end;
if iceproxy_startyr < max(samplein.samplingyr);
    iceproxy_tv = [iceproxy_startyr-max(samplein.samplingyr); iceproxy_tv];
    iceproxyin = [iceproxyin(1); iceproxyin];
end;

% fix isostatic adjustment input
if isfield(samplein,'isostP');
    samplein.isostP = isost_data_load(samplein.isostP);
end;
if isfield(samplein,'isostsubm');
    samplein.isostsubm = isost_data_load(samplein.isostsubm);
end;

% Muon shielding depth vector to ~200,000 g/cm2
z_mu = [(0:3:27) logspace(1.48,5.3,70)];

% glacial erosion matrices
glacErm = [];
glacEim = [];

% for plotting
plnum.d10r = 0; plnum.d26r = 0; plnum.d1026r = 0; plnum.d10i = 0; plnum.d26i = 0; plnum.d1026i = 0;
if plch.depth == 1;
    if glacErate==1 && sum(samplein.N10)>0;
        plnum.d10r = 1;
        if sum(samplein.N26)>0 && plch.combined==1;
            plnum.d26r = plnum.d10r;
        end;
    end;
    if glacErate==1 && sum(samplein.N26)>0 && plnum.d26r==0;
        plnum.d26r = max(cell2mat(struct2cell(plnum)))+1;
    end;
    if glacErate==1 && sum(samplein.N10.*samplein.N26)>0 && nucl1026==1;
        if plch.combined==1; plnum.d1026r = plnum.d10r;
        else; plnum.d1026r = max(cell2mat(struct2cell(plnum)))+1; end;
    end;
    if glacEstep==1 && sum(samplein.N10)>0;
        plnum.d10i = max(cell2mat(struct2cell(plnum)))+1;
        if sum(samplein.N26)>0 && plch.combined==1;
            plnum.d26i = plnum.d10i;
        end;
    end;
    if glacEstep==1 && sum(samplein.N26)>0 && plnum.d26i==0;
        plnum.d26i = max(cell2mat(struct2cell(plnum)))+1;
    end;
    if glacEstep==1 && sum(samplein.N10.*samplein.N26)>0 && nucl1026==1;
        if plch.combined==1; plnum.d1026i = plnum.d10i;
        else; plnum.d1026i = max(cell2mat(struct2cell(plnum)))+1; end;
    end;
end;
pl10r.tm = []; pl10r.tmunc = []; pl10r.mm = []; pl10r.mmunc = []; pl10r.dm = []; pl10r.dmunc = [];
pl10i.tm = []; pl10i.tmunc = []; pl10i.mm = []; pl10i.mmunc = []; pl10i.dm = []; pl10i.dmunc = [];
pl26r.tm = []; pl26r.tmunc = []; pl26r.mm = []; pl26r.mmunc = []; pl26r.dm = []; pl26r.dmunc = [];
pl26i.tm = []; pl26i.tmunc = []; pl26i.mm = []; pl26i.mmunc = []; pl26i.dm = []; pl26i.dmunc = [];
pl1026r.tm = []; pl1026r.tmunc = []; pl1026r.mm = []; pl1026r.mmunc = [];
pl1026r.dm = []; pl1026r.dmunc = [];
pl1026i.tm = []; pl1026i.tmunc = []; pl1026i.mm = []; pl1026i.mmunc = [];
pl1026i.dm = []; pl1026i.dmunc = [];
pl10r.legloc = 'northwest'; pl26r.legloc = 'northwest'; pl1026r.legloc = 'northwest';
pl10i.legloc = 'northwest'; pl26i.legloc = 'northwest'; pl1026i.legloc = 'northwest';

% check if input contains both erosion and glacE
if isfield(samplein,'erosion') && isfield(samplein,'glacE');
    fprintf(1,'Input contains both glacial and nonglacial erosion - ');
    fprintf(1,'there is nothing to calculate!\n');
    return;
end;

% check if erosion calculation should be done for glacial erosion or non-glacial erosion
if isfield(samplein,'glacE');
    glaccalc = 0; % calculate nonglacial erosion
else;
    glaccalc = 1; % calculate glacial erosion
end;

% fix for time depth points headers
for j = 1:numel(tdv);
    if tdv(j) < 1E3;
        tstr = [num2str(tdv(j),'%.0f') 'yr'];
    elseif tdv(j) < 1E6;
        tstr = [num2str(tdv(j).*1E-3,'%.0f') 'ka'];
    elseif tdv(j).*1E-6 == floor(tdv(j).*1E-6);
        tstr = [num2str(tdv(j).*1E-6,'%.0f') 'Ma'];
    else;
        tstr = [num2str(tdv(j).*1E-6,'%.1f') 'Ma'];
    end;
    tdheader(j*3-2:j*3) = {[tstr '(m)'],[tstr '+(m)'],[tstr '-(m)']};
end;

% fix for output and plotting
output(1,1) = {'sample'};
outcol.sample = 1;
% 10Be output
if sum(samplein.N10)>0 && glacErate==1; % if using erosion rate
    output(1,end+1:end+3) = {'E10(mm/ka)','E10+(mm/ka)','E10-(mm/ka)'};
    outcol.E10r = outcol.sample+1; outcol.E10rp = outcol.E10r+1; outcol.E10rn = outcol.E10rp+1;
    if numel(tdv)>0;
        output(1,end+1:end+numel(tdheader)) = tdheader;
        outcol.E10rd1 = outcol.E10rn+1; outcol.E10rd2 = outcol.E10rn+numel(tdheader);
    end;
end;
if sum(samplein.N10)>0 && glacEstep==1; % if using erosion steps
    if glaccalc == 1;
        output(1,end+1:end+3) = {'E10(cm/glac)','E10+(cm/glac)','E10-(cm/glac)'};
    else;
        output(1,end+1:end+3) = {'E10i(mm/ka)','E10i+(mm/ka)','E10i-(mm/ka)'};
    end;
    outcol.E10i = max(cell2mat(struct2cell(outcol)))+1;
    outcol.E10ip = outcol.E10i+1; outcol.E10in = outcol.E10ip+1;
    if numel(tdv)>0;
        output(1,end+1:end+numel(tdheader)) = tdheader;
        outcol.E10id1 = outcol.E10in+1; outcol.E10id2 = outcol.E10in+numel(tdheader);
    end;
end;
% 26Al output
if sum(samplein.N26)>0 && glacErate==1; % if using erosion rate
    output(1,end+1:end+3) = {'E26(mm/ka)','E26+(mm/ka)','E26-(mm/ka)'};
    outcol.E26r = max(cell2mat(struct2cell(outcol)))+1;
    outcol.E26rp = outcol.E26r+1; outcol.E26rn = outcol.E26rp+1;
    if numel(tdv)>0;
        output(1,end+1:end+numel(tdheader)) = tdheader;
        outcol.E26rd1 = outcol.E26rn+1; outcol.E26rd2 = outcol.E26rn+numel(tdheader);
    end;
end;
if sum(samplein.N26)>0 && glacEstep==1; % if using erosion steps
    if glaccalc == 1;
        output(1,end+1:end+3) = {'E26(cm/glac)','E26+(cm/glac)','E26-(cm/glac)'};
    else;
        output(1,end+1:end+3) = {'E26i(mm/ka)','E26i+(mm/ka)','E26i-(mm/ka)'};
    end;
    outcol.E26i = max(cell2mat(struct2cell(outcol)))+1;
    outcol.E26ip = outcol.E26i+1; outcol.E26in = outcol.E26ip+1;
    if numel(tdv)>0;
        output(1,end+1:end+numel(tdheader)) = tdheader;
        outcol.E26id1 = outcol.E26in+1; outcol.E26id2 = outcol.E26in+numel(tdheader);
    end;
end;
% 10Be+26Al output
if sum(samplein.N10.*samplein.N26)>0 && nucl1026==1 && glacErate==1; % if using erosion rate
    output(1,end+1:end+3) = {'E1026(mm/ka)','E1026+(mm/ka)','E1026-(mm/ka)'};
    outcol.E1026r = max(cell2mat(struct2cell(outcol)))+1;
    outcol.E1026rp = outcol.E1026r+1; outcol.E1026rn = outcol.E1026rp+1;
    if numel(tdv)>0;
        output(1,end+1:end+numel(tdheader)) = tdheader;
        outcol.E1026rd1 = outcol.E1026rn+1; outcol.E1026rd2 = outcol.E1026rn+numel(tdheader);
    end;
end;
if sum(samplein.N10.*samplein.N26)>0 && nucl1026==1 && glacEstep==1; % if using erosion steps
    if glaccalc == 1;
        output(1,end+1:end+3) = {'E1026(cm/glac)','E1026+(cm/glac)','E1026-(cm/glac)'};
    else;
        output(1,end+1:end+3) = {'E1026i(mm/ka)','E1026i+(mm/ka)','E1026i-(mm/ka)'};
    end;
    outcol.E1026i = max(cell2mat(struct2cell(outcol)))+1;
    outcol.E1026ip = outcol.E1026i+1; outcol.E1026in = outcol.E1026ip+1;
    if numel(tdv)>0;
        output(1,end+1:end+numel(tdheader)) = tdheader;
        outcol.E1026id1 = outcol.E1026in+1; outcol.E1026id2 = outcol.E1026in+numel(tdheader);
    end;
end;

% pick out samples one by one
for i = 1:numel(samplein.lat);
    sample.sample = samplein.sample(i);
    sample.lat = samplein.lat(i);
    sample.long = samplein.long(i);
    sample.pressure = samplein.pressure(i);
    sample.thick = samplein.thick(i);
    sample.dens = samplein.dens(i);
    if isfield(samplein,'densunc'); sample.densunc = samplein.densunc(i); end;
    sample.shield = samplein.shield(i);
    sample.N10 = samplein.N10(i);
    sample.N10unc = samplein.N10unc(i);
    sample.N26 = samplein.N26(i);
    sample.N26unc = samplein.N26unc(i);
    sample.samplingyr = samplein.samplingyr(i);
    sample.icevalue = samplein.icevalue(i);
    if isfield(samplein,'icevalueunc'); sample.icevalueunc = samplein.icevalueunc(i); end;
    if isfield(samplein,'erosion');
        sample.erosion = samplein.erosion(i);
        if isfield(samplein,'erosionunc'); sample.erosionunc = samplein.erosionunc(i); end;
    else;
        sample.erosion = 0;
    end;
    if isfield(samplein,'glacE');
        sample.glacE = samplein.glacE(i);
        if isfield(samplein,'glacEunc'); sample.glacEunc = samplein.glacEunc(i); end;
        sample.erosion = Etestv;
    else;
        sample.glacE = Etestv;
    end;
    if isfield(samplein,'deglac');
        sample.deglac = samplein.deglac(i);
        if isfield(samplein,'deglacunc'); sample.deglacunc = samplein.deglacunc(i); end;
    end;
    if isfield(samplein,'isostP');
        sample.isostP = samplein.isostP{i};
        sample.elv = samplein.elv(i);
        sample.Pflag = samplein.Pflag(i);
        if strcmp(sample.isostP,'-') || strcmp(sample.Pflag,'pre');
            sample = rmfield(sample,'isostP');
        end;
    end;
    if isfield(samplein,'isostsubm');
        sample.isostsubm = samplein.isostsubm{i};
        sample.elv = samplein.elv(i);
        sample.Pflag = samplein.Pflag(i);
        if isfield(samplein,'isostsubmunc'); sample.isostsubmunc = samplein.isostsubmunc(i); end;
        if strcmp(sample.isostsubm,'-') || strcmp(sample.Pflag,'pre');
            sample = rmfield(sample,'isostsubm');
        end;
    end;
    if isfield(samplein,'burialdepth');
        sample.burialdepth = samplein.burialdepth(i);
    end;
    
    % write sample name to output
    output(i+1,1) = sample.sample;
    
    % Set nucl to 0 for both 10/26 and check if there is N10/N26
    nucl10 = 0; nucl26 = 0;
    if sample.N10 > 0; nucl10 = 1; end;
    if sample.N26 > 0; nucl26 = 1; end;
    
    % if no 10Be or 26Al: move on
    if nucl10+nucl26 == 0;
        continue;
    end;
    
    % Prefs
    Pref10 = consts.Pref10; Pref10unc = consts.Pref10unc;
    Pref26 = consts.Pref26; Pref26unc = consts.Pref26unc;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample{1});
    
    % Age Relative to t0=2010 - LSD tv from LSD_fix
    % tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
    
    % Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
    LSDfix = LSD_fix(sample.lat,sample.long,mt,-1,sample.samplingyr,consts);
    
    % include tv in sample
    sample.tv = LSDfix.tv(:);
    
    % fix ice proxy
    iceproxy = interp1(iceproxy_tv,iceproxyin,sample.tv+iceproxy_startyr-sample.samplingyr);
    
    % Production from muons along shielding depth vector z_mu
    Pmu_z = P_mu_expage(z_mu,sample.pressure,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
    
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
        % change Pref to isostatic calibration values
        Pref10 = consts.Pref10iso; Pref10unc = consts.Pref10isounc;
        Pref26 = consts.Pref26iso; Pref26unc = consts.Pref26isounc;
    end;
    
    % fix submergence vectors if using isostatic adjustment
    if isfield(sample,'isostsubm');
        sample.submelv = isost_elv(sample.isostsubm,sample)';
        sample.delv = sample.submelv - sample.elv;
    end;
    
    % spallation production scaling
    Psp = P_sp_expage(sample.pressure,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,nucl10,nucl26);
    
    % interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
    Lsp = rawattenuationlength(sample.pressure,LSDfix.Rc);
    
    % fix parameters for calculation preparations
    pars.thick = sample.thick;
    pars.Lsp = Lsp(:);
    pars.shield = sample.shield;
    pars.samplingyr = sample.samplingyr;
    pars.iceproxy = iceproxy(:);
    pars.tv = sample.tv;
    pars.z_mu = z_mu;
    if exist('fixice'); pars.fixice = fixice; end;
    if exist('fixE'); pars.fixE = fixE; end;
    if isfield(sample,'burialdepth'); pars.burialdepth = sample.burialdepth; end;
    
    % fix variable parameters
    pars.glacE = sample.glacE;
    pars.icevalue = sample.icevalue;
    pars.erosion = sample.erosion;
    pars.dens = sample.dens;
    if isfield(sample,'deglac'); pars.deglac = sample.deglac; end;
    if isfield(sample,'delv'); pars.delv = sample.delv; pars.elv = sample.elv; end;
    
    % calculate ice cover, erosion, and depth matrices
    prepar = precalc(pars,glacErate,glacEstep);
    
    % nuclide-specific erosion calculations
    if nucl10 == 1; % 10Be
        % fix nuclide specific parameters
        np10.Psp = Psp.sp10;
        np10.Pmu_z = Pmu_z.mu10;
        np10.N = sample.N10; np10.Nunc = sample.N10unc;
        np10.Pref = Pref10; np10.Prefunc = Pref10unc;
        np10.l = l10; np10.lunc = l10unc;
        % calculate erosion and internal unc
        N10E = nuclE(pars,prepar,np10,Etestv,mc,glaccalc);
    end;
    if nucl26 == 1; % 26Al
        % fix nuclide specific parameters
        np26.Psp = Psp.sp26;
        np26.Pmu_z = Pmu_z.mu26;
        np26.N = sample.N26; np26.Nunc = sample.N26unc;
        np26.Pref = Pref26; np26.Prefunc = Pref26unc;
        np26.l = l26; np26.lunc = l26unc;
        % calculate erosion and internal unc
        N26E = nuclE(pars,prepar,np26,Etestv,mc,glaccalc);
    end;
    
    % generate uncertainty parameters for mc iterations
    parsu = uncrnd(pars,sample,mc);
    
    % nuclide-specific uncertainty calculations
    if nucl10 == 1; % 10Be
        % estimate uncertainty
        unc10 = Euncest(N10E,parsu,np10,sample,mc,glaccalc);
        % display output and fill for saving
        if isfield(N10E,'r');
            fprintf(1,'\n10Be: %.2f +/- [%.2f / %.2f] mm/ka',N10E.r,unc10.rpos,unc10.rneg);
            output(i+1,outcol.E10r) = {num2str(N10E.r,'%.2f')};
            output(i+1,outcol.E10rp) = {num2str(unc10.rpos,'%.2f')};
            output(i+1,outcol.E10rn) = {num2str(unc10.rneg,'%.2f')};
            % fix for plotting and depth at time points
            if numel(tdv) > 0 || sum([plch.depth plch.P plch.N]) > 0;
                dPN = depthPconc(pars,parsu,N10E.r,N10E.rpos_int,N10E.rneg_int,unc10.rpos,...
                    unc10.rneg,unc10.Evr,Etestv,np10,plch,tdv,glaccalc,1,0);
                % fill output
                if isfield(dPN,'tdstr');
                    output(i+1,outcol.E10rd1:outcol.E10rd2) = dPN.tdstr;
                end;
                % fix for plotting
                pl10r = fix_plstruct(sample.tv,dPN,pl10r);
            end;
        end
        if isfield(N10E,'i');
            fprintf(1,'\n10Be: %.2f +/- [%.2f / %.2f]',N10E.i,unc10.ipos,unc10.ineg);
            if glaccalc == 1; fprintf(1,' cm/glac'); else; fprintf(1,' mm/ka'); end;
            output(i+1,outcol.E10i) = {num2str(N10E.i,'%.2f')};
            output(i+1,outcol.E10ip) = {num2str(unc10.ipos,'%.2f')};
            output(i+1,outcol.E10in) = {num2str(unc10.ineg,'%.2f')};
            % fix for plotting and depth at time points
            if numel(tdv) > 0 || sum([plch.depth plch.P plch.N]) > 0;
                dPN = depthPconc(pars,parsu,N10E.i,N10E.ipos_int,N10E.ineg_int,unc10.ipos,...
                    unc10.ineg,unc10.Evi,Etestv,np10,plch,tdv,glaccalc,0,1);
                % fill output
                if isfield(dPN,'tdstr');
                    output(i+1,outcol.E10id1:outcol.E10id2) = dPN.tdstr;
                end;
                % fix for plotting
                pl10i = fix_plstruct(sample.tv,dPN,pl10i);
            end;
        end
    end;
    if nucl26 == 1; % 26Al
        % estimate uncertainty
        unc26 = Euncest(N26E,parsu,np26,sample,mc,glaccalc);
        % display output and fill for saving
        if isfield(N26E,'r');
            fprintf(1,'\n26Al: %.2f +/- [%.2f / %.2f] mm/ka',N26E.r,unc26.rpos,unc26.rneg);
            output(i+1,outcol.E26r) = {num2str(N26E.r,'%.2f')};
            output(i+1,outcol.E26rp) = {num2str(unc26.rpos,'%.2f')};
            output(i+1,outcol.E26rn) = {num2str(unc26.rneg,'%.2f')};
            % fix for plotting and depth at time points
            if numel(tdv) > 0 || sum([plch.depth plch.P plch.N]) > 0;
                dPN = depthPconc(pars,parsu,N26E.r,N26E.rpos_int,N26E.rneg_int,unc26.rpos,...
                    unc26.rneg,unc26.Evr,Etestv,np26,plch,tdv,glaccalc,1,0);
                % fill output
                if isfield(dPN,'tdstr');
                    output(i+1,outcol.E26rd1:outcol.E26rd2) = dPN.tdstr;
                end;
                % fix for plotting
                pl26r = fix_plstruct(sample.tv,dPN,pl26r);
            end;
        end
        if isfield(N26E,'i');
            fprintf(1,'\n26Al: %.2f +/- [%.2f / %.2f]',N26E.i,unc26.ipos,unc26.ineg);
            if glaccalc == 1; fprintf(1,' cm/glac'); else; fprintf(1,' mm/ka'); end;
            output(i+1,outcol.E26i) = {num2str(N26E.i,'%.2f')};
            output(i+1,outcol.E26ip) = {num2str(unc26.ipos,'%.2f')};
            output(i+1,outcol.E26in) = {num2str(unc26.ineg,'%.2f')};
            % fix for plotting and depth at time points
            if numel(tdv) > 0 || sum([plch.depth plch.P plch.N]) > 0;
                dPN = depthPconc(pars,parsu,N26E.i,N26E.ipos_int,N26E.ineg_int,unc26.ipos,...
                    unc26.ineg,unc26.Evi,Etestv,np26,plch,tdv,glaccalc,0,1);
                % fill output
                if isfield(dPN,'tdstr');
                    output(i+1,outcol.E26id1:outcol.E26id2) = dPN.tdstr;
                end;
                % fix for plotting
                pl26i = fix_plstruct(sample.tv,dPN,pl26i);
            end;
        end
    end;
    % 10Be + 26Al
    np1026.nofield = 0; % fix for dPN
    if nucl10==1 && nucl26==1 && nucl1026==1 && isfield(N10E,'Nendr') && isfield(N26E,'Nendr');
        % calculate best fit erosion
        E1026r = Ecalc1026(N10E.r,N10E.rpos_int,N10E.rneg_int,N26E.r,N26E.rpos_int,N26E.rneg_int,...
            N10E.Nendr,N26E.Nendr,unc10.Evr,unc26.Evr,N10E.maxEr,N26E.maxEr,sample,Etestv);
        % display output and fill for saving
        if isfield(E1026r,'Epos');
            fprintf(1,'\n10Be+26Al: %.2f +/- [%.2f / %.2f] mm/ka',E1026r.E,E1026r.Epos,E1026r.Eneg);
            output(i+1,outcol.E1026r) = {num2str(E1026r.E,'%.2f')};
            output(i+1,outcol.E1026rp) = {num2str(E1026r.Epos,'%.2f')};
            output(i+1,outcol.E1026rn) = {num2str(E1026r.Eneg,'%.2f')};
            % fix for plotting and depth at time points
            if numel(tdv)>0 || plch.depth==1;
                dPN = depthPconc(pars,parsu,E1026r.E,E1026r.intpos,E1026r.intneg,E1026r.Epos,...
                    E1026r.Eneg,E1026r.Ev,Etestv,np1026,plch,tdv,glaccalc,1,0);
                % fill output
                if isfield(dPN,'tdstr');
                    output(i+1,outcol.E1026rd1:outcol.E1026rd2) = dPN.tdstr;
                end;
                % fix for plotting
                pl1026r = fix_plstruct(sample.tv,dPN,pl1026r);
            end;
        else;
            fprintf(1,'\n10Be+26Al: no overlap');
            output(i+1,outcol.E1026r) = {'nohit'};
            output(i+1,outcol.E1026rp) = {'nohit'};
            output(i+1,outcol.E1026rn) = {'nohit'};
        end;
    end;
    if nucl10==1 && nucl26==1 && nucl1026==1 && isfield(N10E,'Nendi') && isfield(N26E,'Nendi');
        % calculate best fit erosion
        E1026i = Ecalc1026(N10E.i,N10E.ipos_int,N10E.ineg_int,N26E.i,N26E.ipos_int,N26E.ineg_int,...
            N10E.Nendi,N26E.Nendi,unc10.Evi,unc26.Evi,N10E.maxEi,N26E.maxEi,sample,Etestv);
        % display output and fill for saving
        if isfield(E1026i,'Epos');
            fprintf(1,'\n10Be+26Al: %.2f +/- [%.2f / %.2f]',E1026i.E,E1026i.Epos,E1026i.Eneg);
            if glaccalc == 1; fprintf(1,' cm/glac'); else; fprintf(1,' mm/ka'); end;
            output(i+1,outcol.E1026i) = {num2str(E1026i.E,'%.2f')};
            output(i+1,outcol.E1026ip) = {num2str(E1026i.Epos,'%.2f')};
            output(i+1,outcol.E1026in) = {num2str(E1026i.Eneg,'%.2f')};
            % fix for plotting and depth at time points
            if numel(tdv)>0 || plch.depth==1;
                dPN = depthPconc(pars,parsu,E1026i.E,E1026i.intpos,E1026i.intneg,E1026i.Epos,...
                    E1026i.Eneg,E1026i.Ev,Etestv,np1026,plch,tdv,glaccalc,0,1);
                % fill output
                if isfield(dPN,'tdstr');
                    output(i+1,outcol.E1026id1:outcol.E1026id2) = dPN.tdstr;
                end;
                % fix for plotting
                pl1026i = fix_plstruct(sample.tv,dPN,pl1026i);
            end;
        else;
            fprintf(1,'\n10Be+26Al: no overlap');
            output(i+1,outcol.E1026i) = {'nohit'};
            output(i+1,outcol.E1026ip) = {'nohit'};
            output(i+1,outcol.E1026in) = {'nohit'};
        end;
    end;
    
    % new line
    fprintf(1,'\n');
    
    % clear parameters
    clear sample; clear pars; clear np10; clear np26;
    clear N10E; clear N26E;
    clear unc10; clear unc26;
    clear E1026r; clear E1026i;
end;

% plotting ======================================================================================
% do plotting
if plnum.d10r > 0;
    pl10r = do_plot(pl10r,plnum.d10r,plch.uncline,plch.clr10,'^{10}Be');
end;
if plnum.d26r > 0;
    pl26r = do_plot(pl26r,plnum.d26r,plch.uncline,plch.clr26,'^{26}Al');
end;
if plnum.d1026r>0 && numel(pl1026r.tm)>0;
    pl1026r = do_plot(pl1026r,plnum.d1026r,plch.uncline,plch.clr1026,'^{10}Be + ^{26}Al');
end;
if plnum.d10i > 0;
    pl10i = do_plot(pl10i,plnum.d10i,plch.uncline,plch.clr10,'^{10}Be');
end;
if plnum.d26i > 0;
    pl26i = do_plot(pl26i,plnum.d26i,plch.uncline,plch.clr26,'^{26}Al');
end;
if plnum.d1026i>0 && numel(pl1026i.tm)>0;
    pl1026i = do_plot(pl1026i,plnum.d1026i,plch.uncline,plch.clr1026,'^{10}Be + ^{26}Al');
end;
% fix plots
if plnum.d10r > 0;
    if plnum.d10r == plnum.d26r;
        pl10r.leg(end+1) = pl26r.leg; pl10r.legin(end+1) = pl26r.legin;
    end;
    if plnum.d10r==plnum.d1026r && numel(pl1026r.tm)>0;
        pl10r.leg(end+1) = pl1026r.leg; pl10r.legin(end+1) = pl1026r.legin;
    end;
    lastfix_plot(pl10r,plnum.d10r,[0 plch.maxt 0 plch.maxd],1,'Time (ka)','Sample depth (m)');
end;
if plnum.d26r > plnum.d10r;
    lastfix_plot(pl26r,plnum.d26r,[0 plch.maxt 0 plch.maxd],1,'Time (ka)','Sample depth (m)');
end;
if plnum.d1026r>plnum.d10r && numel(pl1026r.tm)>0;
    lastfix_plot(pl1026r,plnum.d1026r,[0 plch.maxt 0 plch.maxd],1,'Time (ka)','Sample depth (m)');
end;
if plnum.d10i > 0;
    if plnum.d10i == plnum.d26i;
        pl10i.leg(end+1) = pl26i.leg; pl10i.legin(end+1) = pl26i.legin;
    end;
    if plnum.d10i==plnum.d1026i && numel(pl1026i.tm)>0;
        pl10i.leg(end+1) = pl1026i.leg; pl10i.legin(end+1) = pl1026i.legin;
    end;
    lastfix_plot(pl10i,plnum.d10i,[0 plch.maxt 0 plch.maxd],1,'Time (ka)','Sample depth (m)');
end;
if plnum.d26i > plnum.d10i;
    lastfix_plot(pl26i,plnum.d26i,[0 plch.maxt 0 plch.maxd],1,'Time (ka)','Sample depth (m)');
end;
if plnum.d1026i>plnum.d10i && numel(pl1026i.tm)>0;
    lastfix_plot(pl1026i,plnum.d1026i,[0 plch.maxt 0 plch.maxd],1,'Time (ka)','Sample depth (m)');
end;
% ===============================================================================================

% fix and save output ======================================
if sum(samplein.N10 + samplein.N26)>0;
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
    out = fopen('out-glacialE.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    fclose(out);
end;
% ==========================================================

toc()
clear;
% end glacialE function ============================================================================


% subfunction precalc ==============================================================================
function out = precalc(pars,glacErate,glacEstep);
% calculate ice cover, erosion, and depth matrices
    % find number of columns to use in matrices
    numcols = max(size(pars.glacE,2),size(pars.erosion,2));
    
    % fix ice cover record
    icem = (pars.iceproxy >= pars.icevalue);
    icem = double(icem); % fix for matlab
    
    % fix deglaciation based on sample.deglac
    if isfield(pars,'deglac');
        icem = fix_deglac(pars.deglac,icem,pars.tv);
    end;
    
    % fix specified ice-cover/noice-cover periods
    if isfield(pars,'fixice');
        icem = fixicef(pars.fixice,icem,pars.samplingyr,pars.tv);
    end;
    
    % fix icetest matrix
    icetest = icem;
    
    % save icev for plotting and fix ice matrix
    if size(icem,2) == 1;
        out.icev = icem;
        icem = repmat(icem(:),1,numcols);
    end;
    
    % fix noice matrix and save to output
    noicem = (icem == 0);
    out.noicem = noicem;
    
    % fix non-glacial erosion matrix
    nonglacEm = repmat(pars.erosion,numel(pars.tv),numcols./numel(pars.erosion));
    
    % fix glacial erosion matrices
    if glacErate == 1;
        glacEm.rr = repmat(pars.glacE,numel(pars.tv),numcols./numel(pars.glacE));
        out.icetestr = icetest;
    end;
    if glacEstep == 1;
        glacEm.ii = repmat(pars.glacE,numel(pars.tv),numcols./numel(pars.glacE));
        out.icetesti = icetest;
    end;
    
    % fix non-glacial and glacial erosion in specific time periods
    if isfield(pars,'fixE');
        [nonglacEm,glacEm,out] = fixEf(pars.fixE,nonglacEm,glacEm,pars.samplingyr,pars.tv,out);
    end;
    
    % fix dens matrix and save to output
    densm = repmat(pars.dens,numel(pars.tv),numcols./numel(pars.dens)); % fix dens matrix
    out.densm = densm;
    
    % fix water depth matrix - water is assumed to have a density of 1 g/cm3
    waterdm = zeros(size(icem));
    if isfield(pars,'delv');
        [waterdm noicem] = get_waterdepth(pars.elv,pars.delv,icem,noicem,pars.tv);
    end;
    
    % calculate non-glacial erosion depth
    nonglacdm = cumtrapz(pars.tv,noicem.*nonglacEm.*1E-4).*densm; % g/cm2
    if isfield(pars,'burialdepth');
        nonglacdm = nonglacdm + pars.burialdepth.*densm; % add burial depth (g/cm2)
    end;
    
    % calculate glacial erosion depth matrix
    if isfield(glacEm,'rr');
        glacdrm = cumtrapz(pars.tv,icem.*glacEm.rr.*1E-4).*densm; % g/cm2
    end;
    if isfield(glacEm,'ri');
        % find deglaciation events
        deglacm = (filter([1 2],1,icem) == 1);
        glacdrm = glacdrm + cumsum(deglacm.*glacEm.ri).*densm; % g/cm2
    end;
    if isfield(glacEm,'ii');
        % find deglaciation events
        deglacm = (filter([1 2],1,icem) == 1);
        glacdim = cumsum(deglacm.*glacEm.ii).*densm; % g/cm2
    end;
    if isfield(glacEm,'ir');
        glacdim = glacdim + cumtrapz(pars.tv,icem.*glacEm.ir.*1E-4).*densm; % g/cm2
    end;
    
    % calculate full depth matrix (glac and nonglac plus water depth) and find deglac index
    if exist('glacdrm');
        out.fulldmr = glacdrm + nonglacdm + waterdm;
        out.dmmr = (out.fulldmr-waterdm)./densm.*1E-2; % depth matrix (m)
    end;
    if exist('glacdim');
        out.fulldmi = glacdim + nonglacdm + waterdm;
        out.dmmi = (out.fulldmi-waterdm)./densm.*1E-2; % depth matrix (m)
    end;
% end subfunction precalc ==========================================================================


% subfunction fix_deglac ===========================================================================
function out = fix_deglac(deglac,icem,tv);
    % fix deglaciation in icem
    out = icem;
    if numel(deglac) > 1;
        tv = tv(:);
        deglac = repmat(deglac(:)',size(tv));
        tv = repmat(tv,1,size(icem,2));
    end;
    predeglac = (tv>=deglac);
    predeglacin = (cumsum(icem)>=1);
    newice = (predeglac-predeglacin == 1);
    out = icem.*predeglac + newice;
% end subfunction fix_deglac =======================================================================


% subfunction fixicef ==============================================================================
function icem = fixicef(fixice,icem,samplingyr,tv);
    % fix specified period(s) with ice cover
    if isfield(fixice,'ice');
        ice = fixice.ice + samplingyr - fixice.fix0yr;
        for j = 1:size(ice,1);
            fixmin = find(tv >= ice(j,1),1,'first');
            fixmax = find(tv <= ice(j,2),1,'last');
            icem(fixmin:fixmax,:) = 1;
        end;
    end;
    % fix specified period(s) with no ice cover
    if isfield(fixice,'noice');
        noice = fixice.noice + samplingyr - fixice.fix0yr;
        for j = 1:size(noice,1);
            fixmin = find(tv >= noice(j,1),1,'first');
            fixmax = find(tv <= noice(j,2),1,'last');
            icem(fixmin:fixmax,:) = 0;
        end;
    end;
% end subfunction fixicef ==========================================================================


% subfunction fixEf ================================================================================
function [nonglacEm,glacEm,out] = fixEf(fixE,nonglacEm,glacEm,samplingyr,tv,out);
    % fix specific nonglacial erosion in specified period(s)
    if isfield(fixE,'nonglacE') && isfield(fixE,'nonglacE_tv');
        nonglacE_tv = fixE.nonglacE_tv + samplingyr - fixE.fix0yr;
        for j = 1:size(nonglacE_tv,1);
            fixmin = find(tv >= nonglacE_tv(j,1),1,'first');
            fixmax = find(tv <= nonglacE_tv(j,2),1,'last');
            nonglacEm(fixmin:fixmax,:) = fixE.nonglacE(j);
        end;
    end;
    % fix specific glacial erosion rate in specified period(s)
    if isfield(fixE,'glacEr') && isfield(fixE,'glacEr_tv');
        if isfield(glacEm,'ii'); glacEm.ir = zeros(size(glacEm.ii)); end;
        glacEr_tv = fixE.glacEr_tv + samplingyr - fixE.fix0yr;
        for j = 1:size(glacEr_tv,1);
            fixmin = min(find(tv >= glacEr_tv(j,1)));
            fixmax = max(find(tv <= glacEr_tv(j,2)));
            if isfield(glacEm,'rr'); glacEm.rr(fixmin:fixmax,:) = fixE.fix_glacEr(j); end;
            if isfield(glacEm,'ii'); glacEm.ir(fixmin:fixmax,:) = fixE.fix_glacEr(j); end;
            out.icetestr(fixmin:fixmax,:) = 0;
        end;
    end;
    % fix specific incremental glacial erosion step in specified period(s)
    if isfield(fixE,'glacEi') && isfield(fixE,'glacEi_tv');
        if isfield(glacEm,'rr'); glacEm.ri = zeros(size(glacEm.rr)); end;
        glacEi_tv = fixE.glacEi_tv + samplingyr - fix0yr;
        for j = 1:size(glacEi_tv,1);
            fixmin = min(find(tv >= glacEi_tv(j,1)));
            fixmax = max(find(tv <= glacEi_tv(j,2)));
            if isfield(glacEm,'ii'); glacEm.ii(fixmin:fixmax,:) = fixE.fix_glacEi(j); end;
            if isfield(glacEm,'rr'); glacEm.ri(fixmin:fixmax,:) = fixE.fix_glacEi(j); end;
            out.icetesti(fixmin:fixmax,:) = 0;
        end;
    end;
% end subfunction fixEf ============================================================================


% subfunction nuclE ================================================================================
function out = nuclE(pars,prepar,np,Evect,mc,glaccalc);
    % calculate end conc and interpolate erosion
    Nmc = normrnd(np.N,np.Nunc,[1 mc]); % N vect for N uncertainty estimation
    if isfield(prepar,'fulldmr'); % glacial erosion rate
        out.Nendr = Ncalc(pars,prepar,np.Psp,np.Pmu_z,np.Pref,np.l,prepar.fulldmr);
        out.r = Ecalc(out.Nendr,Evect,np.N);
        Emcr = Ecalc(out.Nendr,Evect,Nmc);
        Emcrp = Emcr(Emcr>=out.r);
        Emcrn = Emcr(Emcr<=out.r);
        out.rpos_int = std([Emcrp 2*out.r-Emcrp]);
        out.rneg_int = std([Emcrn 2*out.r-Emcrn]);
        if glaccalc == 1;
            Elimr = minmaxE(pars,prepar,np,prepar.fulldmr,prepar.icetestr,out.Nendr,Evect);
            out.maxEr = Elimr.maxE;
            if out.r > Elimr.maxE; out.r = Elimr.maxE; end;
            if isfield(Elimr,'Epos'); out.rpos = Elimr.Epos; end;
            if isfield(Elimr,'Eneg'); out.rneg = Elimr.Eneg; end;
        else;
            out.maxEr = max(Evect);
        end;
    end;
    if isfield(prepar,'fulldmi'); % incremental glacial erosion steps
        out.Nendi = Ncalc(pars,prepar,np.Psp,np.Pmu_z,np.Pref,np.l,prepar.fulldmi);
        out.i = Ecalc(out.Nendi,Evect,np.N);
        Emci = Ecalc(out.Nendi,Evect,Nmc);
        Emcip = Emci(Emci>=out.i);
        Emcin = Emci(Emci<=out.i);
        out.ipos_int = std([Emcip 2*out.i-Emcip]);
        out.ineg_int = std([Emcin 2*out.i-Emcin]);
        if glaccalc == 1;
            Elimi = minmaxE(pars,prepar,np,prepar.fulldmi,prepar.icetesti,out.Nendi,Evect);
            out.maxEi = Elimi.maxE;
            if out.i > Elimi.maxE; out.i = Elimi.maxE; end;
            if isfield(Elimi,'Epos'); out.ipos = Elimi.Epos; end;
            if isfield(Elimi,'Eneg'); out.ineg = Elimi.Eneg; end;
        else;
            out.maxEi = max(Evect);
        end;
    end;
    out.Evect = Evect;
% end subfunction nuclE ============================================================================


% subfunction Ncalc ================================================================================
function Nend = Ncalc(pars,prepar,Psp,Pmu_z,Pref,l,depthm);
    % fix thickness scaling factor
    Lspm = repmat(pars.Lsp,1,size(prepar.noicem,2)); % fix Lsp matrix
    if pars.thick > 0;
        thickSF = (Lspm./(prepar.densm.*pars.thick)) .* ...
            (1 - exp(((-1.*prepar.densm.*pars.thick)./Lspm)));
    else;
        thickSF = 1;
    end;
    % Calculate N(end) including decay and erosion
    Pspm = repmat(Psp(:),1,size(prepar.noicem,2));
    Prefm = repmat(Pref,numel(pars.tv),size(Pspm,2)./numel(Pref));
    Pspm = Pspm.*Prefm.*thickSF.*pars.shield.*prepar.noicem; % surface spal prod matrix
    dpfs = exp(-depthm./Lspm); % spal depth dependence matrix
    Pmum = interp1(pars.z_mu',Pmu_z',depthm+pars.thick.*prepar.densm./2,'pchip') .* ...
        pars.shield .* prepar.noicem; % muon prod matrix
    Pmum(isnan(Pmum)) = 0; % set muon prod to 0 for depths > 2E5
    lm = repmat(l,numel(pars.tv),size(Pspm,2)./numel(l)); % lambda matrix
    tm = repmat(pars.tv,1,size(Pspm,2)); % time matrix
    dcf = exp(-tm.*lm); % decay factor
    Nend = trapz(pars.tv,(Pspm.*dcf.*dpfs + Pmum.*dcf)); % N from full tv
% end subfunction Ncalc ============================================================================


% subfunction Ecalc ================================================================================
function NE = Ecalc(Nend,Evect,N);
    % remove points with same N from Nend and Evect
    rmidx = find(diff(Nend)==0,1,'first');
    if numel(rmidx)>0; Nend(rmidx+1:end) = []; Evect(rmidx+1:end) = []; end;
    % set too low and too high N to limits of Nend
    N(N<min(Nend)) = min(Nend);
    N(N>max(Nend)) = max(Nend);
    % interpolate erosion
    NE = interp1(Nend,Evect,N,'pchip');
% end subfunction Ecalc ============================================================================


% subfunction minmaxE ==============================================================================
function out = minmaxE(pars,prepar,np,depthm,icetest,Nend,Evect);
    out.maxE = max(Evect);
    % if saturated
    if np.N > max(Nend);
        out.Eneg = 0;
        if np.N-np.Nunc > max(Nend);
            out.Epos = 0;
        else;
            out.Epos = Ecalc(Nend,Evect,np.N-np.Nunc);
        end;
    end;
    % fix maxE
    if sum(icetest) > 0 && icetest(1) == 0;
        % find index of last ice cover period and use to cut vectors and matrices
        deglidx = find(icetest==1,1,'first');
        pars.Lsp = pars.Lsp(1:deglidx);
        prepar.noicem = prepar.noicem(1:deglidx,:);
        prepar.densm = prepar.densm(1:deglidx,:);
        np.Psp = np.Psp(1:deglidx);
        pars.tv = pars.tv(1:deglidx);
        depthm = depthm(1:deglidx,:);
        Npostglac = Ncalc(pars,prepar,np.Psp,np.Pmu_z,np.Pref,np.l,depthm);
        Ndeglac = Nend - Npostglac;
        % set max erosion to 1% of surface conc expected from postglacial exposure
        out.maxE = Ecalc(Ndeglac,Evect,Npostglac(1).*0.01);
        % if no inheritance
        if np.N < Npostglac(1);
            out.Epos = max(Evect);
            if np.N+np.Nunc > Npostglac(1).*1.01;
                out.Eneg = out.maxE-Ecalc(Ndeglac,Evect,np.N+np.Nunc-Npostglac(1));
            else;
                out.Eneg = 0;
            end;
        end;
    end;
% end subfunction minmaxE ==========================================================================


% subfunction uncrnd ===============================================================================
function pars = uncrnd(pars,sample,mc);
    % fix parameters with uncertainties
    if isfield(sample,'densunc');
        pars.dens = normrnd(sample.dens,sample.densunc,[1 mc]);
        pars.dens(pars.dens<0) = 0;
    end;
    if isfield(sample,'erosionunc');
        pars.erosion = normrnd(sample.erosion,sample.erosionunc,[1 mc]);
        pars.erosion(pars.erosion<0) = 0;
    else;
        pars.erosion = repmat(pars.erosion,1,mc);
    end;
    if isfield(sample,'glacEunc');
        pars.glacE = normrnd(sample.glacE,sample.glacEunc,[1 mc]);
        pars.glacE(pars.glacE<0) = 0;
    else;
        pars.glacE = repmat(pars.glacE,1,mc);
    end;
    if isfield(sample,'deglacunc');
        pars.deglac = normrnd(sample.deglac,sample.deglacunc,[1 mc]);
    end;
    if isfield(sample,'icevalueunc');
        pars.icevalue = normrnd(sample.icevalue,sample.icevalueunc,[1 mc]);
    end;
    if isfield(sample,'isostsubmunc');
        if strcmp(sample.isostsubm,'-')==0 && isnan(sample.isostsubmunc)==0;
            submuncv = normrnd(1,sample.isostsubmunc,[1 mc]);
            pars.delv = repmat(pars.delv,1,mc) .* repmat(submuncv,size(pars.delv));
        end;
    end;
% end subfunction uncrnd ===========================================================================


% subfunction Euncest ==============================================================================
function out = Euncest(NE,parsu,np,sample,mc,glaccalc);
    % fix Pref and l for mc iterations
    uncpar.Pref = normrnd(np.Pref,np.Prefunc,[1 mc]);
    uncpar.l = normrnd(np.l,np.lunc,[1 mc]);
    
    % calculate uncertainty for glacial erosion rate case
    if isfield(NE,'r');
        % fix erosion parameter to test
        if glaccalc == 1; parsu.glacE = NE.r; else; parsu.erosion = NE.r; end;
        
        % calculate ice cover, erosion, and depth matrices
        prepar = precalc(parsu,1,0);
        
        % fix specific unc parameters
        uncpar.fulldm = prepar.fulldmr;
        uncpar.Nend = NE.Nendr;
        uncpar.Emid = NE.r;
        uncpar.intpos = NE.rpos_int;
        uncpar.intneg = NE.rneg_int;
        uncpar.maxE = NE.maxEr;
        
        % do calculation
        [out.rpos out.rneg out.Evr] = Eunccalc(NE,parsu,prepar,np,uncpar);
        
        % check for unc limits from nuclE
        if isfield(NE,'rpos'); out.rpos = NE.rpos; end;
        if isfield(NE,'rneg'); out.rneg = NE.rneg; end;
        
        % save parsu for plotting and specific time depth
        out.parsr = parsu;
    end;
    
    % calculate uncertainty for incremental depth step case
    if isfield(NE,'i');
        % fix erosion parameter to test
        if glaccalc == 1; parsu.glacE = NE.i; else; parsu.erosion = NE.i; end;
        
        % calculate ice cover, erosion, and depth matrices
        prepar = precalc(parsu,0,1);
        
        % fix specific unc parameters
        uncpar.fulldm = prepar.fulldmi;
        uncpar.Nend = NE.Nendi;
        uncpar.Emid = NE.i;
        uncpar.intpos = NE.ipos_int;
        uncpar.intneg = NE.ineg_int;
        uncpar.maxE = NE.maxEi;
        
        % do calculation
        [out.ipos out.ineg out.Evi] = Eunccalc(NE,parsu,prepar,np,uncpar);
        
        % check for unc limits from nuclE
        if isfield(NE,'ipos'); out.ipos = NE.ipos; end;
        if isfield(NE,'ineg'); out.ineg = NE.ineg; end;
        
        % save parsu for plotting and specific time depth
        out.parsi = parsu;
    end;
% end subfunction Euncest ==========================================================================


% subfunction Eunccalc =============================================================================
function [uncp uncn Ev] = Eunccalc(NE,parsu,prepar,np,uncpar);
    % calculate Nend
    Nend = Ncalc(parsu,prepar,np.Psp,np.Pmu_z,uncpar.Pref,uncpar.l,uncpar.fulldm);
    
    % calculate E based on the central point Nend and Evect
    Ev = Ecalc(uncpar.Nend,NE.Evect,Nend);
    
    % split in positive and negative E and calculate deviation
    Evp = Ev(Ev>=uncpar.Emid);
    Evn = Ev(Ev<=uncpar.Emid);
    extpos = std([Evp 2*uncpar.Emid-Evp]);
    extneg = std([Evn 2*uncpar.Emid-Evn]);
    
    % calculate pos and neg deviation of Nend and check for deviation outside NE.Nend
    Nendp = Nend(Nend>=np.N);
    Nendn = Nend(Nend<=np.N);
    Nendpos = std([Nendp 2*np.N-Nendp]);
    Nendneg = std([Nendn 2*np.N-Nendn]);
    if np.N-Nendneg < min(uncpar.Nend); extpos = max(NE.Evect); end;
    if np.N+Nendpos > max(uncpar.Nend); extneg = uncpar.Emid; end;
    
    % combine internal and external uncertainties
    uncp = sqrt(uncpar.intpos.^2 + extpos.^2);
    uncn = sqrt(uncpar.intneg.^2 + extneg.^2);
    
    % don't allow too large positive unc
    if uncpar.Emid+uncp > uncpar.maxE;
        uncp = max(NE.Evect);
    end;
% end subfunction Eunccalc =========================================================================


% subfunction Ecalc1026 ============================================================================
function out = Ecalc1026(E10,E10p,E10n,E26,E26p,E26n,Nend10,Nend26,Ev10,Ev26,maxE10,maxE26,...
    sample,Evect);
    % calculate weighted E using internal uncertainties and the EVM method
    [out.E,out.intpos,out.intneg] = evm([E10 E26],[E10p E26p],[E10n E26n]);
    
    % calculate theoretical N at combined E
    NE10 = interp1(Evect,Nend10,out.E,'pchip');
    NE26 = interp1(Evect,Nend26,out.E,'pchip');
    
    % check overlap and estimate positive and negative uncertainty
    if abs(NE10-sample.N10)<=sample.N10unc && abs(NE26-sample.N26)<=sample.N26unc;
        % calculate mean E for external uncertainty estimation - use arithmetic mean of 10/26 values
        out.Ev = mean([Ev10;Ev26]);
        
        % split in positive and negative E and calculate deviation
        Evp = out.Ev(out.Ev>=out.E);
        Evn = out.Ev(out.Ev<=out.E);
        extpos = std([Evp 2*out.E-Evp]);
        extneg = std([Evn 2*out.E-Evn]);
        
        % combine internal and external uncertainties
        out.Epos = sqrt(out.intpos.^2 + extpos.^2);
        out.Eneg = sqrt(out.intneg.^2 + extneg.^2);
        
        % don't allow too large positive unc
        if out.E+out.Epos > min(maxE10,maxE26);
            out.Epos = max(Evect);
        end;
    end;
% end subfunction Ecalc1026 ========================================================================


% subfunction depthPconc ===========================================================================
function out = depthPconc(pars,parsv,Emid,uncintp,uncintn,uncextp,uncextn,Ev,Etestv,nuclpar,plch,...
    tdv,glaccalc,glacErate,glacEstep);
    % brief description
    % fix mid-point erosion parameter for depth matrix
    if glaccalc == 1; pars.glacE = Emid; else; pars.erosion = Emid; end;
    % calculate depth for mid-point E
    if glacErate == 1;
        prepar = precalc(pars,1,0);
        dvm = prepar.dmmr;
        deglidx = find(prepar.icetestr==1,1,'first');
    end;
    if glacEstep == 1;
        prepar = precalc(pars,0,1);
        dvm = prepar.dmmi;
        deglidx = find(prepar.icetesti==1,1,'first');
    end;
    % add internal uncertainty to Ev
    Ev_addp = abs(normrnd(0,uncintp,1,numel(Ev)/2));
    Ev_addn = -abs(normrnd(0,uncintn,1,numel(Ev)/2));
    Ev_add = [Ev_addp Ev_addn];
    Ev = Ev + Ev_add;
    % fix erosion parameter for uncertainty depth matrix
    if glaccalc == 1; parsv.glacE = Ev; else; parsv.erosion = Ev; end;
    % calculate depth matrix
    if glacErate == 1;
        preparv = precalc(parsv,1,0);
        dmm = preparv.dmmr;
        deglidxunc = find(sum(preparv.icetestr')>0,1,'first');
    end;
    if glacEstep == 1;
        preparv = precalc(parsv,0,1);
        dmm = preparv.dmmi;
        deglidxunc = find(sum(preparv.icetesti')>0,1,'first');
    end;
    % plot depth or get depth at specific time?
    if plch.depth == 1 || numel(tdv) > 0;
        % calculate pos and neg uncertainties for depth history
        dvmm = repmat(dvm',size(dmm,2),1);
        idxp = (dmm' >= dvmm);
        idxn = (dmm' <= dvmm);
        dvmp = sqrt(sum((dmm'-dvmm).^2.*idxp).*2./(sum(idxp).*2-1));
        dvmn = sqrt(sum((dmm'-dvmm).^2.*idxn).*2./(sum(idxn).*2-1));
        % fix for unconstrained erosion with constant depth at 1 km
        if Emid == max(Etestv); dvm(deglidx:end) = 1E3; end;
        if uncextp == max(Etestv); dvmp(deglidxunc:end) = 1E3; end;
        % save depth and pos and neg uncertainty
        out.depth = dvm;
        out.depthp = cummax(dvm+dvmp'); % cummax to avoid wobbling uncertainty
        out.depthn = flip(cummin(flip(dvm-dvmn'))); % cummin to avoid wobbling uncertainty
        % fix time-depth strings for output
        for j = 1:numel(tdv);
            out.tdstr(j*3-2) = {num2str(interp1(pars.tv,out.depth,tdv(j)),'%.2f')};
            out.tdstr(j*3-1) = {num2str(interp1(pars.tv,out.depthp-out.depth,tdv(j)),'%.2f')};
            out.tdstr(j*3) = {num2str(interp1(pars.tv,out.depth-out.depthn,tdv(j)),'%.2f')};
        end;
    end;
% end subfunction depthPconc =======================================================================


% subfunction fix_plstruct =========================================================================
function out = fix_plstruct(tv,dPN,out);
    % fix plot structure
    num = numel(tv);
    out.tm(1:num,end+1) = tv.*1E-3;
    out.tmunc(1:num*2,end+1) = [flip(tv.*1E-3);tv.*1E-3];
    out.mm(1:num,end+1) = 1;
    out.mmunc(1:num*2,end+1) = 1;
    if isfield(dPN,'depth');
        out.dm(1:num,end+1) = dPN.depth;
        out.dmunc(1:num*2,end+1) = [flip(dPN.depthn);dPN.depthp];
    end;
% end subfunction fix_plstruct =====================================================================


% subfunction do_plot ==============================================================================
function pl = do_plot(pl,num,uncline,clr,nstr);
    figure(num); hold on; box on;
    % fix for samples with varying tv
    pl.tm(pl.mm==0) = NaN; pl.dm(pl.mm==0) = NaN;
    pl.tmunc(pl.mmunc==0) = NaN; pl.dmunc(pl.mmunc==0) = NaN;
    % plot mid line
    pl.leg = plot(pl.tm,pl.dm,'color',clr);
    pl.leg = pl.leg(1); % fix for legend
    pl.legin = {nstr};
    % plot uncertainties
    if uncline == 1;
        plot(pl.tmunc,pl.dmunc,'color',min([clr+0.7;1 1 1]));
    else;
        patch(pl.tmunc,pl.dmunc,clr,'EdgeColor','none','FaceAlpha',0.1);
    end;
    hold off;
% end subfunction do_plot ==========================================================================


% subfunction lastfix_plot =========================================================================
function lastfix_plot(pl,num,axv,yreverse,xlab,ylab);
    figure(num);
    axis(axv);
    if yreverse == 1; axis('ij'); end;
    set(gca(),'xdir','reverse');
    xlabel(xlab); ylabel(ylab);
    legend(pl.leg,pl.legin,'location',pl.legloc);
    set(gca,'layer','top'); % plot axis on top
% end subfunction lastfix_plot =====================================================================


% subfunction get_waterdepth =======================================================================
function [waterdm noicem] = get_waterdepth(elv,delv,icem,noicem,tv);
    % fix postglacial delv matrix
    delvm = repmat(delv,1,size(icem,2)/size(delv,2)) .* (cumsum(icem)==0);
    delvm_min = min(delvm); % max elevation change (assumed to be at deglaciation)
    
    % flipped ice/noice matrices and tv
    flipicem = flip(icem);
    flipnoicem = flip(noicem);
    flipdtv = flip([diff(tv);0]);
    % uplift matrix
    uplm(1,size(icem,2)) = 0;
    % uplift constant
    uplconst = 0.0002;
    % max uplift vector
    maxuplv = repmat(30000,1,size(icem,2));
    
    % calculate uplift through simulation time
    for i = 2:numel(flipdtv);
        uplm(i,:) = min([uplm(i-1,:)+flipicem(i,:).*flipdtv(i);maxuplv]) .* ...
            exp(-uplconst.*flipdtv(i).*flipnoicem(i,:));
    end;
    
    % flip back
    uplm = flip(uplm);
    
    % calculated deglac and present-day uplift
    uplm_degl = repmat(max(uplm.*(cumsum(icem)==0)),size(tv));
    uplm_present = repmat(uplm(1,:),size(tv));
    
    % scale uplm based on delvm_min to get simulated delv
    delvmsim = (uplm-uplm_present)./uplm_degl.*delvm_min.*(cumsum(icem)>0);
    
    % combine postglacial and pre-LGM uplifts and remove uplift in ice cover periods
    delvmfull = (delvm + delvmsim) .* noicem;
    
    % water depth matrix (cm or g/cm2 if assuming a water density of 1)
    waterdm = (-delvmfull>elv).*(-delvmfull-elv).*1E2;
    
    % remove water covered periods from noicem
    noicem = noicem .* (waterdm==0);
% end subfunction get_waterdepth ===================================================================


% subfunction evm ==================================================================================
function [mid,uncp,uncn] = evm(midv,uncpv,uncnv);
    % average and uncertainty estimation using the "expected value method" (Birch and Singh 2014)
    num = numel(midv); % number of input samples
    
    % fix data (make matrices)
    midm = repmat(midv,num,1);
    uncpm = repmat(uncpv,num,1);
    uncnm = repmat(uncnv,num,1);
    xm = repmat(midv',1,num);
    
    % calculate probability matrix
    Mui = sum(sqrt(2./(pi.*(uncpm+uncnm).^2)).*exp(-(xm-midm).^2./(2.*uncnm.^2)).*(xm<=midm) + ...
        sqrt(2./(pi.*(uncpm+uncnm).^2)).*exp(-(xm-midm).^2./(2.*uncpm.^2)).*(xm>midm)) ./ num;
    
    % calculate summed probability for all samples
    Muj = sum(Mui);
    
    % calculate sample weights
    wi = Mui./Muj;
    
    % calculate weighted mean mid
    mid = sum(wi.*midv);
    
    % uncertainty estimation
    uncintp = sqrt(sum(wi.^2.*uncpv.^2));
    uncintn = sqrt(sum(wi.^2.*uncnv.^2));
    uncext = sqrt(sum(wi.*(midv-mid).^2));
    uncp = max(uncintp,uncext);
    uncn = max(uncintn,uncext);
% end subfunction evm ==============================================================================
