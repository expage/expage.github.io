function glacialE();

% 10Be and 26Al glacial erosion calculator.
% Function for quantification/investigation of glacial erosion based on 10Be and/or 26Al conc.
% Long and complicated code that can do some interesting calculations and plots.
% Interpolates the best fitting glacial erosion to yield the 10Be and/or 26Al conc given as input
% with a max duration (mt) of exposure, production during ice-free periods (determined by
% glacialE_LR04.txt) and shielding from non-glacial (given in input) and glacial (output) erosion.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2018 (jakob.heyman@gu.se)

% NOTE! This version works in Octave but for compatibility with Matlab some modification is needed

clear;
close all;
tic();

% What version is this?
ver = '201902';

% =============== MAKE CHOICES HERE ================================================================
% max time (determines the start of the simulation - max 1E7)
mt = 2588000;

% use glacial erosion rate (1) or incremental erosion depth steps (0)?
glacErate = 0;

% which calculation(s) to do? 1 = yes, 0 = no
simpleEcalc = 0;  % simple glacial erosion rate calculation (single nuclides)
range1calc = 0; % calculate possible ranges for parameters (single nuclides)
range2calc = 1; % calculate possible ranges for parameters (10Be + 26Al)

% depth of burial when sampled - for surface samples this is 0
burialdepth = 0; % cm

% fixed parameters for specific time period(s) =======================
% use multiple rows for multiple periods - comment out if not using
%fix_tv_nonglacE = [0 200000]; % yr  [min max]
%fix_nonglacE = 4; % mm/ka - non-glacial erosion rate
%fix_tv_glacEr = [0 30000;100000 300000]; % yr  [min max]
%fix_glacEr = [0;10]; % mm/ka - glacial erosion rate
%fix_tv_glacEi = [10000 11000]; % yr  [min max]
%fix_glacEi = [80]; % cm/glac - incremental glacial erosion steps

% fixed ice-cover and icefree for specific time period(s)
% use multiple rows for multiple periods - comment out if not using
%fix_ice = [40000 70000]; % yr  [min max]
%fix_nonice = [10000 23000]; % yr  [min max]

% simpleEcalc parameters =============================================
% LR04 per mill break value for determining ice cover history
LR04ice = 4.5;
simpleEcalcunc = 1; % plot uncertainties? 1 = yes

% range parameters (range1calc and/or range2calc) ====================
% glacial erosion - unit depending on glacErate above
glacEmin = 0; % mm/ka or cm/glac - min value
glacEmax = 10000; % mm/ka or cm/glac - max value
glacEstlim = 0.01; % maximum step

% non-glacial erosion rate
nonglacEmin = 0; % mm/ka - min value
nonglacEmax = 5; % mm/ka - max value
nonglacEstlim = 0.01; % maximum step

% LR04 per mill break values for determining ice cover history
LR04min = 4.4; % per mil - min value
LR04max = 4.6; % per mil - max value
LR04stlim = 0.01; % maximum step

% set number of parameters in matrix
parn = 8;

% add production rate uncertainty? 1 = yes
Punc = 1;

% random run parameters
rand1lim = 5; % number of random runs for each individual parameter limit
randhitN = 1E3; % minimum number of random hits; if 0 no random runs will be carried out
randrunsz = 150; % size of randrun vectors
randnummax = 1E2; % maximum number of random runs

% plotting choices =====================================================
pl_simpleEcalc_depth = 1;   % plot sample depth for simple glacial erosion (simpleEcalc = 1)
pl_simpleEcalc_buildup = 0; % plot nuclide build-up for simple glacial erosion (simpleEcalc = 1)
pl_range_depth = 1;   % plot sample depth for range calculation(s) (rangeXnulc = 1)
pl_range_buildup = 1; % plot nuclide build-up for range calculation(s) (rangeXnulc = 1)
pl_range_banana = 1;  % plot 26/10 banana for range calculation(s) (range2nulc = 1)
ratiopathN = 10; % number of OK sample 26Al/10Be ratio paths to plot in banana

% max time for plotting (yr)
pl_maxt = 1E5;
if pl_maxt > mt; pl_maxt = mt; end; % don't allow too large maxt
% max depth for plotting (m)
maxmaxy = 10; % max depth within simulated output depth range
%absmaxy = 10; % absolute max depth - comment out to not use

% plot colors
color10 = 'red'; % color for 10Be
color26 = 'blue'; % color for 26Al
color1026 = 'black'; % color for combined 10Be+26Al
% ==================================================================================================

% count number of input columns in line 1 of input
inid = fopen('input.txt','r');
line1 = fgets(inid);
line1 = regexprep(line1,' +',' ');
numcols = numel(strfind(line1,' ')) + numel(strfind(line1,sprintf('\t',''))) + 1;
fclose(inid);

% read input file
if numcols == 16; % if no deglaciation in input
    [samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
        samplein.rho,samplein.othercorr,samplein.E,samplein.N10,samplein.delN10,samplein.be_stds,...
        samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr] = textread...
        ('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n','commentstyle','matlab');
elseif numcols == 17; % if deglaciation in input
    [samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
        samplein.rho,samplein.othercorr,samplein.E,samplein.N10,samplein.delN10,samplein.be_stds,...
        samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr,samplein.deglac] = ...
        textread('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n %n','commentstyle',...
        'matlab');
else; % wrong number of inputs
    fprintf(1,'wrong number of input columns! the input has to contain (deglac is optional):\n');
    fprintf(1,'[sample lat lon elv Pflag thickn shield eros N10 delN10 std10 N26 delN26 samplyr ');
    fprintf(1,'deglac]\n');
    return;
end;

% run and load expage constants
make_consts_expage;
load consts_expage;

% constants
Pref10 = consts.P10_ref_nu; delPref10 = consts.delP10_ref_nu;
Pref26 = consts.P26_ref_nu; delPref26 = consts.delP26_ref_nu;
% Decay constant
l10 = consts.l10; dell10 = consts.dell10;
l26 = consts.l26; dell26 = consts.dell26;

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

% load Lisiecki and Raymo (2005) d18O ice volume record (LR04_tv same as LSD_fix tv to 1E7)
[LR04_tv,d18Oin] = textread('glacialE_LR04.txt','%n %n','commentstyle','matlab');

% glacial erosion vectors
if simpleEcalc==1 && glacErate==1; % if using erosion rate;
    Eglacv = [(0:2:28) logspace(1.48,3,35)]; % mm/ka
elseif simpleEcalc==1 && glacErate==0; % if using incremental erosion depth per glaciation (cm/glac)
    Eglacv = [(0:1:6) (8:2:20) (25:5:60) (70:10:190) (200:25:300) (350:50:500) (600:100:1000) 3000];
end;

% fix parameters for range calculations
if range1calc == 1 || range2calc == 1;
    % input parameter limits
    glacEmin_in = glacEmin;
    glacEmax_in = glacEmax;
    nonglacEmin_in = nonglacEmin;
    nonglacEmax_in = nonglacEmax;
    LR04min_in = LR04min;
    LR04max_in = LR04max;
    
    % set number of parameters
    glacEn = parn;
    nonglacEn = parn;
    LR04n = parn;
    if glacEmin == glacEmax; glacEn = 1; end;
    if nonglacEmin == nonglacEmax; nonglacEn = 1; end;
    if LR04min == LR04max; LR04n = 1; end;
    simn = glacEn*nonglacEn*LR04n;
end;

% for plotting
pl1d10 = []; pl1t10 = []; pl1d10uncmin = []; pl1d10uncmax = []; pl1m10 = []; pl1depth = [];
pl1d26 = []; pl1t26 = []; pl1d26uncmin = []; pl1d26uncmax = []; pl1m26 = []; plice = [];
pl2b10 = []; pl2t10 = []; pl2b10uncmin = []; pl2b10uncmax = []; pl2m10 = [];
pl2b26 = []; pl2t26 = []; pl2b26uncmin = []; pl2b26uncmax = []; pl2m26 = [];
pl3d10 = []; pl3t10 = []; pl3depth = []; pl3m10 = []; plicemin = []; plicemax = [];
pl3d26 = []; pl3t26 = []; pl3m26 = [];
pl3d1026 = []; pl3t1026 = []; pl3m1026 = [];
pl4b10 = []; pl4t10 = []; pl4m10 = [];
pl4b26 = []; pl4t26 = []; pl4m26 = [];
pl5b10 = []; pl5b26 = []; pl5t1026 = []; pl5m1026 = [];
pl6N10y = []; pl6N26y = []; pl6delN10y = []; pl6delN26y = []; pl6del2N10y = []; pl6del2N26y = [];
pl6N10n = []; pl6N26n = []; pl6delN10n = []; pl6delN26n = []; pl6del2N10n = []; pl6del2N26n = [];
pl6N10path = []; pl6N26path = [];
plEl.thick_rho = []; plEl.atm = []; plEl.RcEst = []; plEl.othercorr = [];
plEl.P10sp = []; plEl.P26sp = []; plEl.Lsp = []; plEl.P_mu10 = []; plEl.P_mu26 = [];

% fix for output
output(1,1) = {'sample'};
outn(1) = 1;
if sum(samplein.N10)>0 && simpleEcalc==1;
    if glacErate == 1; % if using erosion rate
        output(1,end+1:end+3) = {'10E(mm/ka)','10E+(mm/ka)','10E-(mm/ka)'};
    else; % if using erosion steps
        output(1,end+1:end+3) = {'10E(cm/glac)','10E+(cm/glac)','10E-(cm/glac)'};
    end;
    outn(1) = max(outn)+1; outn(2) = max(outn)+2;
    if mt >= 1E5;
        output(1,end+1) = {'10E-100ka(m)'};
        outn(3) = max(outn)+1;
    end;
    if mt >= 1E6;
        output(1,end+1) = {'10E-1Ma(m)'};
        outn(4) = max(outn)+1;
    end;
end;
if sum(samplein.N26)>0 && simpleEcalc==1;
    if glacErate == 1; % if using erosion rate
        output(1,end+1:end+3) = {'26E(mm/ka)','26E+(mm/ka)','26E-(mm/ka)'};
    else; % if using erosion steps
        output(1,end+1:end+3) = {'26E(cm/glac)','26E+(cm/glac)','26E-(cm/glac)'};
    end;
    outn(5) = max(outn)+1; outn(6) = max(outn)+2;
    if mt >= 1E5;
        output(1,end+1) = {'26E-100ka(m)'};
        outn(7) = max(outn)+1;
    end;
    if mt >= 1E6;
        output(1,end+1) = {'26E-1Ma(m)'};
        outn(8) = max(outn)+1;
    end;
end;
if sum(samplein.N10)>0 && range1calc == 1;
    if glacErate == 1; % if using erosion rate
        output(1,end+1:end+6) = {'glacEmin10(mm/ka)','glacEmax10(mm/ka)',...
            'nonglacEmin10(mm/ka)','nonglacEmax10(mm/ka)','LR04min10(‰)','LR04max10(‰)'};
    else; % if using erosion steps
        output(1,end+1:end+6) = {'glacEmin10(cm/glac)','glacEmax10(cm/glac)',...
            'nonglacEmin10(mm/ka)','nonglacEmax10(mm/ka)','LR04min10(‰)','LR04max10(‰)'};
    end;
    outn(9) = max(outn)+1; outn(10) = max(outn)+5;
    if mt >= 1E5;
        output(1,end+1:end+2) = {'10Emin-100ka(m)','10Emax-100ka(m)'};
        outn(11) = max(outn)+1; outn(12) = max(outn)+1;
    end;
    if mt >= 1E6;
        output(1,end+1:end+2) = {'10Emin-1Ma(m)','10Emax-1Ma(m)'};
        outn(13) = max(outn)+1; outn(14) = max(outn)+1;
    end;
end;
if sum(samplein.N26)>0 && range1calc == 1;
    if glacErate == 1; % if using erosion rate
        output(1,end+1:end+6) = {'glacEmin26(mm/ka)','glacEmax26(mm/ka)',...
            'nonglacEmin26(mm/ka)','nonglacEmax26(mm/ka)','LR04min26(‰)','LR04max26(‰)'};
    else; % if using erosion steps
        output(1,end+1:end+6) = {'glacEmin26(cm/glac)','glacEmax26(cm/glac)',...
            'nonglacEmin26(mm/ka)','nonglacEmax26(mm/ka)','LR04min26(‰)','LR04max26(‰)'};
    end;
    outn(15) = max(outn)+1; outn(16) = max(outn)+5;
    if mt >= 1E5;
        output(1,end+1:end+2) = {'26Emin-100ka(m)','26Emax-100ka(m)'};
        outn(17) = max(outn)+1; outn(18) = max(outn)+1;
    end;
    if mt >= 1E6;
        output(1,end+1:end+2) = {'26Emin-1Ma(m)','26Emax-1Ma(m)'};
        outn(19) = max(outn)+1; outn(20) = max(outn)+1;
    end;
end;
if sum(samplein.N10.*samplein.N26)>0 && range2calc == 1;
    if glacErate == 1; % if using erosion rate
        output(1,end+1:end+6) = {'glacEmin1026(mm/ka)','glacEmax1026(mm/ka)',...
            'nonglacEmin1026(mm/ka)','nonglacEmax1026(mm/ka)','LR04min1026(‰)','LR04max1026(‰)'};
    else; % if using erosion steps
        output(1,end+1:end+6) = {'glacEmin1026(cm/glac)','glacEmax1026(cm/glac)',...
            'nonglacEmin1026(mm/ka)','nonglacEmax1026(mm/ka)','LR04min1026(‰)','LR04max10(‰)'};
    end;
    outn(21) = max(outn)+1; outn(22) = max(outn)+5;
    if mt >= 1E5;
        output(1,end+1:end+2) = {'1026Emin-100ka(m)','1026Emax-100ka(m)'};
        outn(23) = max(outn)+1; outn(24) = max(outn)+1;
    end;
    if mt >= 1E6;
        output(1,end+1:end+2) = {'1026Emin-1Ma(m)','1026Emax-1Ma(m)'};
        outn(25) = max(outn)+1; outn(26) = max(outn)+1;
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
    sample.N26 = samplein.N26(i);
    sample.delN26 = samplein.delN26(i);
    sample.samplingyr = samplein.samplingyr(i);
    if isfield(samplein,'deglac'); sample.deglac = samplein.deglac(i); end;
    
    % write sample name to output
    output(i+1,1) = sample.sample_name;
    
    % Set nucl to 0 for both 10/26 and check if there is N10/N26
    nucl10 = 0; nucl26 = 0;
    if sample.N10 > 0; nucl10 = 1; end;
    if sample.N26 > 0; nucl26 = 1; end;
    
    if nucl10+nucl26==0 || (simpleEcalc==0 && range1calc==0 && nucl10*nucl26==0);
        continue;
    end;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample_name{1});
    
    % Age Relative to t0=2010 - LSD tv from LSDfix
    % tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
    
    % Fix w,Rc,SPhi, for sp and mu prod rate scaling
    LSDfix = LSD_fix(sample.lat,sample.long,mt+2010-sample.samplingyr,-1,consts);
    
    % adjust tv, Rc, SPhi, and d18O to sampling year
    if sample.samplingyr <= 2010;
        clipidx = min(find(LSDfix.tv > 2010-sample.samplingyr));
        tv = [2010-sample.samplingyr LSDfix.tv(clipidx:end)];
        Rc = interp1(LSDfix.tv,LSDfix.Rc,tv);
        SPhi = interp1(LSDfix.tv,LSDfix.SPhi,tv);
        d18O = interp1(LR04_tv,d18Oin,tv);
        tv = tv - 2010 + sample.samplingyr;
    else; % assume 2010 value for all years >2010
        Rc = [LSDfix.Rc(1) LSDfix.Rc];
        SPhi = [LSDfix.SPhi(1) LSDfix.SPhi];
        tv = [0 (LSDfix.tv + sample.samplingyr - 2010)];
        d18O = interp1(LR04_tv,d18Oin,LSDfix.tv);
        d18O = [d18O(1) d18O];
    end;
    
    % Production from muons
    % Precompute P_mu(z) to ~200,000 g/cm2
    % start at the mid-depth of the sample.
    sample.z_mu = [0 logspace(0,5.3,100)];
    P_mu_z = P_mu_expage(sample.z_mu+(sample.thick.*sample.rho./2),sample.pressure,LSDfix.RcEst,...
        consts.SPhiInf,nucl10,nucl26,consts,'no');
    
    % spallation production scaling
    Psp = LSDspal(sample.pressure,Rc,SPhi,LSDfix.w,nucl10,nucl26,consts);        
    
    % interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,Rc);
    
    % Thickness scaling factor.
    if sample.thick > 0;
        sample.thickSF = (sample.Lsp./(sample.rho.*sample.thick)).*...
            (1 - exp(((-1.*sample.rho.*sample.thick)./sample.Lsp)));
    else;
        sample.thickSF = 1;
    end;
    
    % find idx for deglaciation based on sample.deglac
    if isfield(sample,'deglac');
        d18Oidx = min(find(d18O == 5.02)); % find index of LGM max ice volume (at 18060)
        iceidx = min(find(tv >= sample.deglac));
    end;
    
    % if calculating single nuclide erosion using interpolation ====================================
    if simpleEcalc == 1;
        % fix ice cover record
        icev = (d18O >= LR04ice);
        icev = double(icev); % fix for matlab
        
        % fix deglaciation in icev using iceidx and d18Oidx
        if isfield(sample,'deglac');
            if sample.deglac-sample.samplingyr+2010 < 18060;
                icev(iceidx:d18Oidx) = 1;
            end;
            icev(1:iceidx-1) = 0;
        end;
        
        % fix specified ice-cover periods
        if exist('fix_ice');
            for j = 1:size(fix_ice,1);
                % find tv index
                fixmin = min(find(tv >= fix_ice(j,1)));
                fixmax = max(find(tv <= fix_ice(j,2)));
                icev(fixmin:fixmax) = 1;
            end;
        end;
        % fix specified ice-free periods
        if exist('fix_nonice');
            for j = 1:size(fix_nonice,1);
                % find tv index
                fixmin = min(find(tv >= fix_nonice(j,1)));
                fixmax = max(find(tv <= fix_nonice(j,2)));
                icev(fixmin:fixmax) = 0;
            end;
        end;
        
        icem = repmat(icev',1,numel(Eglacv));
        plice(end+1,1:tv(end)+1) = interp1(tv,icev,(0:1:tv(end)));
        
        % fix noice vector and matrix
        noicev = (icev == 0);
        noicem = (icem == 0);
        
        % fix non-glacial erosion vector
        nonglacEm = repmat(sample.E,numel(tv),numel(Eglacv));
        
        % fix glacial matrix
        if glacErate == 1;
            glacErm = repmat(Eglacv,numel(tv),1);
            glacEim = [];
        else;
            glacEim = repmat(Eglacv,numel(tv),1);
            glacErm = [];
        end;
        
        % fix parameters in specific time periods
        if exist('fix_tv_nonglacE') && exist('fix_nonglacE');
            for j = 1:size(fix_tv_nonglacE,1);
                % find tv index
                fixmin = min(find(tv >= fix_tv_nonglacE(j,1)));
                fixmax = max(find(tv <= fix_tv_nonglacE(j,2)));
                nonglacEm(fixmin:fixmax,:) = fix_nonglacE(j);
            end;
        end;
        if exist('fix_tv_glacEr') && exist('fix_glacEr');
            if numel(glacErm) == 0;
                glacErm = zeros(size(glacEim));
            end;
            for j = 1:size(fix_tv_glacEr,1);
                % find tv index
                fixmin = min(find(tv >= fix_tv_glacEr(j,1)));
                fixmax = max(find(tv <= fix_tv_glacEr(j,2)));
                glacErm(fixmin:fixmax,:) = fix_glacEr(j);
            end;
        end;
        if exist('fix_tv_glacEi') && exist('fix_glacEi');
            if numel(glacEim) == 0;
                glacEim = zeros(size(glacErm));
            end;
            for j = 1:size(fix_tv_glacEi,1);
                % find tv index
                fixmin = min(find(tv >= fix_tv_glacEi(j,1)));
                fixmax = max(find(tv <= fix_tv_glacEi(j,2)));
                glacEim(fixmin:fixmax,:) = fix_glacEi(j);
            end;
        end;
        
        % calculate non-glacial erosion depth
        nonglacdm = cumtrapz(tv',noicem.*nonglacEm.*1E-4).*sample.rho; % g/cm2
        nonglacdm = nonglacdm + burialdepth.*sample.rho; % add burial depth (g/cm2)
        
        % calculate glacial erosion depth matrix
        glacdrm = 0;
        glacdim = 0;
        if numel(glacErm) > 0;
            glacdrm = cumtrapz(tv',icem.*glacErm.*1E-4).*sample.rho; % g/cm2
        end;
        if numel(glacEim) > 0;
            % find deglaciation events
            deglacm = (filter([1 2],1,icem) == 1);
            glacdim = cumsum(deglacm.*glacEim).*sample.rho; % g/cm2
        end;
        
        % full depth matrix (glac and nonglac)
        fulldm = glacdrm + glacdim + nonglacdm;
        
        % depth matrix (cm)
        ddm = fulldm./sample.rho;
        
        % 10Be single nuclide calculations
        if nucl10 == 1;
            % calculate N(t) for matrix
            Nend10 = Ncalc(fulldm,noicem,tv,Psp.sp10,Pref10,P_mu_z.mu10,l10,sample);
            
            % calculate erosion and uncertainties (and display output)
            out10 = fix_output1(sample.N10,sample.delN10,Pref10,delPref10,Nend10,Eglacv,...
                glacErate,'10Be');
            output(i+1,outn(1):outn(2)) = out10.outstr;
            
            % calculate depth history
            dv10 = interp1(Eglacv',ddm',out10.erosion,'pchip')./1E2; % interpolated depth (m)
            if max(tv) >= 1E5;
                output(i+1,outn(3)) = {num2str(interp1(tv,dv10,1E5,'pchip'),'%.2f')}; % d(m) 100 ka
            end;
            if max(tv) >= 1E6;
                output(i+1,outn(4)) = {num2str(interp1(tv,dv10,1E6,'pchip'),'%.2f')}; % d(m) 1 Ma
            end;
            
            % fix sample depth plot matrix
            if pl_simpleEcalc_depth == 1;
                % fill matrix
                pl1d10(1:numel(tv),end+1) = dv10;
                pl1t10(1:numel(tv),end+1) = tv'./1E3;
                pl1m10(1:numel(tv),end+1) = 1;
                pl1depth(end+1) = interp1(tv,dv10,min(pl_maxt,max(tv)),'pchip'); % add d for y axis
                % fix uncertainties
                if simpleEcalcunc==1 && isfield(out10,'mindelE'); % interpret max depth (m)
                    uncmin10 = interp1(Eglacv',ddm',out10.erosion-out10.mindelE,'pchip')./1E2;
                    uncmax10 = interp1(Eglacv',ddm',out10.erosion+out10.maxdelE,'pchip')./1E2;
                    pl1d10uncmin(1:numel(tv),end+1) = uncmin10;
                    pl1d10uncmax(1:numel(tv),end+1) = uncmax10;
                    pl1depth(end+1) = interp1(tv,uncmax10,pl_maxt,'pchip'); % max d for y-axis
                end;
            end
            
            % fix nuclide buildup matrix
            if pl_simpleEcalc_buildup == 1;
                % calculate P through time
                Pfull10 = nucl1_getP(Eglacv,fulldm,out10,sample,Psp.sp10,Pref10,P_mu_z.mu10,noicev);
                % calculate N through time
                Npath = Npathcalc(tv,Pfull10,l10);
                % fill matrix
                Npath10 = Npath.N10./repmat(Npath.N10(1,:),numel(tv),1).*1E2;
                pl2b10(1:numel(tv),end+1) = Npath10(:,1);
                pl2t10(1:numel(tv),end+1) = tv'./1E3;
                pl2m10(1:numel(tv),end+1) = 1;
                % fix uncertainties
                if simpleEcalcunc==1 && isfield(out10,'mindelE');
                    pl2b10uncmax(1:numel(tv),end+1) = Npath10(:,2);
                    pl2b10uncmin(1:numel(tv),end+1) = Npath10(:,3);
                end;
            end;
        end; % end 10Be calculations
        
        % 26Al single nuclide calculations
        if nucl26 == 1;
            % calculate N(t) for matrix
            Nend26 = Ncalc(fulldm,noicem,tv,Psp.sp26,Pref26,P_mu_z.mu26,l26,sample);
            
            % calculate erosion and uncertainties (and display output)
            out26 = fix_output1(sample.N26,sample.delN26,Pref26,delPref26,Nend26,Eglacv,...
                glacErate,'26Al');
            output(i+1,outn(5):outn(6)) = out26.outstr;
            
            % calculate depth history
            dv26 = interp1(Eglacv',ddm',out26.erosion,'pchip')./1E2; % interpolated depth (m)
            if max(tv) >= 1E5;
                output(i+1,outn(7)) = {num2str(interp1(tv,dv26,1E5,'pchip'),'%.2f')}; % d(m) 100 ka
            end;
            if max(tv) >= 1E6;
                output(i+1,outn(8)) = {num2str(interp1(tv,dv26,1E6,'pchip'),'%.2f')}; % d(m) 1 Ma
            end;
            
            % fix sample depth plot matrix
            if pl_simpleEcalc_depth == 1;
                % fill matrix
                pl1d26(1:numel(tv),end+1) = dv26;
                pl1t26(1:numel(tv),end+1) = tv'./1E3;
                pl1m26(1:numel(tv),end+1) = 1;
                pl1depth(end+1) = interp1(tv,dv26,min(pl_maxt,max(tv)),'pchip'); % add d for y axis
                % fix uncertainties
                if simpleEcalcunc==1 && isfield(out26,'mindelE'); % interpret max depth (m)
                    uncmin26 = interp1(Eglacv',ddm',out26.erosion-out26.mindelE,'pchip')./1E2;
                    uncmax26 = interp1(Eglacv',ddm',out26.erosion+out26.maxdelE,'pchip')./1E2;
                    pl1d26uncmin(1:numel(tv),end+1) = uncmin26;
                    pl1d26uncmax(1:numel(tv),end+1) = uncmax26;
                    pl1depth(end+1) = interp1(tv,uncmax26,pl_maxt,'pchip'); % max d for y-axis
                end;
            end
            
            % fix nuclide buildup matrix
            if pl_simpleEcalc_buildup == 1;
                % calculate P through time
                Pfull26 = nucl1_getP(Eglacv,fulldm,out26,sample,Psp.sp26,Pref26,P_mu_z.mu26,noicev);
                % calculate N through time
                Npath = Npathcalc(tv,Pfull26,l26);
                % fill matrix
                Npath26 = Npath.N26./repmat(Npath.N26(1,:),numel(tv),1).*1E2;
                pl2b26(1:numel(tv),end+1) = Npath26(:,1);
                pl2t26(1:numel(tv),end+1) = tv'./1E3;
                pl2m26(1:numel(tv),end+1) = 1;
                % fix uncertainties
                if simpleEcalcunc==1 && isfield(out26,'mindelE');
                    pl2b26uncmax(1:numel(tv),end+1) = Npath26(:,2);
                    pl2b26uncmin(1:numel(tv),end+1) = Npath26(:,3);
                end;
            end;
        end; % end 26Al calculations
    end; % end simpleEcalc
    % ==============================================================================================
    
    % if calculating ranges ========================================================================
    if range1calc==1 || (range2calc==1 && nucl10==1 && nucl26==1);      
        % matrix/vectors to be filled in while loop
        ddmok = []; glacEok = []; nonglacEok = []; LR04ok = []; Ndiffok = []; icemok = [];
        P10ok = []; P26ok = []; N10ok = []; N26ok = []; nhit = 0;
        
        % parameters for random runs
        randrun = 0; randhit = 0; randnum = 0; rand1n = 0;
        
        % parameter tests
        partest(1:6) = 0;
        randtest(1:6) = 0;
        
        % fix looptest parameter
        looptest = [1 1 1];
        if nucl10 == 1 && range1calc == 1;
            looptest(1) = 0;
        end;
        if nucl26 == 1 && range1calc == 1;
            looptest(2) = 0;
        end;
        if nucl10 == 1 && nucl26 == 1 && range2calc == 1;
            looptest(3) = 0;
        end;
        
        % display nuclide(s) calculating ranges
        if min(find(looptest==0)) == 1;
            nuclides = '10Be';
        elseif min(find(looptest==0)) == 2;
            nuclides = '26Al';
        elseif min(find(looptest==0)) == 3;
            nuclides = '10Be+26Al';
        end;
        fprintf(1,'\n%s calculating ranges',nuclides);
        
        % add prodrate uncertainty
        if Punc == 1 && sample.N10 > 0;
            sample.delN10 = sqrt(sample.delN10^2 + (sample.N10*delPref10/Pref10)^2);
        end;
        if Punc == 1 && sample.N26 > 0;
            sample.delN26 = sqrt(sample.delN26^2 + (sample.N26*delPref26/Pref26)^2);
        end;
        
        % while loop
        while prod(looptest) == 0;
            fprintf(1,'.');
            
            if randrun == 0; % all runs except the last
                % fix parameter vectors
                glacEv = repmat(linspace(glacEmin,glacEmax,glacEn),LR04n*nonglacEn,1);
                glacEv = glacEv(:)';
                nonglacEv = repmat(repmat(linspace(nonglacEmin,nonglacEmax,nonglacEn),1,LR04n),1,...
                    glacEn);
                LR04v = repmat(linspace(LR04min,LR04max,LR04n),nonglacEn,1);
                LR04v = repmat(LR04v(:)',1,glacEn);
            else; % random parameters in last runs
                glacEv = glacEmin+rand(1,randrunsz).*(glacEmax-glacEmin);
                nonglacEv = nonglacEmin+rand(1,randrunsz).*(nonglacEmax-nonglacEmin);
                LR04v = LR04min+rand(1,randrunsz).*(LR04max-LR04min);
                randnum = randnum + 1;
            end;
            
            % calculate parameter steps
            glacEst = (glacEmax-glacEmin)/max(glacEn-1,1);
            nonglacEst = (nonglacEmax-nonglacEmin)/max(nonglacEn-1,1);
            LR04st = (LR04max-LR04min)/max(LR04n-1,1);
            
            % ice cover matrix
            icem = (repmat(d18O',size(LR04v)) >= LR04v);
            
            % fix deglaciation in icem using iceidx and d18Oidx
            if isfield(sample,'deglac');
                if sample.deglac-sample.samplingyr+2010 < 18060;
                    icem(iceidx:d18Oidx,:) = 1;
                end;
                icem(1:iceidx-1,:) = 0;
            end;
            
            % fix specified ice-cover periods
            if exist('fix_ice');
                for j = 1:size(fix_ice,1);
                    % find tv index
                    fixmin = min(find(tv >= fix_ice(j,1)));
                    fixmax = max(find(tv <= fix_ice(j,2)));
                    icem(fixmin:fixmax,:) = 1;
                end;
            end;
            % fix specified ice-free periods
            if exist('fix_nonice');
                for j = 1:size(fix_nonice,1);
                    % find tv index
                    fixmin = min(find(tv >= fix_nonice(j,1)));
                    fixmax = max(find(tv <= fix_nonice(j,2)));
                    icem(fixmin:fixmax,:) = 0;
                end;
            end;
            
            % fix noice matrix
            noicem = (icem == 0);
            
            % nonglacE matrix
            nonglacEm = repmat(nonglacEv,numel(tv),1);
            
            % fix glacial matrix
            if glacErate == 1;
                glacErm = repmat(glacEv,numel(tv),1);
                glacEim = [];
            else;
                glacEim = repmat(glacEv,numel(tv),1);
                glacErm = [];
            end;
            
            % fix parameters in specific time periods
            if exist('fix_tv_nonglacE') && exist('fix_nonglacE');
                for j = 1:size(fix_tv_nonglacE,1);
                    % find tv index
                    fixmin = min(find(tv >= fix_tv_nonglacE(j,1)));
                    fixmax = max(find(tv <= fix_tv_nonglacE(j,2)));
                    nonglacEm(fixmin:fixmax,:) = fix_nonglacE(j);
                end;
            end;
            if exist('fix_tv_glacEr') && exist('fix_glacEr');
                if numel(glacErm) == 0;
                    glacErm = zeros(size(glacEim));
                end;
                for j = 1:size(fix_tv_glacEr,1);
                    % find tv index
                    fixmin = min(find(tv >= fix_tv_glacEr(j,1)));
                    fixmax = max(find(tv <= fix_tv_glacEr(j,2)));
                    glacErm(fixmin:fixmax,:) = fix_glacEr(j);
                end;
            end;
            if exist('fix_tv_glacEi') && exist('fix_glacEi');
                if numel(glacEim) == 0;
                    glacEim = zeros(size(glacErm));
                end;
                for j = 1:size(fix_tv_glacEi,1);
                    % find tv index
                    fixmin = min(find(tv >= fix_tv_glacEi(j,1)));
                    fixmax = max(find(tv <= fix_tv_glacEi(j,2)));
                    glacEim(fixmin:fixmax,:) = fix_glacEi(j);
                end;
            end;
            
            % calculate non-glacial erosion depth
            nonglacdm = cumtrapz(tv',noicem.*nonglacEm.*1E-4).*sample.rho; % g/cm2
            nonglacdm = nonglacdm + burialdepth.*sample.rho; % add burial depth (g/cm2)
            
            % calculate glacial erosion depth matrix
            glacdrm = 0;
            glacdim = 0;
            if numel(glacErm) > 0;
                glacdrm = cumtrapz(tv',icem.*glacErm.*1E-4).*sample.rho; % g/cm2
            end;
            if numel(glacEim) > 0;
                % find deglaciation events
                deglacm = (filter([1 2],1,icem) == 1);
                glacdim = cumsum(deglacm.*glacEim).*sample.rho; % g/cm2
            end;
            
            % full depth matrix (glac and nonglac)
            fulldm = glacdrm + glacdim + nonglacdm;
            
            % depth matrix (cm)
            ddm = fulldm./sample.rho;
            
            % do calculation and pick out runs that yield OK conc
            if looptest(1) == 0; % 10Be single nuclide range
                % calculate N(t) for matrix
                Nend = Ncalc(fulldm,noicem,tv,Psp.sp10,Pref10,P_mu_z.mu10,l10,sample);
                
                % pick out concentrations within measured uncertainty
                Ndiff = abs(Nend-sample.N10);
                hitidx = find(Ndiff <= sample.delN10);
                
                % calculate P10 for OK scenarios
                P10m = Pcalc1026(sample,fulldm(:,hitidx),noicem(:,hitidx),tv,Psp.sp10,Pref10,...
                    P_mu_z.mu10,l10);
                % pick out OK P10 vectors
                P10ok(:,end+1:end+numel(hitidx)) = P10m.P10;
            elseif looptest(2) == 0; % 26Al single nuclide range
                % calculate N(t) for matrix
                Nend = Ncalc(fulldm,noicem,tv,Psp.sp26,Pref26,P_mu_z.mu26,l26,sample);
                
                % pick out concentrations within measured uncertainty
                Ndiff = abs(Nend-sample.N26);
                hitidx = find(Ndiff <= sample.delN26);
                
                % calculate P26 for OK scenarios
                P26m = Pcalc1026(sample,fulldm(:,hitidx),noicem(:,hitidx),tv,Psp.sp26,Pref26,...
                    P_mu_z.mu26,l26);
                % pick out OK P26 vectors
                P26ok(:,end+1:end+numel(hitidx)) = P26m.P26;
            elseif looptest(3) == 0; % 26Al/10Be nuclide range
                % calculate N(t) for matrix
                Nend10 = Ncalc(fulldm,noicem,tv,Psp.sp10,Pref10,P_mu_z.mu10,l10,sample);
                Nend26 = Ncalc(fulldm,noicem,tv,Psp.sp26,Pref26,P_mu_z.mu26,l26,sample);
                
                % pick out concentrations within measured uncertainty
                Ndiff = abs(Nend10-sample.N10)./sample.delN10 + ...
                    abs(Nend26-sample.N26)./sample.delN26;
                hitidx = find((abs(Nend10-sample.N10) <= sample.delN10) .* ...
                    (abs(Nend26-sample.N26) <= sample.delN26) == 1);
                
                % calculate P10 and P26 for matrix
                P1026m = Pcalc1026(sample,fulldm(:,hitidx),noicem(:,hitidx),tv,Psp.sp10,Pref10,...
                    P_mu_z.mu10,l10,Psp.sp26,Pref26,P_mu_z.mu26,l26);
                
                % pick out concentrations (for ratio path plot selection)
                N10ok(end+1:end+numel(hitidx)) = Nend10(hitidx);
                N26ok(end+1:end+numel(hitidx)) = Nend26(hitidx);
                
                % pick out OK P10 and P26 vectors
                P10ok(:,end+1:end+numel(hitidx)) = P1026m.P10;
                P26ok(:,end+1:end+numel(hitidx)) = P1026m.P26;
            end;
            
            % pick out OK depth vectors from ddm
            ddmok(:,end+1:end+numel(hitidx)) = ddm(:,hitidx);
            
            % pick out OK ice cover periods from icem
            icemok(:,end+1:end+numel(hitidx)) = icem(:,hitidx);
            
            % pick out parameters yielding OK concentrations + Ncomp
            glacEhit = glacEv(hitidx);
            nonglacEhit = nonglacEv(hitidx);
            LR04hit = LR04v(hitidx);
            Ndiffhit = Ndiff(hitidx);
            
            % fill ok vectors
            glacEok(end+1:end+numel(hitidx)) = glacEhit;
            nonglacEok(end+1:end+numel(hitidx)) = nonglacEhit;
            LR04ok(end+1:end+numel(hitidx)) = LR04hit;
            Ndiffok(end+1:end+numel(hitidx)) = Ndiffhit;
            
            % find parameter with largest relative step (1 2 3)
            [parmax paridx] = max([glacEst/glacEstlim nonglacEst/nonglacEstlim LR04st/LR04stlim]);
                    
            % if any hit
            if numel(hitidx) > 0 && randrun == 0;
                % if first hit
                if nhit == 0;
                    % save hit parameters
                    glacEmin1 = glacEmin;
                    glacEmax1 = glacEmax;
                    nonglacEmin1 = nonglacEmin;
                    nonglacEmax1 = nonglacEmax;
                    LR04min1 = LR04min;
                    LR04max1 = LR04max;
                end;
                
                % test for present parameter
                parcltest = 0;
                
                % test if parameters min/max limits matches input min/max
                if partest(1)==0 && min(glacEok)==glacEmin_in;
                    if min(find(partest==0)) == 1; parcltest = 1; end;
                    partest(1) = 1; randtest(1) = 1;
                end;
                if partest(2)==0 && max(glacEok)==glacEmax_in;
                    if min(find(partest==0)) == 2; parcltest = 1; end;
                    partest(2) = 1; randtest(2) = 1;
                end;
                if partest(3)==0 && min(nonglacEok)==nonglacEmin_in;
                    if min(find(partest==0)) == 3; parcltest = 1; end;
                    partest(3) = 1; randtest(3) = 1;
                end;
                if partest(4)==0 && max(nonglacEok)==nonglacEmax_in;
                    if min(find(partest==0)) == 4; parcltest = 1; end;
                    partest(4) = 1; randtest(4) = 1;
                end;
                if partest(5)==0 && min(LR04ok)==LR04min_in;
                    if min(find(partest==0)) == 5; parcltest = 1; end;
                    partest(5) = 1; randtest(5) = 1;
                end;
                if partest(6)==0 && max(LR04ok)==LR04max_in;
                    if min(find(partest==0)) == 6; parcltest = 1; end;
                    partest(6) = 1; randtest(6) = 1;
                end;
                
                if partest(1) == 0;
                    % test if parameter has been found
                    if abs(mean([glacEmin glacEmax])-min(glacEok))<glacEstlim && ...
                        glacEst<=glacEstlim && nonglacEst<=nonglacEstlim && LR04st<=LR04stlim;
                        partest(1) = 1;
                        parcltest = 1;
                    else; % pick out parameter limit idx
                        limidx = find(glacEok == min(glacEok));
                    end;
                elseif partest(2) == 0;
                    % test if parameter has been found
                    if abs(mean([glacEmin glacEmax])-max(glacEok))<glacEstlim && ...
                        glacEst<=glacEstlim && nonglacEst<=nonglacEstlim && LR04st<=LR04stlim;
                        partest(2) = 1;
                        parcltest = 1;
                    else; % pick out parameter limit idx
                        limidx = find(glacEok == max(glacEok));
                    end;
                elseif partest(3) == 0;
                    % test if parameter has been found
                    if abs(mean([nonglacEmin nonglacEmax])-min(nonglacEok))<nonglacEstlim && ...
                        glacEst<=glacEstlim && nonglacEst<=nonglacEstlim && LR04st<=LR04stlim;
                        partest(3) = 1;
                        parcltest = 1;
                    else; % pick out parameter limit idx
                        limidx = find(nonglacEok == min(nonglacEok));
                    end;
                elseif partest(4) == 0;
                    % test if parameter has been found
                    if abs(mean([nonglacEmin nonglacEmax])-max(nonglacEok))<nonglacEstlim && ...
                        glacEst<=glacEstlim && nonglacEst<=nonglacEstlim && LR04st<=LR04stlim;
                        partest(4) = 1;
                        parcltest = 1;
                    else; % pick out parameter limit idx
                        limidx = find(nonglacEok == max(nonglacEok));
                    end;
                elseif partest(5) == 0;
                    % test if parameter has been found
                    if abs(mean([LR04min LR04max])-min(LR04ok))<LR04stlim && ...
                        glacEst<=glacEstlim && nonglacEst<=nonglacEstlim && LR04st<=LR04stlim;
                        partest(5) = 1;
                        parcltest = 1;
                    else; % pick out parameter limit idx
                        limidx = find(LR04ok == min(LR04ok));
                    end;
                elseif partest(6) == 0;
                    % test if parameter has been found
                    if abs(mean([LR04min LR04max])-max(LR04ok))<LR04stlim && ...
                        glacEst<=glacEstlim && nonglacEst<=nonglacEstlim && LR04st<=LR04stlim;
                        partest(6) = 1;
                        parcltest = 1;
                    else; % pick out parameter limit idx
                        limidx = find(LR04ok == max(LR04ok));
                    end;
                end;
                
                % if the actual parameter limit has not been found
                if parcltest == 0;
                    % find index of best ok parameter limit
                    [Ndiffbest idxb] = min(Ndiffok(limidx));
                    idxb = limidx(idxb);
                    
                    % fix new limits for glacE
                    if paridx == 1; % if glacE has the largest relative step
                        glacEmin = max(glacEok(idxb)-glacEst*(glacEn-1)/4,glacEmin_in);
                        glacEmax = min(glacEok(idxb)+glacEst*(glacEn-1)/4,glacEmax_in);
                    elseif glacEmax-glacEmin < glacEmax_in-glacEmin_in;
                        glacEmin = max(glacEok(idxb)-glacEst*(glacEn-1)/2,glacEmin_in);
                        glacEmax = min(glacEok(idxb)+glacEst*(glacEn-1)/2,glacEmax_in);
                    end;
                    
                    % fix new limits for nonglacE
                    if paridx == 2; % if nonglacE has the largest relative step
                        nonglacEmin = max(nonglacEok(idxb)-nonglacEst*(nonglacEn-1)/4,...
                            nonglacEmin_in);
                        nonglacEmax = min(nonglacEok(idxb)+nonglacEst*(nonglacEn-1)/4,...
                            nonglacEmax_in);
                    elseif nonglacEmax-nonglacEmin < nonglacEmax_in-nonglacEmin_in;
                        nonglacEmin = max(nonglacEok(idxb)-nonglacEst*(glacEn-1)/2,nonglacEmin_in);
                        nonglacEmax = min(nonglacEok(idxb)+nonglacEst*(glacEn-1)/2,nonglacEmax_in);
                    end;
                    
                    % fix new limits for LR04
                    if paridx == 3; % if LR04 has the largest relative step
                        LR04min = max(LR04ok(idxb)-LR04st*(LR04n-1)/4,LR04min_in);
                        LR04max = min(LR04ok(idxb)+LR04st*(LR04n-1)/4,LR04max_in);
                    elseif LR04max-LR04min < LR04max_in-LR04min_in;
                        LR04min = max(LR04ok(idxb)-LR04st*(LR04n-1)/2,LR04min_in);
                        LR04max = min(LR04ok(idxb)+LR04st*(LR04n-1)/2,LR04max_in);
                    end;
                elseif prod(partest) == 0; % if there are more parameters to be found
                    % use first hit parameters
                    glacEmin = glacEmin1;
                    glacEmax = glacEmax1;
                    nonglacEmin = nonglacEmin1;
                    nonglacEmax = nonglacEmax1;
                    LR04min = LR04min1;
                    LR04max = LR04max1;
                end;
                
                % add 1 to nhit
                nhit = nhit + 1;
                
            elseif LR04st>LR04stlim || nonglacEst>nonglacEstlim || glacEst>glacEstlim;
                % if no hit: set parameter min and max values around best fit in last run
                % find which parameter to change based on step size and step size limit
                [parmax paridx] = max([glacEst/glacEstlim nonglacEst/nonglacEstlim ...
                    LR04st/LR04stlim]);
                
                % if no hit: set parameter min and max values around best fit in last run
                [Ndiffbest bestidx] = min(Ndiff);
                
                % do parameter limit change based on bestidx
                if paridx==1 && glacEst>glacEstlim;
                    glacEmin = max(glacEv(bestidx)-glacEst*(glacEn-1)/4,glacEmin_in);
                    glacEmax = min(glacEv(bestidx)+glacEst*(glacEn-1)/4,glacEmax_in);
                elseif paridx==2 && nonglacEst>nonglacEstlim;
                    nonglacEmin = max(nonglacEv(bestidx)-nonglacEst*(nonglacEn-1)/4,nonglacEmin_in);
                    nonglacEmax = min(nonglacEv(bestidx)+nonglacEst*(nonglacEn-1)/4,nonglacEmax_in);
                elseif paridx==3 && LR04st>LR04stlim;
                    LR04min = max(LR04v(bestidx)-LR04st*(LR04n-1)/4,LR04min_in);
                    LR04max = min(LR04v(bestidx)+LR04st*(LR04n-1)/4,LR04max_in);
                end;
            else; % if no hit and all parameter steps are smaller than max steps
                partest(:) = 1; % move on
            end;
            
            % if done with iterative parameter optimization - fix for random runs
            if  prod(partest)==1 && min(randtest)==0; % random runs for specific parameter limits
                if randtest(1)==0 || randtest(2)==0;
                    [sorted sortidx] = sort(glacEok);
                    if randtest(1) == 0;
                        sidx = sortidx(1:round(numel(glacEok)*0.1));
                        glacEmin = max(min(glacEok)*0.9,glacEmin_in);
                        glacEmax = min(glacEok);
                    else;
                        sidx = sortidx(round(numel(glacEok)*0.9):end);
                        glacEmin = max(glacEok);
                        glacEmax = min(max(glacEok)*1.1,glacEmax_in);
                    end;
                    nonglacEmin = min(nonglacEok(sidx));
                    nonglacEmax = max(nonglacEok(sidx));
                    LR04min = min(LR04ok(sidx));
                    LR04max = max(LR04ok(sidx));
                elseif randtest(3)==0 || randtest(4)==0;
                    [sorted sortidx] = sort(nonglacEok);
                    if randtest(3) == 0;
                        sidx = sortidx(1:round(numel(nonglacEok)*0.1));
                        nonglacEmin = max(min(nonglacEok)*0.9,nonglacEmin_in);
                        nonglacEmax = min(nonglacEok);
                    else;
                        sidx = sortidx(round(numel(nonglacEok)*0.9):end);
                        nonglacEmin = max(nonglacEok);
                        nonglacEmax = min(max(nonglacEok)*1.1,nonglacEmax_in);
                    end;
                    glacEmin = min(glacEok(sidx));
                    glacEmax = max(glacEok(sidx));
                    LR04min = min(LR04ok(sidx));
                    LR04max = max(LR04ok(sidx));
                elseif randtest(5)==0 || randtest(6)==0;
                    [sorted sortidx] = sort(LR04ok);
                    if randtest(5) == 0;
                        sidx = sortidx(1:round(numel(LR04ok)*0.1));
                        LR04min = max(min(LR04ok)*0.9,LR04min_in);
                        LR04max = min(LR04ok);
                    else;
                        sidx = sortidx(round(numel(LR04ok)*0.9):end);
                        LR04min = max(LR04ok);
                        LR04max = min(max(LR04ok)*1.1,LR04max_in);
                    end;
                    glacEmin = min(glacEok(sidx));
                    glacEmax = max(glacEok(sidx));
                    LR04min = min(LR04ok(sidx));
                    LR04max = max(LR04ok(sidx));
                end;
                rand1n = rand1n + 1;
                if rand1n == rand1lim;
                    randtest(min(find(randtest==0))) = 1;
                    rand1n = 0;
                end;
                randrun = 1;
            elseif  prod(partest)==1 && randhit<randhitN && nhit>0; % final random runs
                glacEmin = max(min(glacEok)*0.9,glacEmin_in);
                glacEmax = min(max(glacEok)*1.1,glacEmax_in);
                nonglacEmin = max(min(nonglacEok)*0.9,nonglacEmin_in);
                nonglacEmax = min(max(nonglacEok)*1.1,nonglacEmax_in);
                LR04min = max(min(LR04ok)*0.9,LR04min_in);
                LR04max = min(max(LR04ok)*1.1,LR04max_in);
                if randrun == 2;
                    randhit = randhit + numel(hitidx);
                end;
                randrun = 2;
                if randnum >= randnummax; % don't allow too many random runs...
                    randhit = randhitN;
                end;
                randnum = randnum + 1;
            end;
            
            % if done with parameter optimization
            if prod(partest)==1 && (randhit>=randhitN || nhit==0);
                % fix for output and plotting
                if min(find(looptest==0)) == 1;
                    outnr1 = outn(9); outnr2 = outn(10);
                    if mt >= 1E5; outnr3 = outn(11); outnr4 = outn(12); end;
                    if mt >= 1E6; outnr5 = outn(13); outnr6 = outn(14); end;
                elseif min(find(looptest==0)) == 2;
                    outnr1 = outn(15); outnr2 = outn(16);
                    if mt >= 1E5; outnr3 = outn(17); outnr4 = outn(18); end;
                    if mt >= 1E6; outnr5 = outn(19); outnr6 = outn(20); end;
                elseif min(find(looptest==0)) == 3;
                    outnr1 = outn(21); outnr2 = outn(22);
                    if mt >= 1E5; outnr3 = outn(23); outnr4 = outn(24); end;
                    if mt >= 1E6; outnr5 = outn(25); outnr6 = outn(26); end;
                end;
                % write output
                if numel(LR04ok) > 0;
                    output(i+1,outnr1:outnr2) = {num2str(min(glacEok),'%.2f'),...
                        num2str(max(glacEok),'%.2f'),num2str(min(nonglacEok),'%.2f'),...
                        num2str(max(nonglacEok),'%.2f'),num2str(min(LR04ok),'%.2f'),...
                        num2str(max(LR04ok),'%.2f')};
                    % display limits
                    if glacErate == 1; glacEunit = 'mm/ka'; else; glacEunit = 'cm/glac'; end;
                    % if you get a character error here: change the end of the string below and
                    % check that it is ended with a proper quotation mark: '
                    fprintf(1,...
                    '\nglacE: %.2f-%.2f %s   nonglacE: %.2f-%.2f mm/ka   LR04: %.2f-%.2f ‰',...
                    min(glacEok),max(glacEok),glacEunit,min(nonglacEok),max(nonglacEok),...
                    min(LR04ok),max(LR04ok));
                    
                    % calculate and write depth history
                    if max(tv) >= 1E5; % d(m) 100 ka
                        if size(ddmok,2) > 1; % if more than one hit
                            ddm_min = min(ddmok')./1E2; % min sample depth (m)
                            ddm_max = max(ddmok')./1E2; % max sample depth (m)
                        else;
                            ddm_min = ddmok'./1E2;
                            ddm_max = ddmok'./1E2;
                        end;
                        output(i+1,outnr3) = {num2str(interp1(tv,ddm_min,1E5,'pchip'),'%.2f')};
                        output(i+1,outnr4) = {num2str(interp1(tv,ddm_max,1E5,'pchip'),'%.2f')};
                    end;
                    if max(tv) >= 1E6; % d(m) 1 Ma
                        output(i+1,outnr5) = {num2str(interp1(tv,ddm_min,1E6,'pchip'),'%.2f')};
                        output(i+1,outnr6) = {num2str(interp1(tv,ddm_max,1E6,'pchip'),'%.2f')};
                    end;
                else;
                    output(i+1,outnr1:outnr2) = {'nohit'};
                    if max(tv) >= 1E5; output(i+1,outnr3:outnr4) = {'nohit'}; end;
                    if max(tv) >= 1E6; output(i+1,outnr5:outnr6) = {'nohit'}; end;
                    fprintf(1,'\nno hits possible!');
                end;
                
                % fix sample depth plot matrix
                if pl_range_depth==1 && nhit>0 && (range1calc==1 || range2calc==1);
                    if size(ddmok,2) > 1; % if more than one hit
                        ddm_min = min(ddmok');
                        ddm_max = max(ddmok');
                    else;
                        ddm_min = ddmok';
                        ddm_max = ddmok';
                    end;
                    % fill plot matrix
                    if min(find(looptest==0)) == 1;
                        pl3t10(1:numel(tv)*2,end+1) = [flip(tv');tv']./1E3;
                        pl3d10(1:numel(tv)*2,end+1) = [flip(ddm_min');ddm_max']./1E2;
                        pl3m10(1:numel(tv)*2,end+1) = 1;
                    elseif min(find(looptest==0)) == 2;
                        pl3t26(1:numel(tv)*2,end+1) = [flip(tv');tv']./1E3;
                        pl3d26(1:numel(tv)*2,end+1) = [flip(ddm_min');ddm_max']./1E2;
                        pl3m26(1:numel(tv)*2,end+1) = 1;
                    elseif min(find(looptest==0)) == 3;
                        pl3t1026(1:numel(tv)*2,end+1) = [flip(tv');tv']./1E3;
                        pl3d1026(1:numel(tv)*2,end+1) = [flip(ddm_min');ddm_max']./1E2;
                        pl3m1026(1:numel(tv)*2,end+1) = 1;
                    end;
                    % fix max depth for y-axis
                    pl3depth(end+1) = interp1(tv,ddm_max,min(pl_maxt,max(tv)),'pchip')./1E2;
                end;
                
                % fix nuclide buildup matrix
                if pl_range_buildup==1 && nhit>1 && (range1calc==1 || range2calc==1);
                    % calculate nuclide buildup
                    if min(find(looptest==0)) == 1;
                        Npath = Npathcalc(tv,P10ok,l10);
                        N10path = Npath.N10./repmat(Npath.N10(1,:),numel(tv),1).*1E2;
                        if size(N10path,2) == 1;
                            N10path(:,2) = N10path(:,1);
                        end;
                        pl4t10(1:numel(tv)*2,end+1) = [flip(tv');tv']./1E3;
                        pl4b10(1:numel(tv)*2,end+1) = [flip(max(N10path')) min(N10path')]';
                        pl4m10(1:numel(tv)*2,end+1) = 1;
                    elseif min(find(looptest==0)) == 2;
                        Npath = Npathcalc(tv,P26ok,l26);
                        N26path = Npath.N26./repmat(Npath.N26(1,:),numel(tv),1).*1E2;
                        if size(N26path,2) == 1;
                            N26path(:,2) = N26path(:,1);
                        end;
                        pl4t26(1:numel(tv)*2,end+1) = [flip(tv');tv']./1E3;
                        pl4b26(1:numel(tv)*2,end+1) = [flip(max(N26path')) min(N26path')]';
                        pl4m26(1:numel(tv)*2,end+1) = 1;
                    elseif min(find(looptest==0)) == 3;
                        Npath = Npathcalc(tv,P10ok,l10,P26ok,l26);
                        N10path = Npath.N10./repmat(Npath.N10(1,:),numel(tv),1).*1E2;
                        N26path = Npath.N26./repmat(Npath.N26(1,:),numel(tv),1).*1E2;
                        if size(N10path,2) == 1;
                            N10path(:,2) = N10path(:,1);
                            N26path(:,2) = N26path(:,1);
                        end;
                        pl5t1026(1:numel(tv)*2,end+1) = [flip(tv');tv']./1E3;
                        pl5m1026(1:numel(tv)*2,end+1) = 1;
                        pl5b10(1:numel(tv)*2,end+1) = [flip(max(N10path')) min(N10path')]';
                        pl5b26(1:numel(tv)*2,end+1) = [flip(max(N26path')) min(N26path')]';
                    end;
                end;
                
                % fix min and max ice cover matrix
                if nhit>0&&(pl_range_depth==1||pl_range_buildup==1)&&(range1calc==1||range2calc==1);
                    if size(icemok,2)<2;
                        icemok(:,end+1) = icemok(:,1);
                    end;
                    plicemin(end+1,1:max(tv)+1) = interp1(tv,prod(icemok'),(0:1:max(tv)));
                    plicemax(end+1,1:max(tv)+1) = interp1(tv,double(sum(icemok')>0),(0:1:max(tv)));
                end;
                
                % fix banana plot matrix
                if min(find(looptest==0))==3 && pl_range_banana==1;
                    % calculate normalized P
                    Psp10 = Psp.sp10.*Pref10.*sample.thickSF.*sample.othercorr; % surface spal P10
                    Pmu10 = P_mu_z.mu10(1).*sample.othercorr; % surface muon P10
                    Psp26 = Psp.sp26.*Pref26.*sample.thickSF.*sample.othercorr; % surface spal P26
                    Pmu26 = P_mu_z.mu26(1).*sample.othercorr; % surface muon P26
                    P0avg10 = trapz(tv,(Psp10+Pmu10).*exp(-l10.*tv))/trapz(tv,exp(-l10.*tv));
                    P0avg26 = trapz(tv,(Psp26+Pmu26).*exp(-l26.*tv))/trapz(tv,exp(-l26.*tv));
                    % normalized N
                    N10n = sample.N10/P0avg10;
                    N26n = sample.N26/P0avg26;
                    delN10n = sample.delN10/P0avg10;
                    delN26n = sample.delN26/P0avg26;
                    if Punc == 1; % fix for prod rate uncertainty
                        if nhit > 0;
                            pl6del2N10y(end+1) = delN10n;
                            pl6del2N26y(end+1) = delN26n;
                        else;
                            pl6del2N10n(end+1) = delN10n;
                            pl6del2N26n(end+1) = delN26n;
                        end;
                        delN10n = sqrt(sample.delN10^2 - (sample.N10*delPref10/Pref10)^2) / P0avg10;
                        delN26n = sqrt(sample.delN26^2 - (sample.N26*delPref26/Pref26)^2) / P0avg26;
                    end;
                    if nhit > 0;
                        pl6delN10y(end+1) = delN10n;
                        pl6delN26y(end+1) = delN26n;
                        pl6N10y(end+1) = N10n;
                        pl6N26y(end+1) = N26n;
                    else;
                        pl6delN10n(end+1) = delN10n;
                        pl6delN26n(end+1) = delN26n;
                        pl6N10n(end+1) = N10n;
                        pl6N26n(end+1) = N26n;
                    end;
                    % calculate N paths
                    if nhit > 0; % if any hits
                        % pick out case with min and max ratio
                        [rmin rmini] = min(N26ok./N10ok);
                        [rmax rmaxi] = max(N26ok./N10ok);
                        P10r = P10ok(:,[rmini rmaxi]);
                        P26r = P26ok(:,[rmini rmaxi]);
                        % remove from Pok
                        P10ok(:,[rmini rmaxi]) = [];
                        P26ok(:,[rmini rmaxi]) = [];
                        % pick out up to ratiopathN random paths
                        rnum = min(ratiopathN-2,size(P10ok,2)); % number of further P to pick out
                        ridx = ceil(rand(1,rnum).*size(P10ok,2)); % idx to pick out
                        P10r(:,end+1:end+rnum) = P10ok(:,ridx);
                        P26r(:,end+1:end+rnum) = P26ok(:,ridx);
                        % calculate normalized N26/N10 ratios
                        Npath = Npathcalc(tv,P10r,l10,P26r,l26);
                        N10path = Npath.N10./P0avg10;
                        N26path = Npath.N26./P0avg26;
                        pl6N10path(1:size(N10path,1),end+1:end+size(N10path,2)) = N10path;
                        pl6N26path(1:size(N26path,1),end+1:end+size(N26path,2)) = N26path;
                    end;
                    % fix for Eline
                    plEl.thick_rho(end+1,1) = sample.thick * sample.rho;
                    plEl.atm(end+1,1) = sample.pressure;
                    plEl.RcEst(end+1,1) = LSDfix.RcEst;
                    plEl.othercorr(end+1,1) = sample.othercorr;
                    plEl.P10sp(end+1,1) = trapz(tv,Psp10.*exp(-l10.*tv))/trapz(tv,exp(-l10.*tv));
                    plEl.P26sp(end+1,1) = trapz(tv,Psp26.*exp(-l26.*tv))/trapz(tv,exp(-l26.*tv));
                    plEl.Lsp(end+1,1) = trapz(tv,sample.Lsp.*exp(-l10.*tv))/trapz(tv,exp(-l10.*tv));
                    plEl.P_mu10(end+1,1) = Pmu10;
                    plEl.P_mu26(end+1,1) = Pmu26;
                end;
                
                % reset parameters
                partest(1:6) = 0; randtest(1:6) = 0; % reset parameter tests
                looptest(min(find(looptest==0))) = 1; % reset looptest
                ddmok = []; glacEok = []; nonglacEok = []; LR04ok = []; Ndiffok = []; icemok = [];
                P10ok = []; P26ok = [];
                nhit = 0; % reset nhit
                randrun = 0; randhit = 0; randnum = 0; % reset random run parameters
                
                % reset input parameter limits
                glacEmin = glacEmin_in;
                glacEmax = glacEmax_in;
                nonglacEmin = nonglacEmin_in;
                nonglacEmax = nonglacEmax_in;
                LR04min = LR04min_in;
                LR04max = LR04max_in;
                
                % display for next while loop round
                if min(find(looptest==0)) == 2;
                    fprintf(1,'\n26Al calculating ranges');
                elseif min(find(looptest==0)) == 3;
                    fprintf(1,'\n10Be+26Al calculating ranges');
                end;
            end;
        end; % end while loop
    end; % end range calculation
    % ==============================================================================================
    
    fprintf(1,'\n');
    clear sample;
end; % end sample loop
% ==================================================================================================

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

% fix plotting =====================================================================================
fprintf(1,'plotting...');
% plot simple erosion sample depth
if numel(pl1d10)+numel(pl1d26) > 0;
    figure,set(gcf,'visible','off'); hold on; box on;
    leg = []; legin = {};
    % plot 10Be
    if numel(pl1d10) > 0;
        % fix for samples with varying tv
        pl1t10(pl1m10==0) = NaN; pl1d10(pl1m10==0) = NaN;
        % plot sample depth
        legtemp = plot(pl1t10,pl1d10,'color',color10);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{10}Be'};
        % plot uncertainties
        if numel(pl1d10uncmax) > 0;
            % fix for samples with varying tv
            pl1d10uncmin(pl1m10==0) = NaN; pl1d10uncmax(pl1m10==0) = NaN;
            % plot sample depth uncertainties
            plot(pl1t10,pl1d10uncmin,'linestyle','--','color',color10);
            plot(pl1t10,pl1d10uncmax,'linestyle','--','color',color10);
        end;
    end;
    % plot 26Al
    if numel(pl1d26) > 0;
        % fix for samples with varying tv
        pl1t26(pl1m26==0) = NaN; pl1d26(pl1m26==0) = NaN;
        % plot sample depth
        legtemp = plot(pl1t26,pl1d26,'color',color26);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{26}Al'};
        % plot uncertainties
        if numel(pl1d26uncmax) > 0;
            % fix for samples with varying tv
            pl1d26uncmin(pl1m26==0) = NaN; pl1d26uncmax(pl1m26==0) = NaN;
            % plot sample depth uncertainties
            plot(pl1t26,pl1d26uncmin,'linestyle','--','color',color26);
            plot(pl1t26,pl1d26uncmax,'linestyle','--','color',color26);
        end;
    end;
    % determine maxy for axis
    if max(pl1depth) > maxmaxy; maxy = maxmaxy; else; maxy = max(pl1depth); end;
    if maxy == 0; maxy = 1; end; % fix for zero erosion case
    if exist('absmaxy'); maxy = absmaxy; end;
    % plot ice cover
    if size(plice,1) > 1; plicev = (sum(plice==1) > 0); else; plicev = (plice == 1); end;
    leg(end+1) = plot_ice(plicev,maxy*0.97,[0.8 0.8 1]);
    legin(end+1) = {'ice cover'};
    % fix figure
    axis([0 pl_maxt/1E3 0 maxy]);
    axis('ij');
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('Sample depth (m)');
    legend(leg,legin,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off; set(gcf,'visible','on');
end;
% ========================================================================
% plot simple erosion nuclide buildup
if numel(pl2b10)+numel(pl2b26) > 0;
    figure,set(gcf,'visible','off'); hold on; box on;
    leg = []; legin = {};
    % plot 10Be
    if numel(pl2b10) > 0;
        % fix for samples with varying tv
        pl2t10(pl2m10==0) = NaN; pl2b10(pl2m10==0) = NaN;
        % plot buildup
        legtemp = plot(pl2t10,pl2b10,'color',color10);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{10}Be'};
        % plot uncertainties
        if numel(pl2b10uncmax) > 0;
            % fix for samples with varying tv
            pl2b10uncmin(pl2m10==0) = NaN; pl2b10uncmax(pl2m10==0) = NaN;
            plot(pl2t10,pl2b10uncmin,'linestyle','--','color',color10);
            plot(pl2t10,pl2b10uncmax,'linestyle','--','color',color10);
        end;
    end;
    % plot 26Al
    if numel(pl2b26) > 0;
        % fix for samples with varying tv
        pl2t26(pl2m26==0) = NaN; pl2b26(pl2m26==0) = NaN;
        % plot buildup
        legtemp = plot(pl2t26,pl2b26,'color',color26);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{26}Al'};
        % plot uncertainties
        if numel(pl2b26uncmax) > 0;
            % fix for samples with varying tv
            pl2b26uncmin(pl2m26==0) = NaN; pl2b26uncmax(pl2m26==0) = NaN;
            plot(pl2t26,pl2b26uncmin,'linestyle','--','color',color26);
            plot(pl2t26,pl2b26uncmax,'linestyle','--','color',color26);
        end;
    end;
    % plot ice cover
    if size(plice,1) > 1; plicev = (sum(plice==1) > 0); else; plicev = (plice == 1); end;
    leg(end+1) = plot_ice(plicev,3,[0.8 0.8 1]);
    legin(end+1) = {'ice cover'};
    % fix figure
    axis([0 pl_maxt/1E3 0 100]);
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('Cosmogenic concentration (%)');
    legend(leg,legin,'location','northwest');
    set(gca,'layer','top'); % plot axis on top
    hold off; set(gcf,'visible','on');
end;
% ========================================================================
% fix for ice cover plotting
if size(plicemin,1)>1; pliceminv = (prod(plicemin)==1); else; pliceminv = (plicemin==1); end;
if size(plicemax,1)>1; plicemaxv = (sum(plicemax==1)>0); else; plicemaxv = (plicemax==1); end;
% ========================================================================
% plot range erosion sample depth
if numel(pl3d10)+numel(pl3d26)+numel(pl3d1026) > 0;
    figure,set(gcf,'visible','off'); hold on; box on;
    leg = []; legin = {};
    % plot 10Be
    if numel(pl3d10) > 0;
        % fix for samples with varying tv
        pl3t10(pl3m10==0) = NaN; pl3d10(pl3m10==0) = NaN;
        % plot sample depth
        legtemp = plot(pl3t10,pl3d10,'color',color10);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{10}Be'};
    end;
    % plot 26Al
    if numel(pl3d26) > 0;
        % fix for samples with varying tv
        pl3t26(pl3m26==0) = NaN; pl3d26(pl3m26==0) = NaN;
        % plot sample depth
        legtemp = plot(pl3t26,pl3d26,'color',color26);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{26}Al'};
    end;
    % plot combined 10Be+26Al
    if numel(pl3d1026) > 0;
        % fix for samples with varying tv
        pl3t1026(pl3m1026==0) = NaN; pl3d1026(pl3m1026==0) = NaN;
        % plot sample depth
        legtemp = plot(pl3t1026,pl3d1026,'color',color1026);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{10}Be+^{26}Al'};
    end;
    % determine maxy for axis
    if max(pl3depth) > maxmaxy; maxy = maxmaxy; else; maxy = max(pl3depth); end;
    if exist('absmaxy'); maxy = absmaxy; end;
    % plot ice cover
    leg(end+1) = plot_ice(pliceminv,maxy*0.95,[0.6 0.6 1]);
    legin(end+1) = {'min ice cover'};
    leg(end+1) = plot_ice(plicemaxv,maxy*0.97,[0.2 0.2 1]);
    legin(end+1) = {'max ice cover'};
    % fix figure
    axis([0 pl_maxt/1E3 0 maxy]);
    axis('ij');
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('Sample depth (m)');
    legend(leg,legin,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off; set(gcf,'visible','on');
end;
% ========================================================================
% plot range single nuclide buildup
if numel(pl4b10)+numel(pl4b26) > 0;
    figure,set(gcf,'visible','off'); hold on; box on;
    leg = []; legin = {};
    % plot 10Be
    if numel(pl4b10) > 0;
        % fix for samples with varying tv
        pl4t10(pl4m10==0) = NaN; pl4b10(pl4m10==0) = NaN;
        % plot buildup
        legtemp = plot(pl4t10,pl4b10,'color',color10);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{10}Be'};
    end;
    % plot 26Al
    if numel(pl4b26) > 0;
        % fix for samples with varying tv
        pl4t26(pl4m26==0) = NaN; pl4b26(pl4m26==0) = NaN;
        % plot buildup
        legtemp = plot(pl4t26,pl4b26,'color',color26);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{26}Al'};
    end;
    % plot ice cover
    leg(end+1) = plot_ice(pliceminv,5,[0.6 0.6 1]);
    legin(end+1) = {'min ice cover'};
    leg(end+1) = plot_ice(plicemaxv,3,[0.2 0.2 1]);
    legin(end+1) = {'max ice cover'};
    % fix figure
    axis([0 pl_maxt/1E3 0 100]);
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('Nuclide concentration (%) - single nuclide simulation');
    legend(leg,legin,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off; set(gcf,'visible','on');
end;
% ========================================================================
% plot range combined nuclide buildup
if numel(pl5t1026) > 0;
    figure,set(gcf,'visible','off'); hold on; box on;
    leg = []; legin = {};
    % fix for samples with varying tv
    pl5t1026(pl5m1026==0) = NaN;
    % plot 10Be
    if numel(pl5b10) > 0;
        % fix for samples with varying tv
        pl5b10(pl5m1026==0) = NaN;
        % plot buildup
        legtemp = plot(pl5t1026,pl5b10,'color',color10);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{10}Be'};
    end;
    % plot 26Al
    if numel(pl5b26) > 0;
        % fix for samples with varying tv
        pl5b26(pl5m1026==0) = NaN;
        % plot buildup
        legtemp = plot(pl5t1026,pl5b26,'color',color26);
        leg(end+1) = legtemp(1); legin(end+1) = {'^{26}Al'};
    end;
    % plot ice cover
    leg(end+1) = plot_ice(pliceminv,5,[0.6 0.6 1]);
    legin(end+1) = {'min ice cover'};
    leg(end+1) = plot_ice(plicemaxv,3,[0.2 0.2 1]);
    legin(end+1) = {'max ice cover'};
    % fix figure
    axis([0 pl_maxt/1E3 0 100]);
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('Nuclide concentration (%) - combined nuclide simulation');
    legend(leg,legin,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off; set(gcf,'visible','on');
end;
% ========================================================================
% plot banana
if numel(pl6N10y)+numel(pl6N10n) > 0;
    figure,set(gcf,'visible','off'); hold on; box on;
    leg = []; legin = {};
    color_nohit = [0.7 0.7 0.7];
    % plot points
    semilogx(pl6N10y,pl6N26y./pl6N10y,'.','color',color1026,'markersize',15);
    semilogx(pl6N10n,pl6N26n./pl6N10n,'.','color',color_nohit,'markersize',15);
    % calculate and plot sigma lines
    for i = 1:numel(pl6N10y);
        sigmal = sigmaline(pl6N10y(i),pl6N26y(i),pl6delN10y(i),pl6delN26y(i));
        % plot sigmaline
        semilogx(sigmal.x,sigmal.y,'color',color1026);
        if Punc == 1;
            sigmal = sigmaline(pl6N10y(i),pl6N26y(i),pl6del2N10y(i),pl6del2N26y(i));
            semilogx(sigmal.x,sigmal.y,'color',color1026,'linestyle',':');
        end;
    end;
    for i = 1:numel(pl6N10n);
        sigmal = sigmaline(pl6N10n(i),pl6N26n(i),pl6delN10n(i),pl6delN26n(i));
        % plot sigmaline
        semilogx(sigmal.x,sigmal.y,'color',color_nohit);
        if Punc == 1;
            sigmal = sigmaline(pl6N10n(i),pl6N26n(i),pl6del2N10n(i),pl6del2N26n(i));
            semilogx(sigmal.x,sigmal.y,'color',color_nohit,'linestyle',':');
        end;
    end;
    % fix for N10 = 0 case
    pl6N10path(pl6N10path==0) = NaN; pl6N26path(pl6N10path==0) = NaN;
    % plot 26/10 ratios against N10
    semilogx(pl6N10path,pl6N26path./pl6N10path,'color',color1026);
    % create data for the simple-exposure line including ratio uncertainties
    tempt = logspace(0,7,100);
    be = (1/consts.l10)*(1-exp(-consts.l10*tempt));
    al = (1/consts.l26)*(1-exp(-consts.l26*tempt));
    al_be = al./be;
    % calculate Eline
    fprintf(1,' calculating erosion end-point line...');
    [bee,ale_bee] = Eline(plEl,consts);
    % plot simple exposure line and erosion end-point line
    semilogx(be,al_be,'color','black');
    semilogx(bee,ale_bee,'linestyle','--','color','black');
    axis([1000 3000000 0.2 1.2]);
    xlabel('[^{10}Be]*');
    ylabel('[^{26}Al]*/[^{10}Be]*');
    set(gca,'layer','top'); % plot axis on top
    set(gca,'XScale','log'); % fix for matlab
    hold off; set(gcf,'visible','on');
end;
fprintf(1,' done!\n');
toc()
% end glacialE function ============================================================================


% subfunction Ncalc ================================================================================
function Nend = Ncalc(fulldm,noicem,tv,Psp,Pref,P_mu_z,l,sample);
    % Calculate N(t) including decay and erosion
    Pspv = Psp.*Pref.*sample.thickSF.*sample.othercorr; % surface spal prod
    Pspm = repmat(Pspv',1,size(fulldm,2)) .* noicem; % spal prod matrix
    dpfs = exp(-bsxfun(@rdivide,fulldm,sample.Lsp')); % spal depth dependence matrix
    dcf = repmat(exp(-tv.*l)',1,size(fulldm,2)); % decay factor;
    Pmum = interp1(sample.z_mu',P_mu_z',fulldm,'pchip') .* sample.othercorr .* noicem; % muon P matr
    Pmum(isnan(Pmum)) = 0; % set muon prod to 0 for depths > 2E5
    Nend = trapz(tv',(Pspm.*dcf.*dpfs + Pmum.*dcf)); % potential N back in time
% end subfunction Ncalc ============================================================================


% subfunction nucl1_getP ===========================================================================
function Pfull = nucl1_getP(Eglacv,fulldm,Eout,sample,Psp0,Pref,P_mu_z,noicev);
    % calculate P through time
    fulldv = interp1(Eglacv',fulldm',Eout.erosion,'pchip'); % interpolate depth (g/cm2)
    dpfs = exp(-fulldv./sample.Lsp); % spallation depth dependence
    Pspv = Psp0 .* Pref .* sample.thickSF .* sample.othercorr .* dpfs .* noicev; % spallation prod
    Pmuv = interp1(sample.z_mu,P_mu_z,fulldv,'pchip') .* sample.othercorr .* noicev;
    Pfull = Pspv' + Pmuv'; % full prod
    % same for min and max uncertainties
    if isfield(Eout,'mindelE');
        mindv = interp1(Eglacv',fulldm',Eout.erosion-Eout.mindelE,'pchip'); % interp depth (g/cm2)
        mindpfs = exp(-mindv./sample.Lsp); % spallation depth dependence
        minPspv = Psp0.*Pref.*sample.thickSF.*sample.othercorr.*mindpfs.*noicev; % spallation prod
        minPmuv = interp1(sample.z_mu,P_mu_z,mindv,'pchip') .* sample.othercorr .* noicev;
        Pfull(:,2) = minPspv' + minPmuv'; % full prod
        maxdv = interp1(Eglacv',fulldm',Eout.erosion+Eout.maxdelE,'pchip'); % (g/cm2)
        maxdpfs = exp(-maxdv./sample.Lsp); % spallation depth dependence
        maxPspv = Psp0.*Pref.*sample.thickSF.*sample.othercorr.*maxdpfs.*noicev; % spallation prod
        maxPmuv = interp1(sample.z_mu,P_mu_z,maxdv,'pchip') .* sample.othercorr .* noicev;
        Pfull(:,3) = maxPspv' + maxPmuv'; % full prod
    end;
% end subfunction nucl1_getP =======================================================================


% subfunction Pcalc1026 ============================================================================
function out = Pcalc1026(sample,fulldm,noicem,tv,Psp10,Pref10,P_mu_z10,l10,Psp26,Pref26,P_mu_z26,...
    l26);
    % fix if only one nuclide input
    if nargin == 8;
        Psp26 = Psp10; Pref26 = Pref10; P_mu_z26 = P_mu_z10; l26 = l10;
    end;
    
    % spal depth dependence matrix
    dpfs = exp(-bsxfun(@rdivide,fulldm,sample.Lsp'));
    
    % Calculate P10(t) including erosion
    Pspv10 = Psp10.*Pref10.*sample.thickSF.*sample.othercorr; % surface spal prod
    Pspm10 = repmat(Pspv10',1,size(fulldm,2)).*noicem.*dpfs; % spal prod matrix
    Pmum10 = interp1(sample.z_mu',P_mu_z10',fulldm,'pchip').*sample.othercorr.*noicem; % muon P matr
    Pmum10(isnan(Pmum10)) = 0; % set muon prod to 0 for depths > 2E5
    out.P10 = Pspm10 + Pmum10;
    
    % Calculate P26(t) including erosion
    Pspv26 = Psp26.*Pref26.*sample.thickSF.*sample.othercorr; % surface spal prod
    Pspm26 = repmat(Pspv26',1,size(fulldm,2)).*noicem.*dpfs; % spal prod matrix
    Pmum26 = interp1(sample.z_mu',P_mu_z26',fulldm,'pchip').*sample.othercorr.*noicem; % muon P matr
    Pmum26(isnan(Pmum26)) = 0; % set muon prod to 0 for depths > 2E5
    out.P26 = Pspm26 + Pmum26;
% end subfunction Pcalc1026 ========================================================================


% subfunction Npathcalc ============================================================================
function out = Npathcalc(tv,P10,l10,P26,l26);
    % fix if only one nuclide input
    if nargin == 3;
        P26 = P10; l26 = l10;
    end;
    
    % fix flipped full prod matrices and tv matrix
    P10m = flip(P10);
    P26m = flip(P26);
    tvflip = flip(tv);
    
    Pmean10 = (P10m(1:end-1,:)+P10m(2:end,:))./2;
    Pmean26 = (P26m(1:end-1,:)+P26m(2:end,:))./2;
    dtv = tvflip(1:end-1) - tvflip(2:end);
    fidx = min(find(dtv(1:end-1)-dtv(2:end)==0));
    
    N10m(1,size(P10m,2)) = 0;
    N26m(1,size(P26m,2)) = 0;
    
    % for loop used for first part with varying time steps
    if fidx > 1;
        for i = 1:fidx-1;
            N10m(i+1,:) = (N10m(end,:) + Pmean10(i,:).*(tvflip(i)-tvflip(i+1))) .* ...
                exp(-(tvflip(i)-tvflip(i+1)).*l10);
            N26m(i+1,:) = (N26m(end,:) + Pmean26(i,:).*(tvflip(i)-tvflip(i+1))) .* ...
                exp(-(tvflip(i)-tvflip(i+1)).*l26);
        end;
        % remove calculated steps
        tvflip(1:fidx-1) = [];
        dtv(1:fidx-1) = [];
        Pmean10(1:fidx-1,:) = [];
        Pmean26(1:fidx-1,:) = [];
    end;
    
    % while loop used for filtering parts with constant time step
    while numel(tvflip) > 1;
        widx = max(find(dtv==dtv(1)));
        Ptemp10 = Pmean10(1:widx,:).*dtv(1);
        Ptemp10(2:end+1,:) = Ptemp10;
        Ptemp10(1,:) = N10m(end,:);
        decay10 = -exp(-l10.*dtv(1));
        Ntemp10 = filter(1,[1,decay10],Ptemp10);
        N10m(end+1:end+widx,:) = Ntemp10(2:end,:);
        Ptemp26 = Pmean26(1:widx,:).*dtv(1);
        Ptemp26(2:end+1,:) = Ptemp26;
        Ptemp26(1,:) = N26m(end,:);
        decay26 = -exp(-l26.*dtv(1));
        Ntemp26 = filter(1,[1,decay26],Ptemp26);
        N26m(end+1:end+widx,:) = Ntemp26(2:end,:);
        % remove calculated steps
        tvflip(1:widx) = [];
        dtv(1:widx) = [];
        Pmean10(1:widx,:) = [];
        Pmean26(1:widx,:) = [];
    end;
    
    % flip back
    out.N10 = flip(N10m);
    out.N26 = flip(N26m);
% end subfunction Npathcalc ========================================================================


% subfunction fix_output1 ==========================================================================
function results = fix_output1(N,delN,Pref,delPref,Nend,Eglacv,glacErate,nuclstr);
    % fix unit
    if glacErate == 1;
        unitstr = 'mm/ka';
    else;
        unitstr = 'cm/glac';
    end;
    
    if N >= min(Nend) && N <= max(Nend); % if N is within range
        % Interpolate erosion
        erosion = interp1(Nend,Eglacv,N,'pchip');
        
        % uncertainty estimation: interpolate max and min erosion based on concentration unc
        if N-delN >= min(Nend);
            maxdel1 = interp1(Nend,Eglacv,N-delN,'pchip');
            maxdelE = sqrt((maxdel1-erosion)^2 + (erosion.*delPref./Pref)^2);
        else;
            maxdelE = max(Eglacv); % set uncertainty to max Eglacv
        end;
        if N+delN <= max(Nend);
            mindel1 = interp1(Nend,Eglacv,N+delN,'pchip');
            mindelE = sqrt((erosion-mindel1)^2 + (erosion.*delPref./Pref)^2);
        else;
            mindelE = erosion; % set uncertainty to erosion
        end;
        
        % fix outstr
        results.outstr(1) = {num2str(erosion,'%.1f')};
        results.outstr(2) = {num2str(maxdelE,'%.1f')};
        results.outstr(3) = {num2str(mindelE,'%.1f')};
        
        % fix erosion and uncertainty
        results.erosion = erosion;
        results.maxdelE = maxdelE;
        results.mindelE = mindelE;
        
        % display output
        fprintf(1,' \t%s = %s ± %s/%s %s',nuclstr,results.outstr{1},results.outstr{2},...
            results.outstr{3},unitstr);
    else; % if N is not within range
        if N < min(Nend); % if too low N
            if N+delN < min(Nend); % if no overlap
                results.outstr(1) = {strcat('>',num2str(max(Eglacv),'%.0f'))};
                fprintf(1,' \t%s: %s %s (no prior exposure!)',nuclstr,results.outstr{1},unitstr);
            else; % if one-sided overlap with Eglacv
                erosion = interp1(Nend,Eglacv,N+delN,'pchip');
                results.outstr(1) = {strcat('>',num2str(erosion,'%.1f'))};
                fprintf(1,' \t%s: %s %s (not enough prior exposure!)',nuclstr,results.outstr{1},...
                    unitstr);
            end;
            erosion = max(Eglacv); % set erosion to max Eglacv
        else; % if too high N
            if N-delN > max(Nend); % if no overlap
                results.outstr(1) = {'-'};
                fprintf(1,' \t%s: too much prior exposure!',nuclstr);
            else; % if one-sided prior exposure
                erosion = interp1(Nend,Eglacv,N-delN,'pchip');
                results.outstr(1) = {strcat('<',num2str(erosion,'%.1f'))};
                fprintf(1,' \t%s: %s %s (too much prior exposure!)',nuclstr,results.outstr{1},...
                    unitstr);
            end;
            erosion = min(Eglacv); % set erosion to min Eglacv
        end;
        % fix erosion
        results.erosion = erosion;
    end;
% end subfunction fix_output1 ======================================================================


% subfunction sigmaline ============================================================================
function out = sigmaline(N10n,N26n,delN10n,delN26n);
    % Estimate range and create mesh
    R = (N26n/N10n);
    delR = sqrt((delN26n/N26n)^2 + (delN10n/N10n)^2);
    
    [x,y] = meshgrid((N10n-4*delN10n):(0.1*delN10n):(N10n+4*delN10n),...
        (R*(1-4*delR)):(0.1*R*delR):(R*(1+4*delR)));
    
    % calculate PDF
    Prob = x.*exp(-0.5.*((((y.*x) - N26n)./delN26n).^2 + ((x - N10n)./delN10n).^2));
    
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
    
    % Sometimes contourc returns several contours for one level - grid size issue?
    % This is spurious, so plot only the major one.
    contourStarts = find(cmat(1,:) == sigma1);
    contourSizes = cmat(2,contourStarts);
    contourToPlot = find(contourSizes == max(contourSizes));
    
    out.x = cmat(1,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + ...
        contourSizes(contourToPlot)));
    out.y = cmat(2,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + ...
        contourSizes(contourToPlot)));
% end subfunction sigmaline ========================================================================


% subfunction Eline ================================================================================
function [x,y] = Eline(plEl,consts);
    % Precompute P_mu(z) to ~200,000 g/cm2 - start at the average sample mid-depth.
    z_sp = [0 logspace(0,5.3,100)];
    z_mu = z_sp+(mean(plEl.thick_rho)./2);
    P_mu_z = P_mu_expage(z_mu,mean(plEl.atm),mean(plEl.RcEst),consts.SPhiInf,1,1,consts,'no');
    P_mu_z10 = P_mu_z.mu10.*mean(plEl.othercorr);
    P_mu_z26 = P_mu_z.mu26.*mean(plEl.othercorr);
    P10sp_z = mean(plEl.P10sp).*exp(-z_sp./mean(plEl.Lsp));
    P26sp_z = mean(plEl.P26sp).*exp(-z_sp./mean(plEl.Lsp));
    
    tempe = logspace(-5,1,100);
    for j = 1:numel(tempe);
        tv = z_sp./tempe(j);
        dcf10 = exp(-tv.*consts.l10);
        dcf26 = exp(-tv.*consts.l26);
        N10tv = cumtrapz(tv,(P10sp_z.*dcf10 + P_mu_z10.*dcf10));
        N26tv = cumtrapz(tv,(P26sp_z.*dcf26 + P_mu_z26.*dcf26));
        bee(j+1,1) = N10tv(end)./(mean(plEl.P10sp) + mean(plEl.P_mu10));
        ale(j+1,1) = N26tv(end)./(mean(plEl.P26sp) + mean(plEl.P_mu26));
    end;
    bee(1,1) = (1/consts.l10)*(1-exp(-consts.l10*1E7));
    ale(1,1) = (1/consts.l26)*(1-exp(-consts.l26*1E7));
    ale_bee = ale./bee;
    x = bee;
    y = ale_bee;
% end subfunction Eline ============================================================================


% subfunction plot_ice =============================================================================
function legh = plot_ice(plicev,y,clr);
    iceend = (filter([1 2],1,plicev) == 1); % find end of ice cover periods
    icestart = flip(filter([1 2],1,flip(plicev)) == 1); % find start of ice cover periods
    iceendidx = find(iceend == 1)-1; % time of iceend
    icestartidx = find(icestart == 1)-1; % time of icestart
    x = [];
    if numel(iceendidx) >= 1;
        legh = plot([iceendidx(1) icestartidx(1)+1]./1E3,[y y],'color',clr,'LineWidth',2);
    end;
    if numel(iceendidx) >= 2;
        xm = [iceendidx(2:end);icestartidx(2:end)]./1E3; ym = repmat(y,size(xm));
        plot(xm,ym,'color',clr,'LineWidth',2);
    end;
% end subfunction plot_ice =========================================================================
