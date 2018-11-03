function glacialE();

% 10Be and 26Al glacial erosion calculator.
% Function for quantification/investigation of glacial erosion based on 10Be and/or 26Al conc.
% Long and complicated code that can do some interesting calculations and plots.
% Interpolates the best fitting glacial erosion to yield the 10Be and/or 26Al conc given as input
% with a max duration (mt) of exposure, production during ice-free periods (determined by
% glacialE_LR04.txt) and shielding from non-glacial (given in input) and glacial (output) erosion.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2018 (jakob.heyman@gu.se)

clear;
close all;
tic();

% What version is this?
ver = '201810';

% =============== MAKE CHOICES HERE ================================================================
% max time (determines the start of the simulation - max 1E7)
mt = 2588000;

% use glacial erosion rate (1) or incremental erosion depth steps (0)?
glacErate = 1;

% which calculation(s) to do? 1 = yes, 0 = no
simpleEcalc = 1;  % simple glacial erosion rate calculation (single nuclides)
range1calc = 0; % calculate possible ranges for parameters (single nuclides)
range2calc = 0; % calculate possible ranges for parameters (10Be + 26Al)

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
glacEmax = 1000; % mm/ka or cm/glac - max value
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
parn = 10;

% add production rate uncertainty? 1 = yes
Punc = 1;

% plotting choices ===================================================
color10 = 'red'; % color for 10Be
color26 = 'blue'; % color for 26Al

% plot choices
pl_simpleEcalc_buildup = 1; % plot nuclide build-up for simple glacial erosion (simpleEcalc = 1)
pl_simpleEcalc_depth = 1;   % plot sample depth for simple glacial erosion (simpleEcalc = 1)
pl_range_depth = 1;   % plot sample depth for range calculation(s) (rangeXnulc = 1)
pl_range_banana = 1;  % plot 26/10 banana for range calculation(s) (range2nulc = 1)
ratiopathN = 100; % number of OK sample 26Al/10Be ratio paths to plot in banana

% max time for plotting (yr)
pl_maxt = 1E6;
if pl_maxt > mt; pl_maxt = mt; end; % don't allow too large maxt
% max depth for plotting (m)
maxmaxy = 10;
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
pl(1:6) = 0;
if pl_simpleEcalc_buildup==1 && simpleEcalc==1; pl(1) = max(pl)+1; end;
if pl_simpleEcalc_depth==1 && simpleEcalc==1; pl(2) = max(pl)+1; end;
if pl_range_depth==1 && range1calc==1 && sum(samplein.N10)>0; pl(3) = max(pl)+1; end;
if pl_range_depth==1 && range1calc==1 && sum(samplein.N26)>0; pl(4) = max(pl)+1; end;
if pl_range_depth==1 && range2calc==1 && sum(samplein.N10.*samplein.N26)>0; pl(5) = max(pl)+1; end;
if pl_range_banana==1 && range2calc==1 && sum(samplein.N10.*samplein.N26)>0; pl(6) = max(pl)+1; end;
pl10b = 1; pl10d = 1; pl26b = 1; pl26d = 1;
pldepth2 = []; pldepth3 = []; pldepth4 = []; pldepth5 = [];
plicemin3 = []; plicemax3 = []; plicemin4 = []; plicemax4 = []; plicemin5 = []; plicemax5 = [];
legh1 = []; legin1 = {}; legh2 = []; legin2 = {};
plr = [];
plice = [];
plEl.thick_rho = []; plEl.atm = []; plEl.RcEst = []; plEl.othercorr = [];
plEl.P10sp = []; plEl.P26sp = []; plEl.Lsp = []; plEl.P_mu10 = []; plEl.P_mu26 = [];
clrv = ...
    [0 0 0;1 0 0;0 0 1;0 1 0;1 1 0;0.75 0 0.75;1 0.4 0;0.5 0 0;0 0 0.5;0 0.5 0;1 0.5 0.5;0.4 0.4 1];

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
    
    % pick color for plotting
    clr = clrv(i-(ceil(i/size(clrv,1))-1)*size(clrv,1),:);
    
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
            
            % plot nuclide buildup
            if pl(1) > 0;
                if isfield(out10,'mindelE');
                    N10v = nucl1_buildup(Eglacv,fulldm,out10.erosion,sample,Psp.sp10,P_mu_z.mu10,...
                        noicev,tv,l10,out10.mindelE,out10.maxdelE);
                else;
                    N10v = nucl1_buildup(Eglacv,fulldm,out10.erosion,sample,Psp.sp10,P_mu_z.mu10,...
                        noicev,tv,l10);
                end;
                
                figure(pl(1)),set(pl(1),'visible','off');
                hold on;
                if pl10b == 0;
                    plot(tv'./1E3,N10v.mid,'color',color10);
                else;
                    legh1(end+1) = plot(tv'./1E3,N10v.mid,'color',color10);
                    legin1(end+1) = {'^{10}Be'};
                    pl10b = 0;
                end;
                % plot uncertainties
                if simpleEcalcunc==1 && isfield(N10v,'unc');
                    plot(tv'./1E3,N10v.unc,'color',color10,'linestyle','--');
                end;
            end;
            
            % calculate depth history
            dv10 = interp1(Eglacv',ddm',out10.erosion,'pchip')./1E2; % interpolated depth (m)
            if max(tv) >= 1E5;
                output(i+1,outn(3)) = {num2str(interp1(tv,dv10,1E5,'pchip'),'%.2f')}; % d(m) 100 ka
            end;
            if max(tv) >= 1E6;
                output(i+1,outn(4)) = {num2str(interp1(tv,dv10,1E6,'pchip'),'%.2f')}; % d(m) 1 Ma
            end;
            
            % plot depth history
            if pl(2) > 0;
                pldepth2(end+1) = interp1(tv,dv10,min(pl_maxt,max(tv)),'pchip'); % add d for y axis
                figure(pl(2)),set(pl(2),'visible','off');
                hold on;
                if pl10d == 0;
                    plot(tv./1E3,dv10,'color',color10);
                else;
                    legh2(end+1) = plot(tv./1E3,dv10,'color',color10);
                    legin2(end+1) = {'^{10}Be'};
                    pl10d = 0;
                end;
                % plot uncertainties
                if simpleEcalcunc==1 && isfield(out10,'mindelE'); % interpret max depth (m)
                    pluncmin10 = interp1(Eglacv',ddm',out10.erosion-out10.mindelE,'pchip')./1E2;
                    pluncmax10 = interp1(Eglacv',ddm',out10.erosion+out10.maxdelE,'pchip')./1E2;
                    pldepth2(end+1) = interp1(tv,pluncmax10,pl_maxt,'pchip'); % max d for y-axis
                    plot(tv./1E3,pluncmin10,'linestyle','--','color',color10); % plot min unc
                    plot(tv./1E3,pluncmax10,'linestyle','--','color',color10); % plot max unc
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
            
            % plot nuclide buildup
            if pl(1) > 0;
                if isfield(out26,'mindelE');
                    N26v = nucl1_buildup(Eglacv,fulldm,out26.erosion,sample,Psp.sp26,P_mu_z.mu26,...
                        noicev,tv,l26,out26.mindelE,out26.maxdelE);
                else;
                    N26v = nucl1_buildup(Eglacv,fulldm,out26.erosion,sample,Psp.sp26,P_mu_z.mu26,...
                        noicev,tv,l26);
                end;
                
                figure(pl(1)),set(pl(1),'visible','off');
                hold on;
                if pl26b == 0;
                    plot(tv'./1E3,N26v.mid,'color',color26);
                else;
                    legh1(end+1) = plot(tv'./1E3,N26v.mid,'color',color26);
                    legin1(end+1) = {'^{26}Be'};
                    pl26b = 0;
                end;
                % plot uncertainties
                if simpleEcalcunc==1 && isfield(N26v,'unc');
                    plot(tv'./1E3,N26v.unc,'color',color26,'linestyle','--');
                end;
            end;
            
            %{
            % plot nuclide buildup
            if pl(1) > 0;
                N26v = nucl1_buildup(Eglacv,fulldm,out26.erosion,sample,Psp.sp26,P_mu_z.mu26,...
                    noicev,tv,l26);
                
                figure(pl(1)),set(pl(1),'visible','off');
                hold on;
                if pl26b == 0;
                    plot(tv./1E3,N26v,'color',color26);
                else;
                    legh1(end+1) = plot(tv./1E3,N26v,'color',color26);
                    legin1(end+1) = {'^{26}Al'};
                    pl26b = 0;
                end;
            end;
            %}
            
            % calculate depth history
            dv26 = interp1(Eglacv',ddm',out26.erosion,'pchip')./1E2; % interpolated depth (m)
            if max(tv) >= 1E5;
                output(i+1,outn(7)) = {num2str(interp1(tv,dv26,1E5,'pchip'),'%.2f')}; % d(m) 100 ka
            end;
            if max(tv) >= 1E6;
                output(i+1,outn(8)) = {num2str(interp1(tv,dv26,1E6,'pchip'),'%.2f')}; % d(m) 1 Ma
            end;
            
            % plot depth history
            if pl(2) > 0;
                pldepth2(end+1) = interp1(tv,dv26,min(pl_maxt,max(tv)),'pchip'); % add d for y axis
                figure(pl(2)),set(pl(2),'visible','off');
                hold on;
                if pl26d == 0;
                    plot(tv./1E3,dv26,'color',color26);
                else;
                    legh2(end+1) = plot(tv./1E3,dv26,'color',color26);
                    legin2(end+1) = {'^{26}Al'};
                    pl26d = 0;
                end;
                % plot uncertainties
                if simpleEcalcunc==1 && isfield(out26,'mindelE'); % interpret max depth (m)
                    pluncmin26 = interp1(Eglacv',ddm',out26.erosion-out26.mindelE,'pchip')./1E2;
                    pluncmax26 = interp1(Eglacv',ddm',out26.erosion+out26.maxdelE,'pchip')./1E2;
                    pldepth2(end+1) = interp1(tv,pluncmax26,pl_maxt,'pchip'); % max d for y-axis
                    plot(tv./1E3,pluncmin26,'linestyle','--','color',color26); % plot min unc
                    plot(tv./1E3,pluncmax26,'linestyle','--','color',color26); % plot max unc
                end;
            end;
        end; % end 26Al calculations
    end; % end simpleEcalc
    % ==============================================================================================
    
    % if calculating ranges ========================================================================
    if range1calc==1 || (range2calc==1 && nucl10==1 && nucl26==1);      
        % matrix/vectors to be filled in while loop
        ddmok = []; glacEok = []; nonglacEok = []; LR04ok = []; Ndiffok = [];
        P10ok = []; P26ok = []; N10ok = []; N26ok = []; nhit = 0;
        
        % parameter test
        partest(1:6) = 0;
        
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
            
            % fix parameter vectors
            glacEv1 = repmat(linspace(glacEmin,glacEmax,glacEn),LR04n*nonglacEn,1);
            glacEv = glacEv1(:)';
            nonglacEv = repmat(repmat(linspace(nonglacEmin,nonglacEmax,nonglacEn),1,LR04n),1,...
                glacEn);
            LR04v1 = repmat(linspace(LR04min,LR04max,LR04n),nonglacEn,1);
            LR04v = repmat(LR04v1(:)',1,glacEn);
            
            % calculate parameter steps
            glacEst = (glacEmax-glacEmin)/max(glacEn-1,1);
            nonglacEst = (nonglacEmax-nonglacEmin)/max(nonglacEn-1,1);
            LR04st = (LR04max-LR04min)/max(LR04n-1,1);
            
            % ice cover matrix
            icem = (repmat(d18O',1,simn) >= LR04v);
            
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
                
                % pick out OK depth vectors from ddm
                ddmok(:,end+1:end+numel(hitidx)) = ddm(:,hitidx);
            elseif looptest(2) == 0; % 26Al single nuclide range
                % calculate N(t) for matrix
                Nend = Ncalc(fulldm,noicem,tv,Psp.sp26,Pref26,P_mu_z.mu26,l26,sample);
                
                % pick out concentrations within measured uncertainty
                Ndiff = abs(Nend-sample.N26);
                hitidx = find(Ndiff <= sample.delN26);
                
                % pick out OK depth vectors from ddm
                ddmok(:,end+1:end+numel(hitidx)) = ddm(:,hitidx);
            elseif looptest(3) == 0; % 26Al/10Be nuclide range
                % calculate N(t) for matrix
                Nend10 = Ncalc(fulldm,noicem,tv,Psp.sp10,Pref10,P_mu_z.mu10,l10,sample);
                Nend26 = Ncalc(fulldm,noicem,tv,Psp.sp26,Pref26,P_mu_z.mu26,l26,sample);
                
                % calculate P10 and P26 for matrix
                P1026m = Pcalc1026(fulldm,noicem,tv,Psp.sp10,Psp.sp26,Pref10,Pref26,P_mu_z.mu10,...
                    P_mu_z.mu26,l10,l26,sample);
                
                % pick out concentrations within measured uncertainty
                Ndiff = abs(Nend10-sample.N10)./sample.delN10 + ...
                    abs(Nend26-sample.N26)./sample.delN26;
                hitidx = find((abs(Nend10-sample.N10) <= sample.delN10) .* ...
                    (abs(Nend26-sample.N26) <= sample.delN26) == 1);
                
                % pick out concentrations (for ratio path plot selection)
                N10ok(end+1:end+numel(hitidx)) = Nend10(hitidx);
                N26ok(end+1:end+numel(hitidx)) = Nend26(hitidx);
                
                % pick out OK depth vectors from ddm
                ddmok(:,end+1:end+numel(hitidx)) = ddm(:,hitidx);
                
                % pick out OK P10 and P26 vectors
                P10ok(:,end+1:end+numel(hitidx)) = P1026m.P10(:,hitidx);
                P26ok(:,end+1:end+numel(hitidx)) = P1026m.P26(:,hitidx);
            end;
            
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
            if numel(hitidx) > 0;
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
                
                % test if the actual parameter has been found
                parcltest = 0;
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
                
                % test if parameters min/max limits matches input min/max
                if partest(1)==0 && min(glacEok)==glacEmin_in;
                    if min(find(partest==0)) == 1; parcltest = 1; end;
                    partest(1) = 1;
                end;
                if partest(2)==0 && max(glacEok)==glacEmax_in;
                    if min(find(partest==0)) == 2; parcltest = 1; end;
                    partest(2) = 1;
                end;
                if partest(3)==0 && min(nonglacEok)==nonglacEmin_in;
                    if min(find(partest==0)) == 3; parcltest = 1; end;
                    partest(3) = 1;
                end;
                if partest(4)==0 && max(nonglacEok)==nonglacEmax_in;
                    if min(find(partest==0)) == 4; parcltest = 1; end;
                    partest(4) = 1;
                end;
                if partest(5)==0 && min(LR04ok)==LR04min_in;
                    if min(find(partest==0)) == 5; parcltest = 1; end;
                    partest(5) = 1;
                end;
                if partest(6)==0 && max(LR04ok)==LR04max_in;
                    if min(find(partest==0)) == 6; parcltest = 1; end;
                    partest(6) = 1;
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
                    %elseif min(partest==0) <= 2;  % if glacE is the actual parameter: center limits
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
                    %elseif min(partest==0)==3 || min(partest==0)==4; % if nonglacE is actual...
                    elseif nonglacEmax-nonglacEmin < nonglacEmax_in-nonglacEmin_in;
                        nonglacEmin = max(nonglacEok(idxb)-nonglacEst*(glacEn-1)/2,nonglacEmin_in);
                        nonglacEmax = min(nonglacEok(idxb)+nonglacEst*(glacEn-1)/2,nonglacEmax_in);
                    end;
                    
                    % fix new limits for LR04
                    if paridx == 3; % if LR04 has the largest relative step
                        LR04min = max(LR04ok(idxb)-LR04st*(LR04n-1)/4,LR04min_in);
                        LR04max = min(LR04ok(idxb)+LR04st*(LR04n-1)/4,LR04max_in);
                    %elseif min(partest==0) >= 5; % if LR04 is the actual parameter: center limits
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
            
            % if done with parameter optimization
            if  prod(partest) == 1;
                % fix for output and plotting
                if min(find(looptest==0)) == 1;
                    outnr1 = outn(9); outnr2 = outn(10);
                    if mt >= 1E5; outnr3 = outn(11); outnr4 = outn(12); end;
                    if mt >= 1E6; outnr5 = outn(13); outnr6 = outn(14); end;
                    plr = pl(3); % fix for plotting
                elseif min(find(looptest==0)) == 2;
                    outnr1 = outn(15); outnr2 = outn(16);
                    if mt >= 1E5; outnr3 = outn(17); outnr4 = outn(18); end;
                    if mt >= 1E6; outnr5 = outn(19); outnr6 = outn(20); end;
                    plr = pl(4); % fix for plotting
                elseif min(find(looptest==0)) == 3;
                    outnr1 = outn(21); outnr2 = outn(22);
                    if mt >= 1E5; outnr3 = outn(23); outnr4 = outn(24); end;
                    if mt >= 1E6; outnr5 = outn(25); outnr6 = outn(26); end;
                    plr = pl(5); % fix for plotting
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
                        ddm_min = min(ddmok')./1E2; % min sample depth (m)
                        ddm_max = max(ddmok')./1E2; % max sample depth (m)
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
                
                % plot sample depth history
                if pl_range_depth==1 && nhit>1 && (range1calc==1 || range2calc==1);
                    ddm_min = min(ddmok');
                    ddm_max = max(ddmok');
                    if ishandle(plr);
                        set(0,'currentfigure',plr);
                    else;
                        figure(plr),set(plr,'visible','off');
                        hold on;
                    end;
                    plot(tv./1E3,ddm_min./1E2,'color',clr);
                    plot(tv./1E3,ddm_max./1E2,'color',clr);
                    
                    % find min and max icev
                    icevmin = (d18O >= max(LR04ok));
                    icevmax = (d18O >= min(LR04ok));
                    icevmin = double(icevmin); % fix for matlab
                    icevmax = double(icevmax); % fix for matlab
                    
                    % fix max depth for y-axis and min and max ice cover
                    if min(find(looptest==0)) == 1;
                        pldepth3(end+1) = interp1(tv,ddm_max,min(pl_maxt,max(tv)),'pchip')./1E2;
                        plicemin3(end+1,1:max(tv)+1) = interp1(tv,icevmin,(0:1:max(tv)));
                        plicemax3(end+1,1:max(tv)+1) = interp1(tv,icevmax,(0:1:max(tv)));
                    elseif min(find(looptest==0)) == 2;
                        pldepth4(end+1) = interp1(tv,ddm_max,min(pl_maxt,max(tv)),'pchip')./1E2;
                        plicemin4(end+1,1:max(tv)+1) = interp1(tv,icevmin,(0:1:max(tv)));
                        plicemax4(end+1,1:max(tv)+1) = interp1(tv,icevmax,(0:1:max(tv)));
                    elseif min(find(looptest==0)) == 3;
                        pldepth5(end+1) = interp1(tv,ddm_max,min(pl_maxt,max(tv)),'pchip')./1E2;
                        plicemin5(end+1,1:max(tv)+1) = interp1(tv,icevmin,(0:1:max(tv)));
                        plicemax5(end+1,1:max(tv)+1) = interp1(tv,icevmax,(0:1:max(tv)));
                    end;
                end;
                
                % plot banana
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
                    
                    % calculate sigma line
                    sigmal = sigmaline(N10n,N26n,delN10n,delN26n);
                    
                    % plot banana
                    if ishandle(pl(6));
                        set(0,'currentfigure',pl(6));
                    else;
                        figure(pl(6)),set(pl(6),'visible','off');
                        hold on;
                    end
                    % plot point
                    semilogx(N10n,N26n/N10n,'.','color',clr,'markersize',15);
                    
                    if Punc == 1;
                        % plot sigmaline
                        semilogx(sigmal.x,sigmal.y,'color',clr,'linestyle','--');
                        % remove prodrate uncertainty
                        delN10n = sqrt(sample.delN10^2 - (sample.N10*delPref10/Pref10)^2) / P0avg10;
                        delN26n = sqrt(sample.delN26^2 - (sample.N26*delPref26/Pref26)^2) / P0avg26;
                        % calculate sigma line
                        sigmal = sigmaline(N10n,N26n,delN10n,delN26n);
                    end;
                    
                    % plot sigmaline
                    semilogx(sigmal.x,sigmal.y,'color',clr);
                    
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
                        [N10path N26path] = N1026calc(tv,P10r,P26r,l10,l26,P0avg10,P0avg26);
                        
                        % plot 26/10 ratios against N10
                        semilogx(N10path,N26path./N10path,'color',clr);
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
                
                partest(1:6) = 0; % reset parameter test
                looptest(min(find(looptest==0))) = 1; % reset looptest
                ddmok = []; glacEok = []; nonglacEok = []; LR04ok = []; Ndiffok = []; % reset
                nhit = 0; % reset nhit
                
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
    
    fprintf(1,'\n')
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
figh = findobj('type','figure');
if pl(1)>0 && any(figh==pl(1));
    figure(pl(1)),set(pl(1),'visible','off');
    
    % plot ice cover
    if size(plice,1) > 1; plicev = (sum(plice==1) > 0); else; plicev = (plice == 1); end;
    [plicex,plicey] = plot_ice(plicev,0,5);
    if numel(plicex) > 0;
        legh1(end+1) = patch(plicex./1E3,plicey,[0.8 0.8 1],'EdgeColor',[0.8 0.8 1]);
        legin1(end+1) = {'ice cover'};
    end;
    
    % fix and display plot
    figure(pl(1));
    box on;
    axis([0 pl_maxt/1E3 0 100]);
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('Cosmogenic concentration (%)');
    legend(legh1,legin1,'location','northwest');
    set(gca,'layer','top'); % plot axis on top
    hold off;
end;
if pl(2)>0 && any(figh==pl(2));
    figure(pl(2)),set(pl(2),'visible','off');
    
    % determine maxy for axis
    if max(pldepth2) > maxmaxy; maxy = maxmaxy; else; maxy = max(pldepth2); end;
    if maxy == 0; maxy = 1; end; % fix for zero erosion case
    
    % plot ice cover
    if size(plice,1) > 1; plicev = (sum(plice==1) > 0); else; plicev = (plice == 1); end;
    [plicex,plicey] = plot_ice(plicev,maxy,maxy*0.95);
    if numel(plicex) > 0;
        legh2(end+1) = patch(plicex./1E3,plicey,[0.8 0.8 1],'EdgeColor',[0.8 0.8 1]);
        legin2(end+1) = {'ice cover'};
    end;
    
    figure(pl(2));
    box on;
    axis([0 pl_maxt/1E3 0 maxy]);
    axis('ij');
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('Sample depth (m)');
    legend(legh2,legin2,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off;
end;
if pl(3)>0 && numel(pldepth3)>0 && any(figh==pl(3));
    figure(pl(3)),set(pl(3),'visible','off');
    
    % determine maxy for axis
    if max(pldepth3) > maxmaxy; maxy = maxmaxy; else; maxy = max(pldepth3); end;
    
    % plot ice cover
    if size(plicemin3,1)>1; pliceminv = (prod(plicemin3)==1); else; pliceminv = (plicemin3==1); end;
    if size(plicemax3,1)>1; plicemaxv = (sum(plicemax3==1)>0); else; plicemaxv=(plicemax3==1); end;
    [pliceminx,pliceminy] = plot_ice(pliceminv,maxy,maxy*0.95);
    [plicemaxx,plicemaxy] = plot_ice(plicemaxv,maxy,maxy*0.95);
    if numel(plicemaxx) > 0;
        legh3(1) = patch(plicemaxx./1E3,plicemaxy,[0.8 0.8 1],'EdgeColor',[0.8 0.8 1]);
        legin3(1) = {'max ice cover'};
    end;
    if numel(pliceminx) > 0;
        legh3(2) = patch(pliceminx./1E3,pliceminy,[0.5 0.5 1],'EdgeColor',[0.5 0.5 1]);
        legin3(2) = {'min ice cover'};
    end;
    
    figure(pl(3));
    box on;
    axis([0 pl_maxt/1E3 0 maxy]);
    axis('ij');
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('^{10}Be Sample depth (m)');
    legend(legh3,legin3,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off;
end;
if pl(4)>0 && any(figh==pl(4));
    figure(pl(4)),set(pl(4),'visible','off');
    
    % determine maxy for axis
    if max(pldepth4) > maxmaxy; maxy = maxmaxy; else; maxy = max(pldepth4); end;
    
    % plot ice cover
    if size(plicemin4,1)>1; pliceminv = (prod(plicemin4)==1); else; pliceminv = (plicemin4==1); end;
    if size(plicemax4,1)>1; plicemaxv = (sum(plicemax4==1)>0); else; plicemaxv=(plicemax4==1); end;
    [pliceminx,pliceminy] = plot_ice(pliceminv,maxy,maxy*0.95);
    [plicemaxx,plicemaxy] = plot_ice(plicemaxv,maxy,maxy*0.95);
    if numel(plicemaxx) > 0;
        legh4(1) = patch(plicemaxx./1E3,plicemaxy,[0.8 0.8 1],'EdgeColor',[0.8 0.8 1]);
        legin4(1) = {'max ice cover'};
    end;
    if numel(pliceminx) > 0;
        legh4(2) = patch(pliceminx./1E3,pliceminy,[0.5 0.5 1],'EdgeColor',[0.5 0.5 1]);
        legin4(2) = {'min ice cover'};
    end;
    
    figure(pl(4));
    box on;
    axis([0 pl_maxt/1E3 0 maxy]);
    axis('ij');
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('^{26}Al Sample depth (m)');
    legend(legh4,legin4,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off;
end;
if pl(5)>0 && any(figh==pl(5));
    figure(pl(5)),set(pl(5),'visible','off');
    
    % determine maxy for axis
    if max(pldepth5) > maxmaxy; maxy = maxmaxy; else; maxy = max(pldepth5); end;
    
    % plot ice cover
    if size(plicemin5,1)>1; pliceminv = (prod(plicemin5)==1); else; pliceminv = (plicemin5==1); end;
    if size(plicemax5,1)>1; plicemaxv = (sum(plicemax5==1)>0); else; plicemaxv=(plicemax5==1); end;
    [pliceminx,pliceminy] = plot_ice(pliceminv,maxy,maxy*0.95);
    [plicemaxx,plicemaxy] = plot_ice(plicemaxv,maxy,maxy*0.95);
    if numel(plicemaxx) > 0;
        legh5(1) = patch(plicemaxx./1E3,plicemaxy,[0.8 0.8 1],'EdgeColor',[0.8 0.8 1]);
        legin5(1) = {'max ice cover'};
    end;
    if numel(pliceminx) > 0;
        legh5(2) = patch(pliceminx./1E3,pliceminy,[0.5 0.5 1],'EdgeColor',[0.5 0.5 1]);
        legin5(2) = {'min ice cover'};
    end;
    
    figure(pl(5));
    box on;
    axis([0 pl_maxt/1E3 0 maxy]);
    axis('ij');
    set(gca(),'xdir','reverse');
    xlabel('Time (ka)');
    ylabel('^{10}Be+^{26}Al Sample depth (m)');
    legend(legh5,legin5,'location','northwest')
    set(gca,'layer','top'); % plot axis on top
    hold off;
end;
if pl(6)>0 && any(figh==pl(6));
    % create data for the simple-exposure line including ratio uncertainties
    tempt = logspace(0,7,100);
    be = (1/consts.l10)*(1-exp(-consts.l10*tempt));
    al = (1/consts.l26)*(1-exp(-consts.l26*tempt));
    al_be = al./be;
    
    % calculate Eline
    [bee,ale_bee] = Eline(plEl,consts);
    figure(pl(6));
    box on;
    % plot simple exposure line and erosion end-point line
    semilogx(be,al_be,'color','black');
    semilogx(bee,ale_bee,'linestyle','--','color','black');
    axis([1000 3000000 0.2 1.2]);
    xlabel('[^{10}Be]*');
    ylabel('[^{26}Al]*/[^{10}Be]*');
    set(gca,'layer','top'); % plot axis on top
    set(gca,'XScale','log'); % fix for matlab
    hold off;
end;

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


% subfunction Pcalc1026 ============================================================================
function out = Pcalc1026(fulldm,noicem,tv,Psp10,Psp26,Pref10,Pref26,P_mu_z10,P_mu_z26,l10,l26,...
    sample);
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


% subfunction N1026calc ============================================================================
function [N10m N26m] = N1026calc(tv,P10,P26,l10,l26,P0avg10,P0avg26);
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
            N10m(i+1,:) = (N10m(end,:) + Pmean10(i).*(tvflip(i)-tvflip(i+1))) .* ...
                exp(-(tvflip(i)-tvflip(i+1)).*l10);
            N26m(i+1,:) = (N26m(end,:) + Pmean26(i).*(tvflip(i)-tvflip(i+1))) .* ...
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
    
    % flip back and divide with average surface P
    N10m = flip(N10m);
    N26m = flip(N26m);
    N10m = N10m./P0avg10;
    N26m = N26m./P0avg26;
    N10m(end,:) = []; % skip end row to avoid division by 0
    N26m(end,:) = []; % skip end row to avoid division by 0
% end subfunction N1026calc ========================================================================


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


% subfunction nucl1_buildup ========================================================================
function out = nucl1_buildup(Eglacv,fulldm,eros,sample,Psp0,P_mu_z,noicev,tv,l,mindelE,maxdelE);
    % calculate P through time
    fulldv = interp1(Eglacv',fulldm',eros,'pchip'); % interpolate depth (g/cm2)
    dpfs = exp(-fulldv./sample.Lsp); % spallation depth dependence
    Pspv = Psp0 .* sample.othercorr .* dpfs .* noicev; % spallation prod
    Pmuv = interp1(sample.z_mu,P_mu_z,fulldv,'pchip') .* sample.othercorr .* noicev;
    Pfull = Pspv' + Pmuv'; % full prod
    % same for min and max uncertainties
    if exist('mindelE');
        mindv = interp1(Eglacv',fulldm',eros-mindelE,'pchip'); % interpolate depth (g/cm2)
        mindpfs = exp(-mindv./sample.Lsp); % spallation depth dependence
        minPspv = Psp0 .* sample.othercorr .* mindpfs .* noicev; % spallation prod
        minPmuv = interp1(sample.z_mu,P_mu_z,mindv,'pchip') .* sample.othercorr .* noicev;
        Pfull(:,2) = minPspv' + minPmuv'; % full prod
        maxdv = interp1(Eglacv',fulldm',eros+maxdelE,'pchip'); % (g/cm2)
        maxdpfs = exp(-maxdv./sample.Lsp); % spallation depth dependence
        maxPspv = Psp0 .* sample.othercorr .* maxdpfs .* noicev; % spallation prod
        maxPmuv = interp1(sample.z_mu,P_mu_z,maxdv,'pchip') .* sample.othercorr .* noicev;
        Pfull(:,3) = maxPspv' + maxPmuv'; % full prod
    end;
    
    tvflip = flip(tv');
    Pflip = flip(Pfull);
    Pmean = (Pflip(1:end-1,:)+Pflip(2:end,:))./2;
    dtv = tvflip(1:end-1) - tvflip(2:end);
    fidx = min(find(dtv(1:end-1)-dtv(2:end)==0));
    
    N(1,size(Pflip,2)) = 0;
    
    % for loop used for first part with varying time steps
    if fidx > 1;
        for i = 1:fidx-1;
            N(i+1,:) = (N(end,:) + Pmean(i).*(tvflip(i)-tvflip(i+1))) .* ...
                exp(-(tvflip(i)-tvflip(i+1)).*l);
        end;
        % remove calculated steps
        tvflip(1:fidx-1) = [];
        Pmean(1:fidx-1,:) = [];
        dtv(1:fidx-1) = [];
    end;
    
    % while loop used for filtering parts with constant time step
    while numel(tvflip) > 1;
        widx = max(find(dtv==dtv(1)));
        Ptemp = Pmean(1:widx,:).*dtv(1);
        Ptemp(2:end+1,:) = Ptemp;
        Ptemp(1,:) = N(end,:);
        decay = -exp(-l.*dtv(1));
        Ntemp = filter(1,[1,decay],Ptemp);
        N(end+1:end+widx,:) = Ntemp(2:end,:);
        % remove calculated steps
        tvflip(1:widx) = [];
        Pmean(1:widx,:) = [];
        dtv(1:widx) = [];
    end;
    
    N = flip(N);
    out.mid = N(:,1)./N(1,1).*1E2;
    if exist('mindelE');
        out.unc = N(:,2:3)./repmat(N(1,2:3),size(tvflip)).*1E2;
    end;
% end subfunction nucl1_buildup ====================================================================


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


% subfunction plot_ice =============================================================================
function [x,y] = plot_ice(plicev,y1,y2);
    pliceend = (filter([1 2],1,plicev) == 1); % find end of ice cover periods
    plicestart = flip(filter([1 2],1,flip(plicev)) == 1); % find start of ice cover periods
    pliceendidx = find(pliceend == 1)-1; % time of iceend
    plicestartidx = find(plicestart == 1)-1; % time of icestart
    x = []; y = [];
    for j = 1:numel(pliceendidx);
        x(end+1:end+4) = [pliceendidx(j) pliceendidx(j) plicestartidx(j) plicestartidx(j)+1];
        y(end+1:end+4) = [y1 y2 y2 y1];
    end;
% end subfunction plot_ice =========================================================================


% subfunction Eline ================================================================================
function [x,y] = Eline(plEl,consts);
    fprintf(1,'calculating erosion end-point line...');
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
    fprintf(1,' done!\n');
% end subfunction Eline ============================================================================
