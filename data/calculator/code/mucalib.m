function mucalib()

% function for calibrating muogenic 10Be and 26Al production parameters based on depth profile data
% based on calbhcore.m from CRONUScalc: https://bitbucket.org/cronusearth/cronus-calc
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
% Jakob Heyman - 2016-2024 (jakob.heyman@gu.se)

clear;
close all;

% What version is this?
ver = '202403';

% parameters for cluster test / outlier rejection
cl.minpvalue = 0.05; % lower limit for chi-square p-value
cl.maxoutratio = 1/3; % maximum outlier ratio

% Choose site. Options: Beacon, LeymonHigh, LeymonLow, LaCiotat
site = {'Beacon','LeymonHigh','LeymonLow','LaCiotat'};

% display correlation matrix for fitted parameters? (1 = yes)
display_corrmatrix = 1;

% run and load expage constants
make_consts_expage;
load consts_expage;

% declare fstarv and sigma0v vectors
fstar10v = []; fstar10uncv = [];
sigma010v = []; sigma010uncv = [];
fstar26v = []; fstar26uncv = [];
sigma026v = []; sigma026uncv = [];

% declare ff struct for figures
ff = struct;

% fix for output
output(1,1:18) = {'site','fstar10','fstar10unc','sigma010','sigma010unc','Lsp10(g/cm2)',...
    'Lsp10unc','fstar26','fstar26unc','sigma026','sigma026unc','Lsp26(g/cm2)','Lsp26unc',...
    'E(g/cm2/ka)','Eunc','chi2','p-value','outliers(n)'};

for sn = 1:numel(site);
    % display site
    fprintf(1,'\nSite: %s\n',site{sn});

    % write site name to output
    output(sn+1,1) = site(sn);

    % load data. All data is derived from Balco (2017):
    % http://hess.ess.washington.edu/repository/muons2016/
    mucal = mucalib_data(site{sn});

    % convert 10Be concentrations according to standards
    [testi,stdi] = ismember(mucal.std10,consts.std10); % find index of standard conversion factors
    mult10 = consts.std10_cf(stdi); % pick out conversion factor
    mucal.N10 = mucal.N10 .* mult10;
    mucal.N10unc = mucal.N10unc .* mult10;

    % convert 26Al concentrations according to standards
    [testi,stdi] = ismember(mucal.std26,consts.std26); % find index of standard conversion factors
    mult26 = consts.std26_cf(stdi); % pick out conversion factor
    mucal.N26 = mucal.N26 .* mult26;
    mucal.N26unc = mucal.N26unc .* mult26;

    % set uncertainty to minimum 2.9% (10Be) or 4.9% (26Al)
    mucal.N10unc(mucal.N10unc./mucal.N10<0.029) = mucal.N10(mucal.N10unc./mucal.N10<0.029) .* 0.029;
    mucal.N26unc(mucal.N26unc./mucal.N26<0.049) = mucal.N26(mucal.N26unc./mucal.N26<0.049) .* 0.049;

    % declare rm variables (for potential plotting)
    mucal.sample_name10rm = {}; mucal.Nz10rm = []; mucal.N10rm = []; mucal.N10uncrm = [];
    mucal.sample_name26rm = {}; mucal.Nz26rm = []; mucal.N26rm = []; mucal.N26uncrm = [];

    % calculate atmospheric pressure
    mucal.atm = ERA40atm(mucal.lat,mucal.lon,mucal.elv);

    % We have 7 parameters to fit.
    mucal.npars = 7;

    % initial guess at the parameters
    pinit = mucal.pinit;
    % Parameters in pinit are:
    %   p(1)    erosion rate ([g/cm2]/ka)
    %   p(2)    attenuation length for 10-Be (g/cm2)
    %   p(3)    attenuation length for 26-Al (g/cm2)
    %   p(4)    fstar10 (scaled by 1.0e-3)
    %   p(5)    sigma010 (scaled by 1.0e-31)
    %   p(6)    fstar26 (scaled by 1.0e-2)
    %   p(7)    sigma026 (scaled by 1.0e-30)

    % pre-calculated muon P at depth (without fstar and sigma0)
    mucalpre = mucalib_depthcalc_precalc(site{sn});
    mucal.dz = mucalpre.dz;
    mucal.Pfast10d = mucalpre.Pfast10d;
    mucal.Pneg10d = mucalpre.Pneg10d;
    mucal.Pfast26d = mucalpre.Pfast26d;
    mucal.Pneg26d = mucalpre.Pneg26d;

    % fit parameters
    [mucal,pstar,rstar,Jstar,chi2,pvalue] = fit_params(mucal,consts);

    % check pvalue and remove outliers
    if pvalue < cl.minpvalue;
        [mucal,pstar,rstar,Jstar,chi2,pvalue] = get_cluster(mucal,consts,pstar,pvalue,cl);
    end;

    % Compute the covariance matrix for the fitted parameters.
    covp=inv(Jstar'*Jstar);
    sigmapstar=sqrt(diag(covp));

    % Compute the correlation matrix for the fitted parameters.
    for i=1:mucal.npars 
        for j=1:mucal.npars
            corp(i,j)=covp(i,j)/(sqrt(covp(i,i))*sqrt(covp(j,j)));
        end
    end
    if display_corrmatrix == 1;
        disp('Correlations between fitted parameters');
        corp
    end;

    % get pstar parameters with names multiplied with exponents and as strings
    pst = get_pstar_names(pstar);
    pstunc = get_pstar_names(sigmapstar);

    % display the fitted parameters.
    fprintf(1,'Chi2 = %f   p-value = %f\n',[chi2; pvalue]);
    fprintf(1,'Erosion Rate = %s ± %s ((g/cm2)/ka)\n',pst.erosion_s,pstunc.erosion_s);
    fprintf(1,'10Be Lsp = %s ± %s (g/cm2)\n',pst.Lsp10_s,pstunc.Lsp10_s);
    fprintf(1,'26Al Lsp = %s ± %s (g/cm2)\n',pst.Lsp26_s,pstunc.Lsp26_s);
    fprintf(1,'fstar10 = %s ± %s E-3\n',pst.fstar10_s,pstunc.fstar10_s);
    fprintf(1,'sigma010 = %s ± %s E-31\n',pst.sigma010_s,pstunc.sigma010_s);
    fprintf(1,'fstar26 = %s ± %s E-2\n',pst.fstar26_s,pstunc.fstar26_s);
    fprintf(1,'sigma026 = %s ± %s E-30\n',pst.sigma026_s,pstunc.sigma026_s);

    % fill fstar and sigma0 vectors
    fstar10v(end+1) = pst.fstar10 / 1e-3;
    fstar10uncv(end+1) = pstunc.fstar10 / 1e-3;
    sigma010v(end+1) = pst.sigma010 / 1e-31;
    sigma010uncv(end+1) = pstunc.sigma010 / 1e-31;
    fstar26v(end+1) = pst.fstar26 / 1e-2;
    fstar26uncv(end+1) = pstunc.fstar26 / 1e-2;
    sigma026v(end+1) = pst.sigma026 / 1e-30;
    sigma026uncv(end+1) = pstunc.sigma026 / 1e-30;

    % fill output
    output(sn+1,2:3) = {[pst.fstar10_s 'e-3'],[pstunc.fstar10_s 'e-3']};
    output(sn+1,4:5) = {[pst.sigma010_s 'e-31'],[pstunc.sigma010_s 'e-31']};
    output(sn+1,6:7) = {pst.Lsp10_s,pstunc.Lsp10_s};
    output(sn+1,8:9) = {[pst.fstar26_s 'e-2'],[pstunc.fstar26_s 'e-2']};
    output(sn+1,10:11) = {[pst.sigma026_s 'e-30'],[pstunc.sigma026_s 'e-30']};
    output(sn+1,12:13) = {pst.Lsp26_s,pstunc.Lsp26_s};
    output(sn+1,14:15) = {pst.erosion_s,pstunc.erosion_s};
    output(sn+1,16:17) = {num2str(chi2,'%.3f'),num2str(pvalue,'%.3f')};
    output(sn+1,18) ={num2str(numel(mucal.N10rm)+numel(mucal.N26rm),'%.0f')};

    % save mucal and pst for plotting
    pl.mucal.(site{sn}) = mucal;
    pl.pst.(site{sn}) = pst;

    % clear mucal and pst
    clear mucal;
    clear pst;
end;

% calculate error-weighted fstar and sigma0
if numel(site) > 1;
    % calculate weighted mean and unc
    [wfstar10,wfstar10unc] = w_mean_unc(fstar10v,fstar10uncv);
    [wsigma010,wsigma010unc] = w_mean_unc(sigma010v,sigma010uncv);
    [wfstar26,wfstar26unc] = w_mean_unc(fstar26v,fstar26uncv);
    [wsigma026,wsigma026unc] = w_mean_unc(sigma026v,sigma026uncv);
    % fill output
    output(sn+2,2:3) = {[num2str(wfstar10,'%.3f') 'e-3'],[num2str(wfstar10unc,'%.3f') 'e-3']};
    output(sn+2,4:5) = {[num2str(wsigma010,'%.3f') 'e-31'],[num2str(wsigma010unc,'%.3f') 'e-31']};
    output(sn+2,8:9) = {[num2str(wfstar26,'%.3f') 'e-2'],[num2str(wfstar26unc,'%.3f') 'e-2']};
    output(sn+2,10:11) = {[num2str(wsigma026,'%.3f') 'e-30'],[num2str(wsigma026unc,'%.3f') 'e-30']};
    output(sn+2,[1 6 7 (12:1:18)]) = {'Wmean±Wunc','-','-','-','-','-','-','-','-','-'};
    % display error-weighted fstar and sigma0 values
    fprintf(1,'\nWeighted mean parameters:\n');
    fprintf(1,'fstar10 = %.3f ± %.3f E-3\n',wfstar10,wfstar10unc);
    fprintf(1,'sigma010 = %.3f ± %.3f E-31\n',wsigma010,wsigma010unc);
    fprintf(1,'fstar26 = %.3f ± %.3f E-2\n',wfstar26,wfstar26unc);
    fprintf(1,'sigma026 = %.3f ± %.3f E-30\n',wsigma026,wsigma026unc);
end;

% plot depth profiles
for sn = 1:numel(site);
    % fix depth vector and calculate best fit depth profile conc
    plotz = (0:10:7000);
    conc_fit10 = get_z_conc(pl.mucal.(site{sn}),plotz,pl.pst.(site{sn}),'10');
    conc_fit26 = get_z_conc(pl.mucal.(site{sn}),plotz,pl.pst.(site{sn}),'26');

    % use error-weighted fstar and sigma0 values and calculate depth profile conc
    % (here we use the best fit erosion rate and attenuation length)
    if numel(site) > 1;
        pl.pst.(site{sn}).fstar10 = wfstar10 .* 1e-3;
        pl.pst.(site{sn}).sigma010 = wsigma010 .* 1e-31;
        pl.pst.(site{sn}).fstar26 = wfstar26 .* 1e-2;
        pl.pst.(site{sn}).sigma026 = wsigma026 .* 1e-30;
    else;
        pl.pst.(site{sn}).fstar10 = fstar10v .* 1e-3;
        pl.pst.(site{sn}).sigma010 = sigma010v .* 1e-31;
        pl.pst.(site{sn}).fstar26 = fstar26v .* 1e-2;
        pl.pst.(site{sn}).sigma026 = sigma026v .* 1e-30;
    end;
    conc_w10 = get_z_conc(pl.mucal.(site{sn}),plotz,pl.pst.(site{sn}),'10');
    conc_w26 = get_z_conc(pl.mucal.(site{sn}),plotz,pl.pst.(site{sn}),'26');

    % plot depth profiles
    plot_depth_profiles(pl.mucal.(site{sn}),conc_fit10,conc_w10,plotz,site{sn},[1E3 1E8 0 7000],...
        '10','Be');
    plot_depth_profiles(pl.mucal.(site{sn}),conc_fit26,conc_w26,plotz,site{sn},[1E4 5E8 0 7000],...
        '26','Al');
end;

% fix and save output
outstr = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
% write out-expage.txt
out = fopen('out-mucalib.txt','w');
for i = 1:size(output,1);
    fprintf(out,outstr,output{i,:});
end;
fclose(out);
% end mucalib function =============================================================================


% subfunction fit_params ===========================================================================
function [mucal,pstar,rstar,Jstar,chi2,pvalue] = fit_params(mucal,consts);
    % Get the number of samples.
    nsamples10 = length(mucal.N10);
    nsamples26 = length(mucal.N26);
    mucal.nresids = nsamples10+nsamples26;

    % Age Relative to t0=2010 - LSD tv from LSD_fix
    % tv = [0:10:50 60:100:2960 3060:200:74860 75060:1000:799060 800060:2000:2000060 1E7];

    % Fix w,Rc,SPhi, for sp and nu prod rate scaling 10 Ma back in time
    LSDfix = LSD_fix(mucal.lat,mucal.lon,1E7,-1,mucal.samplingyr,consts);

    mucal.tv = LSDfix.tv;

    % decay factors
    mucal.dcf10 = exp(-mucal.tv.*consts.l10);
    mucal.dcf26 = exp(-mucal.tv.*consts.l26);

    % surface spallation production
    Psp = P_sp_expage(mucal.atm,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,1,1,0);
    mucal.Psp10 = Psp.sp10 .* consts.Pref10;
    mucal.Psp26 = Psp.sp26 .* consts.Pref26;

    % muon production (without sigma0 and fstar!)
    % This is pre-computed and the data is included in mucalib_depthcalc_precalc.m to save time
    %
    %Pmud = mucalib_Pmu(mucal.dz,mucal.atm,LSDfix.RcEst,consts.SPhiInf,1,1,consts);
    %mucal.Pfast10d = Pmud.P_fast10';
    %mucal.Pneg10d = Pmud.P_neg10';
    %mucal.Pfast26d = Pmud.P_fast26';
    %mucal.Pneg26d = Pmud.P_neg26';

    % Use LM to find optimal values of erosion rate and attenuation length.
    [pstar,iter]=mucalib_lm('mucalib_fun','mucalib_jac',mucal.pinit,1.0e-5,100,mucal);

    % Compute the residual and J at the optimal parameters.
    rstar=mucalib_fun(pstar,mucal);
    Jstar=mucalib_jac(pstar,mucal);

    % Compute Chi2 and pvalue.
    chi2=norm(rstar,2)^2;
    pvalue=1-chi2cdf(chi2,mucal.nresids-mucal.npars);
% end subfunction fit_params =======================================================================


% subfunction get_pstar_names ======================================================================
function pst = get_pstar_names(pstar);
    % multiply values with exponents
    pst.erosion = pstar(1) .* 1e-3;
    pst.Lsp10 = pstar(2);
    pst.Lsp26 = pstar(3);
    pst.fstar10 = pstar(4) .* 1e-3;
    pst.sigma010 = pstar(5) .* 1e-31;
    pst.fstar26 = pstar(6) .* 1e-2;
    pst.sigma026 = pstar(7) .* 1e-30;
    % fix as strings
    pst.erosion_s = num2str(pstar(1),'%.3f');
    pst.Lsp10_s = num2str(pstar(2),'%.2f');
    pst.Lsp26_s = num2str(pstar(3),'%.2f');
    pst.fstar10_s = num2str(pstar(4),'%.3f');
    pst.sigma010_s = num2str(pstar(5),'%.3f');
    pst.fstar26_s = num2str(pstar(6),'%.3f');
    pst.sigma026_s = num2str(pstar(7),'%.3f');
% end subfunction get_pstar_names ==================================================================


% subfunction get_cluster ==========================================================================
function [mucal,pstar,rstar,Jstar,chi2,pvalue] = get_cluster(mucal,consts,pstar,pvalue,cl);
    % save input mucal
    mucal_in = mucal;
    
    % define maximum number of outliers
    remove10 = floor(numel(mucal.N10)*cl.maxoutratio);
    remove26 = floor(numel(mucal.N26)*cl.maxoutratio);
    r10 = 1; % outlier counter
    r26 = 1; % outlier counter

    % outlier counter
    outn = 0;

    % loop for outlier removal
    while pvalue<cl.minpvalue && remove10>=r10 && remove26>=r26;
        % add 1 to counter
        outn = outn + 1;
        
        % get pstar parameters with names multiplied with exponents and as strings
        pst = get_pstar_names(pstar);

        % calculate expected concentration for all sample depths
        conc10 = get_z_conc(mucal,mucal.Nz10',pst,'10');
        conc26 = get_z_conc(mucal,mucal.Nz26',pst,'26');
        
        % calculate deviation for all individual samples
        chidev10 = ((mucal.N10-conc10.Nfull')./mucal.N10unc).^2;
        chidev26 = ((mucal.N26-conc26.Nfull')./mucal.N26unc).^2;

        % find index of sample with largest dev
        [value10,rm_idx10] = max(chidev10);
        [value26,rm_idx26] = max(chidev26);

        % remove outlier and save in rm variables
        if value10 > value26;
            fprintf(1,'outlier %.0f removed: %s - 10Be\n',outn,mucal.sample_name10{rm_idx10});
            mucal = remove_outlier(mucal,rm_idx10,'10');
        else;
            fprintf(1,'outlier %.0f removed: %s - 26Al\n',outn,mucal.sample_name26{rm_idx26});
            mucal = remove_outlier(mucal,rm_idx26,'26');
        end;

        % fit parameters
        [mucal,pstar,rstar,Jstar,chi2,pvalue] = fit_params(mucal,consts);
    end;

    % if no cluster: redo parameter fitting and remove all samples as outliers
    if pvalue < cl.minpvalue;
        fprintf(1,'\nno clustering for site!\n');
        mucal = mucal_in;
        [mucal,pstar,rstar,Jstar,chi2,pvalue] = fit_params(mucal,consts);
        mucal = remove_outlier(mucal,(1:1:numel(mucal.N10)),'10');
        mucal = remove_outlier(mucal,(1:1:numel(mucal.N26)),'26');
    end;
% end subfunction get_cluster ======================================================================


% subfunction get_z_conc ===========================================================================
function conc = get_z_conc(mucal,zv,pst,nucl);
    % fix depth vector and depth matrix for tv
    tvzv = pst.erosion .* mucal.tv; % depth vector for tv (g/cm2)
    tvzm = bsxfun(@plus,tvzv',zv); % depth matrix for tv (one col per depth sample)

    % sp production
    Psprm = exp(-tvzm./pst.(['Lsp' nucl])); % spallation depth prod ratio
    Pspm = bsxfun(@times,Psprm,mucal.(['Psp' nucl])'); % spallation production matrix

    % muon prod (Pfast and Pneg surface prod)
    Pfastm = interp1(mucal.dz,mucal.(['Pfast' nucl 'd']),tvzm,'pchip') .* pst.(['sigma0' nucl]);
    Pnegm = interp1(mucal.dz,mucal.(['Pneg' nucl 'd']),tvzm,'pchip') .* pst.(['fstar' nucl]);

    % full production plus decay
    Pfullm = Pspm + Pfastm + Pnegm; % full P matrix
    Pfull_lm = bsxfun(@times,Pfullm,mucal.(['dcf' nucl])'); % full P including decay matrix

    % calculate full N by integration
    conc.Nfull = trapz(mucal.tv',Pfull_lm); % calculated N

    % calculate sp production and N
    Psp_lm = bsxfun(@times,Pspm,mucal.(['dcf' nucl])'); % Psp including decay matrix
    conc.Nsp = trapz(mucal.tv',Psp_lm); % calculated N

    % calculate mu production and N
    Pfast_lm = bsxfun(@times,Pfastm,mucal.(['dcf' nucl])'); % Pfast including decay matrix
    Pneg_lm = bsxfun(@times,Pnegm,mucal.(['dcf' nucl])'); % Pneg including decay matrix
    conc.Nfast = trapz(mucal.tv',Pfast_lm); % calculated N
    conc.Nneg = trapz(mucal.tv',Pneg_lm); % calculated N
% end subfunction get_z_conc =======================================================================


% subfunction remove_outlier =======================================================================
function mucal = remove_outlier(mucal,rm_idx,nucl);
    % save sample to remove for plotting
    mucal.(['sample_name' nucl 'rm'])(end+1) = mucal.(['sample_name' nucl])(rm_idx);
    mucal.(['Nz' nucl 'rm'])(end+1) = mucal.(['Nz' nucl])(rm_idx);
    mucal.(['N' nucl 'rm'])(end+1) = mucal.(['N' nucl])(rm_idx);
    mucal.(['N' nucl 'uncrm'])(end+1) = mucal.(['N' nucl 'unc'])(rm_idx);

    % remove sample
    mucal.(['sample_name' nucl])(rm_idx) = [];
    mucal.(['Nz' nucl])(rm_idx) = [];
    mucal.(['N' nucl])(rm_idx) = [];
    mucal.(['N' nucl 'unc'])(rm_idx) = [];
% end subfunction remove_outlier ===================================================================


% subfunction plot_depth_profiles ==================================================================
function plot_depth_profiles(mucal,conc_fit,conc_w,plotz,site,axisv,nn,nucl);
    % plot depth profile
    figure('name',[site ' ' nn nucl],'NumberTitle','off'); hold on; box on;
    legin = {'spallation','slow muons','fast muons','site fit','weighted fit',['sample N' nn]};
    leg = semilogx(conc_fit.Nsp,plotz,'color','red');
    leg(end+1) = semilogx(conc_fit.Nneg,plotz,'color','green');
    leg(end+1) = semilogx(conc_fit.Nfast,plotz,'color','blue');
    leg(end+1) = semilogx(conc_fit.Nfull,plotz,'color','black');
    leg(end+1) = semilogx(conc_w.Nfull,plotz,'color',[1 0.8 0]);
    if numel(mucal.(['N' nn 'rm'])) > 0;
        % fix and plot uncertainties
        Nuncrmv = [mucal.(['N' nn 'rm'])-mucal.(['N' nn 'uncrm']);mucal.(['N' nn 'rm'])+...
            mucal.(['N' nn 'uncrm'])];
        Nzrmv = [mucal.(['Nz' nn 'rm']);mucal.(['Nz' nn 'rm'])];
        semilogx(Nuncrmv,Nzrmv,'color',[0.7 0.7 0.7]);
        % plot points
        leg(end+1) = semilogx(mucal.(['N' nn 'rm']),mucal.(['Nz' nn 'rm']),'.','color',...
            [0.7 0.7 0.7],'markersize',15);
        legin = {'spallation','slow muons','fast muons','site fit','weighted fit','outlier',...
            ['sample N' nn]};
    end;
    % fix and plot uncertainties
    Nuncv = [mucal.(['N' nn])'-mucal.(['N' nn 'unc'])';mucal.(['N' nn])'+mucal.(['N' nn 'unc'])'];
    Nzv = [mucal.(['Nz' nn])';mucal.(['Nz' nn])'];
    semilogx(Nuncv,Nzv,'color','black');
    % plot points
    leg(end+1) = semilogx(mucal.(['N' nn]),mucal.(['Nz' nn]),'.','color','black','markersize',15);
    set(gca,'xaxislocation','top');
    set(gca,'XScale','log'); % fix for matlab
    xlabel(['^{' nn '}' nucl ' (atoms/g)']);
    ylabel('Shielding depth (g/cm^{2})');
    legend(leg,legin,'location','southeast');
    axis(axisv,'ij');
    hold off;
% end subfunction plot_depth_profiles ==============================================================
