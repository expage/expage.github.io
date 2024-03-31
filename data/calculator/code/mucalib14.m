function mucalib14()

% function for calibrating muogenic 14C production parameters based on depth profile data
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

% Choose site. Options: LeymonHigh14
site = {'LeymonHigh14'};

% display correlation matrix for fitted parameters? (1 = yes)
display_corrmatrix = 1;

% run and load expage constants
make_consts_expage;
load consts_expage;

% declare fstarv and sigma0v vectors
fstar14v = []; fstar14uncv = [];
sigma014v = []; sigma014uncv = [];

% declare ff struct for figures
ff = struct;

% fix for output
output(1,1:12) = {'site','fstar14','fstar14unc','sigma014','sigma014unc','Lsp14(g/cm2)',...
    'Lsp14unc','E(g/cm2/ka)','Eunc','chi2','p-value','outliers(n)'};

for sn = 1:numel(site);
    % display site
    fprintf(1,'\nSite: %s\n',site{sn});

    % write site name to output
    output(sn+1,1) = site(sn);

    % load data. 14C data is derived from Lupker et al. (2015)
    mucal = mucalib_data(site{sn});

    % set uncertainty to minimum 5% (14C)
    mucal.N14unc(mucal.N14unc./mucal.N14<0.05) = mucal.N14(mucal.N14unc./mucal.N14<0.05) .* 0.05;

    % declare rm variables (for potential plotting)
    mucal.sample_name14rm = {}; mucal.Nz14rm = []; mucal.N14rm = []; mucal.N14uncrm = [];

    % calculate atmospheric pressure
    mucal.atm = ERA40atm(mucal.lat,mucal.lon,mucal.elv);

    % We have 4 parameters to fit.
    mucal.npars = 4;

    % initial guess at the parameters
    pinit = mucal.pinit;
    % Parameters in pinit are:
    %   p(1)    erosion rate ([g/cm2]/ka)
    %   p(2)    attenuation length for 14-C (g/cm2)
    %   p(3)    fstar14 (scaled by 1.0e-1)
    %   p(4)    sigma014 (scaled by 1.0e-31)

    % pre-calculated muon P at depth (without fstar and sigma0)
    mucalpre = mucalib_depthcalc_precalc(site{sn});
    mucal.dz = mucalpre.dz;
    mucal.Pfast14d = mucalpre.Pfast14d;
    mucal.Pneg14d = mucalpre.Pneg14d;

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
    fprintf(1,'14C Lsp = %s ± %s (g/cm2)\n',pst.Lsp14_s,pstunc.Lsp14_s);
    fprintf(1,'fstar14 = %s ± %s E-1\n',pst.fstar14_s,pstunc.fstar14_s);
    fprintf(1,'sigma014 = %s ± %s E-31\n',pst.sigma014_s,pstunc.sigma014_s);

    % fill fstar and sigma0 vectors
    fstar14v(end+1) = pst.fstar14 / 1e-1;
    fstar14uncv(end+1) = pstunc.fstar14 / 1e-1;
    sigma014v(end+1) = pst.sigma014 / 1e-31;
    sigma014uncv(end+1) = pstunc.sigma014 / 1e-31;

    % fill output
    output(sn+1,2:3) = {[pst.fstar14_s 'e-1'],[pstunc.fstar14_s 'e-1']};
    output(sn+1,4:5) = {[pst.sigma014_s 'e-31'],[pstunc.sigma014_s 'e-31']};
    output(sn+1,6:7) = {pst.Lsp14_s,pstunc.Lsp14_s};
    output(sn+1,8:9) = {pst.erosion_s,pstunc.erosion_s};
    output(sn+1,10:11) = {num2str(chi2,'%.3f'),num2str(pvalue,'%.3f')};
    output(sn+1,12) ={num2str(numel(mucal.N14rm),'%.0f')};

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
    [wfstar14,wfstar14unc] = w_mean_unc(fstar14v,fstar14uncv);
    [wsigma014,wsigma014unc] = w_mean_unc(sigma014v,sigma014uncv);
    % fill output
    output(sn+2,2:3) = {[num2str(wfstar14,'%.3f') 'e-1'],[num2str(wfstar14unc,'%.3f') 'e-1']};
    output(sn+2,4:5) = {[num2str(wsigma014,'%.3f') 'e-31'],[num2str(wsigma014unc,'%.3f') 'e-31']};
    output(sn+2,[1 (6:1:12)]) = {'Wmean±Wunc','-','-','-','-','-','-'};
    % display error-weighted fstar and sigma0 values
    fprintf(1,'\nWeighted mean parameters:\n');
    fprintf(1,'fstar14 = %.3f ± %.3f E-1\n',wfstar14,wfstar14unc);
    fprintf(1,'sigma014 = %.3f ± %.3f E-31\n',wsigma014,wsigma014unc);
end;

% plot depth profiles
for sn = 1:numel(site);
    % fix depth vector and calculate best fit depth profile conc
    plotz = (0:10:5000);
    conc_fit14 = get_z_conc(pl.mucal.(site{sn}),plotz,pl.pst.(site{sn}),'14');

    % use error-weighted fstar and sigma0 values and calculate depth profile conc
    % (here we use the best fit erosion rate and attenuation length)
    if numel(site) > 1;
        pl.pst.(site{sn}).fstar14 = wfstar14 .* 1e-1;
        pl.pst.(site{sn}).sigma014 = wsigma014 .* 1e-31;
    else;
        pl.pst.(site{sn}).fstar14 = fstar14v .* 1e-1;
        pl.pst.(site{sn}).sigma014 = sigma014v .* 1e-31;
    end
    conc_w14 = get_z_conc(pl.mucal.(site{sn}),plotz,pl.pst.(site{sn}),'14');

    % plot depth profiles
    plot_depth_profiles(pl.mucal.(site{sn}),conc_fit14,conc_w14,plotz,site{sn},[1E3 5E5 0 5000],...
        '14','C');
end;

% fix and save output
outstr = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
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
    nsamples14 = length(mucal.N14);
    mucal.nresids = nsamples14;

    % Age Relative to t0=2010 - LSD tv from LSD_fix
    % tv = [0:10:50 60:100:2960 3060:200:74860 75060:1000:799060 800060:2000:2000060 1E7];

    % Fix w,Rc,SPhi, for sp and nu prod rate scaling 10 Ma back in time
    LSDfix = LSD_fix(mucal.lat,mucal.lon,1E7,-1,mucal.samplingyr,consts);

    mucal.tv = LSDfix.tv;

    % decay factor
    mucal.dcf14 = exp(-mucal.tv.*consts.l14);

    % surface spallation production
    Psp = P_sp_expage(mucal.atm,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,0,0,1);
    mucal.Psp14 = Psp.sp14 .* consts.Pref14;

    % muon production (without sigma0 and fstar!)
    % This is pre-computed and the data is included in mucalib_depthcalc_precalc.m to save time
    %
    %Pmud = mucalib_Pmu(mucal.dz,mucal.atm,LSDfix.RcEst,consts.SPhiInf,0,0,1,consts);
    %mucal.Pfast14d = Pmud.P_fast14';
    %mucal.Pneg14d = Pmud.P_neg14';

    % Use LM to find optimal values of erosion rate and attenuation length.
    [pstar,iter]=mucalib_lm('mucalib14_fun','mucalib14_jac',mucal.pinit,1.0e-5,50,mucal);

    % Compute the residual and J at the optimal parameters.
    rstar=mucalib14_fun(pstar,mucal);
    Jstar=mucalib14_jac(pstar,mucal);

    % Compute Chi2 and pvalue.
    chi2=norm(rstar,2)^2;
    pvalue=1-chi2cdf(chi2,mucal.nresids-mucal.npars);
% end subfunction fit_params =======================================================================


% subfunction get_pstar_names ======================================================================
function pst = get_pstar_names(pstar);
    % multiply values with exponents
    pst.erosion = pstar(1) .* 1e-3;
    pst.Lsp14 = pstar(2);
    pst.fstar14 = pstar(3) .* 1e-1;
    pst.sigma014 = pstar(4) .* 1e-31;
    % fix as strings
    pst.erosion_s = num2str(pstar(1),'%.3f');
    pst.Lsp14_s = num2str(pstar(2),'%.2f');
    pst.fstar14_s = num2str(pstar(3),'%.3f');
    pst.sigma014_s = num2str(pstar(4),'%.3f');
% end subfunction get_pstar_names ==================================================================


% subfunction get_cluster ==========================================================================
function [mucal,pstar,rstar,Jstar,chi2,pvalue] = get_cluster(mucal,consts,pstar,pvalue,cl);
    % save input mucal
    mucal_in = mucal;
    
    % define maximum number of outliers
    remove14 = floor(numel(mucal.N14)*cl.maxoutratio);
    r14 = 1; % outlier counter

    % outlier counter
    outn = 0;

    % loop for outlier removal
    while pvalue<cl.minpvalue && remove14>=r14;
        % add 1 to counter
        outn = outn + 1;
        
        % get pstar parameters with names multiplied with exponents and as strings
        pst = get_pstar_names(pstar);

        % calculate expected concentration for all sample depths
        conc14 = get_z_conc(mucal,mucal.Nz14',pst,'14');
        
        % calculate deviation for all individual samples
        chidev14 = ((mucal.N14-conc14.Nfull')./mucal.N14unc).^2;

        % find index of sample with largest dev
        [value14,rm_idx14] = max(chidev14);

        % remove outlier and save in rm variables
        fprintf(1,'outlier %.0f removed: %s - 14C\n',outn,mucal.sample_name14{rm_idx14});
        mucal = remove_outlier(mucal,rm_idx14,'14');

        % fit parameters
        [mucal,pstar,rstar,Jstar,chi2,pvalue] = fit_params(mucal,consts);
    end;

    % if no cluster: redo parameter fitting and remove all samples as outliers
    if pvalue < cl.minpvalue;
        fprintf(1,'\nno clustering for site!\n');
        mucal = mucal_in;
        [mucal,pstar,rstar,Jstar,chi2,pvalue] = fit_params(mucal,consts);
        mucal = remove_outlier(mucal,(1:1:numel(mucal.N14)),'14');
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
    % plot 14C figure
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
        Nuncrmv(Nuncrmv<1) = 1;
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
    Nuncv(Nuncv<1) = 1;
    Nzv = [mucal.(['Nz' nn])';mucal.(['Nz' nn])'];
    semilogx(Nuncv,Nzv,'color','black');
    % plot points
    leg(end+1) = semilogx(mucal.(['N' nn]),mucal.(['Nz' nn]),'o','color','black','markersize',2);
    set(gca,'xaxislocation','top');
    set(gca,'XScale','log'); % fix for matlab
    xlabel(['^{' nn '}' nucl ' (atoms/g)']);
    ylabel('Shielding depth (g/cm^{2})');
    legend(leg,legin,'location','southeast');
    axis(axisv,'ij','square');
    hold off;
% end subfunction plot_depth_profiles ==============================================================
