function [wmean,wunc,rchi2,Pvalue] = w_mean_unc(mid,unc,chisquare)

% Calculates weighted mean and estimates an unbiased weighted uncertainty assuming each mid and unc
% is sampled from the same distribution with mean wmean and variance wunc^2.
% The reduced chi-square value and P-value is calculated if a third input chisquare is included.

wmean = sum(mid./unc.^2) / sum(1./unc.^2);
wunc = sqrt(sum(1./unc.^2.*(mid-wmean).^2) / (sum(1./unc.^2)-(sum(1./unc.^4)/sum(1./unc.^2))));

if exist('chisquare');
    rchi2 = 1/(numel(mid)-1) .* sum(((mid-wmean)./unc).^2);
    Pvalue = 1 - chi2cdf(rchi2.*(numel(mid)-1),numel(mid)-1);
end;
