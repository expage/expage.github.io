function [pstar,iter]=mucalib_lm(func,jac,p0,tol,maxiter,mucal)
% [xstar,iter]=lm(func,jac,p0,tol,maxiter,mucal)
%
% Use the Levenberg-Marquardt algorithm to minimize 
%
%  f(x)=sum(F_i(x)^2)
%
% Input Parameters:
%    func - name of the function F(x)
%    jac - name of the Jacobian function J(x)
%    p0 - initial guess
%    tol - stopping tolerance
%    maxiter - maximum number of iterations allowed
%    mucal - structure with data for muon calibration
%
% Output Parameters:
%    xstar - best solution found. 
%    iter - Iteration count.
%
% Based on lm.m from CRONUScalc: https://bitbucket.org/cronusearth/cronus-calc
% Jakob Heyman - 2016-2018 (jakob.heyman@gu.se)
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).

%
%  Initialize p and oldp.  
%
p=p0;
fp=1.0e30;             % nonsense value will go away after 1st iteration
oldp=p0*2;
oldfp=fp*2;
%
% Initialize lambda.
%
lambda=0.0001;
%
% The main loop.  While the current solution isn't good enough, keep
% trying...  Stop after maxiter iterations in the worst case.
%
iter=0;
while (1==1)
iter
%
% Compute the current f values and Jacobian.
%
  f=feval(func,p,mucal);
  fp=norm(f,2)^2;
  J=feval(jac,p,mucal);
%
% Check the termination criteria.
%
  rhs=-J'*f;
  if ((norm(rhs,2)< sqrt(tol)*(1+abs(fp))) && ...
      (abs(oldfp-fp)<tol*(1+abs(fp))) && ...
      (norm(oldp-p,2)<sqrt(tol)*(1+norm(p,2))))
     pstar=p;
     return;
  end
%
% Not yet optimal.
% 
%
% Compute rhs=-J'*f
%
%    rhs=-J'*feval(func,p);
%
% We use a clever trick here.  The least squares problem
%
%  min || [ J              ] s - [ -F ] ||
%      || [ sqrt(lambda)*I ]     [ 0  ] ||
%
% Has the normal equations solution
%
%  s=-inv(J'*J+lambda*I)*J'*F
%
% which is precisely the LM step.  We can solve this least squares problem
% more accurately using the QR factorization then by computing
% inv(J'*J+lambda*I) explicitly.
%
  myrhs=[-f; zeros(length(p),1)];
  s=[J; sqrt(lambda)*eye(length(p))]\myrhs;
%
% See whether this improves chisq or not.
%
  fnew=feval(func,p+s,mucal);
  fpnew=norm(fnew,2)^2;
%
% If this does not improve the objective value, then increase
% lambda and try again.
%
  while (fpnew > fp)
    iter=iter+1;
    if (iter < maxiter)
      lambda=lambda*2.5;
      s=[J; sqrt(lambda)*eye(length(p))]\myrhs;
      fnew=feval(func,p+s,mucal);
      fpnew=norm(fnew,2)^2;
    else
      warning('Maximum iterations exceeded');
      pstar=p;
      return;
    end
  end
%
% If we got here, then we have an improved function value.
%
  oldp=p;
  oldfp=fp;
  p=p+s;
  fp=fpnew;
  lambda=lambda/2;
  if (lambda <10^(-12))
    lambda=1.0e-12;
  end
  iter=iter+1;
%
% Check for maximum iterations exceeded.
%
   if (iter > maxiter)
     warning('Maximum iterations exceeded');
     pstar=p;
     return;
   end
end
