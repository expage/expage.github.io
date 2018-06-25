% Copyright 2013, New Mexico Tech. Written by Brian Brochers, 2011.
% borchers@nmt.edu
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 3,
% as published by the Free Software Foundation (www.fsf.org).

%
% yout=interpolate(xin,yin,x)
%
% Given an ordered set of (xin(i),yin(i)) pairs, with 
%   xin(1)< xin(2) < ... <xin(n)
%
% and an ordered set of x values
%
%   x(1) < x(2) < x(3) < ... < x(m)
%
% use linear interpolation to produce yout(1), yout(2), ..., yout(m).
%
% Returns NaN for any input x(i) where x(i) < xin(1) or
% x(i)>xin(n).
%
function yout=interpolate(xin,yin,x)
%
% Check the sizes of xin and yin.
%
if (length(xin) ~= length(yin))
  error('xin and yin must have the same size!');
else
  n=length(xin);
end
%
% setup yout.
%
yout=zeros(size(x));
m=length(x);
%
% Work through x.
%
i=1;
j=1;
while (i<=m)
  if (x(i) < xin(1))
    yout(i)=NaN;
    i=i+1;
  else
    if (x(i) > xin(n))
      yout(i)=NaN;
      i=i+1;
    else
      while (xin(j+1) < x(i))
	j=j+1;
      end
%
% Now, xin(j)<=x(i)<xin(j+1)
%
      yout(i)=yin(j)+(yin(j+1)-yin(j))*(x(i)-xin(j))/(xin(j+1)-xin(j));
      i=i+1;
    end
  end
end
