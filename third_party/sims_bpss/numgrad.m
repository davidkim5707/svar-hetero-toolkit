function [g, badg] = numgrad(fcn,x,varargin)
% Source: Christopher Sims optimization suite (csminwel); lightly adapted

delta = 1e-5;  % ← Changed from 1e-6 to match R code
n=length(x);
tvec=delta*eye(n);
g=zeros(n,1);

[f0,cost_flag] = feval(fcn, x, varargin{:});

badg=0;
goog=1;
scale=1;

for i=1:n
   if size(x,1)>size(x,2)
      tvecv=tvec(i,:);
   else
      tvecv=tvec(:,i);
   end
   
   [fh,cost_flag] = feval(fcn, x+scale*transpose(tvecv), varargin{:});
   
   if cost_flag
       g0 = (fh - f0) / (scale*delta);
   else
       [fh,cost_flag] = feval(fcn, x-scale*transpose(tvecv), varargin{:});
       if cost_flag
           g0 = (f0-fh) / (scale*delta);
       else
           goog=0;
       end
   end
   
   if goog && abs(g0)< 1e15
      g(i)=g0;
   else
      disp('bad gradient ------------------------')
      g(i)=0;
      badg=1;
   end
end