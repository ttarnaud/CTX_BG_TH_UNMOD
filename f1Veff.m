function [Out] = f1Veff(Qv,US,varargin)
% US is boolean: 1 - US on, 0 - US off
% Can be called with single cell length-2 varargin (nargin=1) 'Qv' with
% Qv{1}=Q, Qv{2} = US. This is done to circumvent the split-regexp in
% Dynasim at comma's (i.e., f1Veff(Q,V) will error in Dynasim due to the
% comma, while f1Veff({Q V}) does not)
persistent QmRange VecVeff0 VecVeffPa
if isempty(QmRange)
QmRange = varargin{1};
end
if isempty(VecVeff0)
VecVeff0 = varargin{2};
end
if isempty(VecVeffPa)
VecVeffPa = varargin{3};
end

if (nargin == 1) && length(Qv) == 2
Q = Qv{1}; US = Qv{2};
else
Q = Qv;
end

if isrow(Q), Q = Q'; end

if ~~(US) 
Out = nakeinterp1(QmRange',VecVeffPa,Q)';
else
Out = nakeinterp1(QmRange',VecVeff0,Q)';
end
end
