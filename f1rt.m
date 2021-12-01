function [Out] = f1rt(Q,SONICrate,US,varargin)
% US is boolean: 1 - US on, 0 - US off
persistent QmRange Vecrt0 VecrtPa 
if isempty(QmRange)
QmRange = varargin{1};
end
if isempty(Vecrt0)
Vecrt0 = varargin{2};
end
if isempty(VecrtPa)
VecrtPa = varargin{3};
end

if isrow(Q), Q = Q'; end

if ~isempty(SONICrate)
if ~~(US) 
Out = nakeinterp1(QmRange',VecrtPa.(SONICrate),Q)';
else
Out = nakeinterp1(QmRange',Vecrt0.(SONICrate),Q)';
end
end
end
