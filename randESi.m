function Out=randESi(t,varargin)
persistent rA rPS rPD rF rPW rsdF rsdPW t_end ...
    nRands rndF randESiPS randESiPW 
if isempty(rA) || isempty(rPS) || isempty(rPD) || isempty(rF) || isempty(rPW) ...
        || isempty(rsdF) || isempty(rsdPW) || isempty(t_end) 
if nargin ~= 9
    error('randESi called but not initialised');
end
rA = varargin{1}; rPS = varargin{2}; rPD = varargin{3}; 
rF = varargin{4}; rPW = varargin{5}; rsdF = varargin{6}; rsdPW = varargin{7};
t_end = varargin{8};
end

if isempty(nRands) || isempty(rndF) || isempty(randESiPS) || isempty(randESiPW)
nRands = 50; rndF = Inf;
while sum(1./rndF) < t_end
nRands=2*nRands; 
rndF = max(rF+rsdF.*randn(nRands,1),0);
end
randESiPS = [0;cumsum(1./rndF)];        % randESi pulse starts
randESiPW = max(rPW+rsdPW.*randn(nRands+1,1),0); % randESi pulse widths
end

Out = rA*(t>=rPS&(t<=(rPS+rPD)))...
    .*(mod(t,randESiPS(find(t>=randESiPS,1,'last')))<=...
    randESiPW(find(t>=randESiPS,1,'last')));
end