function Out = USstep(t,varargin)
persistent USdc USpstart USpd USprf
if isempty(USdc)||isempty(USpstart)||isempty(USpd)||isempty(USprf)
USpstart = varargin{1}; USpd = varargin{2}; USdc = varargin{3}; USprf = varargin{4};
end
if USdc == 1
USstep = double(t>=USpstart&t<=USpd+USpstart); 
else
USprp = (1/USprf);      % Pulse repetition period (s)
USstep = double(mod(t-USpstart,USprp)<=USdc*USprp).*double(t>=USpstart&t<=USpd+USpstart);
end
Out = USstep;
end