function Out = funtauu(Vv,Vx)
if (nargin == 1)
V = Vv{1}; Vx = Vv{2};
else
V = Vv;
end
if (V+Vx) < -80
    Out = (1/3.7)*exp((V+Vx+467)/66.6);
else
    Out = (1/3.7)*(exp(-(V+Vx+22)/10.5)+28);
end
end