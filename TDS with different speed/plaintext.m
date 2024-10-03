function [sys,x0,str,ts,simStateCompliance] = impulsive(t,x,u,flag,a,b)

switch flag,

  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  case 1,
    sys=mdlDerivatives(t,x,u,a,b);

  case 2,
    sys=mdlUpdate(t,x,u,a,b);

  case 3,
    sys=mdlOutputs(t,x,u,a,b);

  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u,a,b);

  case 9,
    sys=mdlTerminate(t,x,u,a,b);

  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
;

x0  = [];

str = [];

ts  = [0 0];

simStateCompliance = 'UnknownSimState';

function sys=mdlDerivatives(t,x,u,a,b)

sys = [];

function sys=mdlUpdate(t,x,u,a,b)

sys = [];

function sys=mdlOutputs(t,x,u,a,b)

if mod(t,a)==0
    x=rand(1,1);
    if x>=0.5
        x=b;
    else
        x=0;
    end
end
    
sys = [x];

function sys=mdlGetTimeOfNextVarHit(t,x,u,a,b)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

function sys=mdlTerminate(t,x,u,a,b)

sys = [];
