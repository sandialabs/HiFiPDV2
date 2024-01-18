% This function calculates the 10-50-90 rise times from a velocity profile
%   Author:  Travis J. Voorhees [email: tjvoorh@sandia.gov]
% =========================================================================
% How does this script work?:
% - This script finds the point in time where the sharpest (and greatest)
%   increase in velocity occurs, then crops your data to that section to
%   find the 10-50-90 points.
% - 
% =========================================================================
% Assumptions:
% - The sharpest and greatest velocity increase occurs in shock front.
% - 
% =========================================================================
% When will this code not work?:
% - When any of the assumptions above are not true
% - 
% =========================================================================


function [Rise10,Rise50,Rise90,PeakVelocity] = CalculateRiseTimes(Time,Velocity)
% check that time and velocity were input.  If not, push error
narginchk(2,2)

% Determine time-scale of input data
%  Calculate the order of magnitude of point-to-point time resolution
tscale = floor(log10(abs(mean(diff(Time)))));
switch tscale
    case {0,1} % nanoseconds
        %disp('Time-scale determined as nanoseconds.')
        ConvFactor = 1;
    case {-3,-2} % microseconds
        %disp('Time-scale determined as microseconds.')
        ConvFactor = 1e-3;
    case {-6,-5} % miliseconds
        %disp('Time-scale determined as milliseconds.')
        ConvFactor = 1e-6;
    case {-9,-8} % seconds
        %disp('Time-scale determined as seconds.')
        ConvFactor = 1e-9;
    otherwise % not detected
        CF = input(['Please enter the time-scale of your input velocity' ...
            ' profile (ns,us,ms,s,etc.):  ']);
        switch CF
            case {'seconds','s','sec',1}
                %disp('Time-scale determined as seconds.')
                ConvFactor = 1e-9;
            case {'milliseconds','miliseconds','ms',1e-3}
                %disp('Time-scale determined as milliseconds.')
                ConvFactor = 1e-6;
            case {'micro','microseconds','us','mus','\mus',1e-6}
                %disp('Time-scale determined as microseconds.')
                ConvFactor = 1e-3;
            case {'nanoseconds','ns','nano',1e-9}
                %disp('Time-scale determined as nanoseconds.')
                ConvFactor = 1;
            case {'picoseconds','ps',1e-12}
                %disp('Time-scale determined as picoseconds.')
                ConvFactor = 1e3;
            otherwise
                error('User input unusable time-scale');
        end    
end

% Find where the sharpest increase in velocity occurs
%  calculate derivative
tshift = (Time(1:(end-1))+Time(2:end))./2;
dV = diff(Velocity);
dT = diff(Time);
Accel = dV./dT;
[~,MaxInd] = max(Accel);


% Find indices of data 100 ns to left and right of biggest jump
[~,LeftCrop] = min(abs(tshift-tshift(MaxInd)+(100*ConvFactor)));
[~,RightCrop] = min(abs(tshift-tshift(MaxInd)-(100*ConvFactor)));


% Find the peak velocity following that increase in velocity
%   Assume:  Peak vel. occurs when accel decreases to consistent low value
tshift2 = (tshift(1:(end-1))+tshift(2:end))./2;
dA = diff(Accel);
dT2 = diff(tshift);
Jerk = dA./dT2;
[~,MaxDecrease] = max(-Jerk(LeftCrop:(RightCrop-1)));
MaxDecrease = MaxDecrease+LeftCrop;

% Find the peak velocity and its index
if MaxInd<26
    L = 1;
else
    L = MaxInd-25;
end
if MaxInd+25>length(Velocity)
    R = RightCrop-1;
else
    R = MaxInd+25;
end
[PeakVelocity,PeakInd]= max(movmean(Velocity(L:R),5));
PeakInd = PeakInd+MaxInd-25;

% Calculate the 10-50-90 rise times
%   First find when (velocity-%rise) crosses zero
Cross10 = diff(sign(Velocity-0.1*PeakVelocity));
Cross50 = diff(sign(Velocity-0.5*PeakVelocity));
Cross90 = diff(sign(Velocity-0.9*PeakVelocity));

[~,r10] = max(Cross10);
[~,r50] = max(Cross50);
[~,r90] = max(Cross90);

Rise10 = interp1(Velocity((r10):(r10+1)),Time((r10):(r10+1)),0.1*PeakVelocity,'linear');
Rise50 = interp1(Velocity((r50):(r50+1)),Time((r50):(r50+1)),0.5*PeakVelocity,'linear');
Rise90 = interp1(Velocity((r90):(r90+1)),Time((r90):(r90+1)),0.9*PeakVelocity,'linear');


end





























