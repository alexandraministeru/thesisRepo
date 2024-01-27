function [PitCom,intErr] = getPIcom(ref,measOutput,prevPitCom,prevIntErr,Ts)
KP = 15;
KI = 5;

% KP = 23;
% KI = 17;
% KK      = 0.1099965;
KK = 10; 

% Compute the gain scheduling correction factor based on the previously commanded pitch angle for blade 1:
GK = 1.0/( 1.0 + prevPitCom/KK );
% GK = 1;

uhat_max = deg2rad(10);
minPit  = deg2rad(-10)/uhat_max;
maxPit  = deg2rad(10)/uhat_max;
maxRate = deg2rad(5)/uhat_max;

err = measOutput - ref;
intErr = prevIntErr + err*Ts;
intErr = min( max( intErr, minPit/( GK*KI ) ), maxPit/( GK*KI )); % Saturate integral error

% Compute the pitch commands associated with the proportional and integral gains:
PitComP   = GK*KP*err;                                                            % Proportional term
PitComI   = GK*KI*intErr;                                                         % Integral term (saturated)

% Superimpose the individual commands to get the total pitch command;
%  saturate the overall command using the pitch angle limits:
PitComT   = PitComP + PitComI;                                                                        % Overall command (unsaturated)
PitComT   = min( max( PitComT, minPit ), maxPit);                  % Saturate the overall command using the pitch angle limits

% Saturate the overall commanded pitch using the pitch rate limit:
% NOTE: Since the current pitch angle may be different for each blade
% (depending on the type of actuator implemented in the structural
% dynamics model), this pitch rate limit calculation and the
% resulting overall pitch angle command may be different for each blade.
PitRate   = ( PitComT - prevPitCom)/Ts;                                                 % Pitch rate of the blades (unsaturated)
PitRate   = min( max( PitRate, -maxRate ), maxRate );                 % Saturate the pitch rate of the blade using its maximum absolute value

PitCom    = prevPitCom + PitRate*Ts;                                                    % Saturate the overall command of the blade using the pitch rate limit
PitCom    = min( max( PitCom, minPit ), maxPit );                   % Saturate the overall command using the pitch angle limits

end