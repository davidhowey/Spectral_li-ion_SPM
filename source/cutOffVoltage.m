function [val,term,dir] = cutOffVoltage(t,y,I,data,matrices_spm,V_limit)
% CUTOFFVOLTAGE detects when voltage limits are exceeded in order to 
% stop time integration
% see http://uk.mathworks.com/help/matlab/math/ode-event-location.html
%
% INPUTS
% t             The time at which the cut off condition has to be
%               calculated.
% y             The state at which the cut off condition has to be
%               calculated.
% I             The applied current.
% data          Structure with parameters of the battery.
% params        Structure with parameters of the model.
% matrices_spm  Structure with the differentiation matrices.
% V_limit       Vector defining the cut off voltages [V_min V_max].
%
% OUTPUTS
% val           Value of the event. event(i) occurs if val(i) = 0.
% term          term(i) indicates if event(i) is terminal or not (i.e. if
%               time integration should be stopped if the event occurs or   
%               not).
% dir           dir(i) indicates the direction in which event (i) happens.
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

    % Calculate the voltage using the measurement function get_measurements
    [V,~] = get_measurements(t,y,I,data,matrices_spm);
    V = real(V);

    % Calculate the voltage 'window' left before the cut off voltage is
    % reached. (When the window is <= 0, the V limits are exceeded)
    Vspare_dis = V - V_limit(1); % V 'left' on discharge
    Vspare_cha = V_limit(2)-V; % V 'left' on charge

    % If we reach the voltage limits, Vspare will become 0, and an event
    % occurs
    val = [Vspare_dis ; Vspare_cha];
    % if term(i)=1, event(i) is terminal and time integration is halted
    % when it happens. 
    term = [1;1];
    % event occurs only if 0 is reached from direction(i) (-1 means 'val'
    % is decreasing and event occurs when val=0, 0 = from any direction, 1
    % means 'val' is increasing and event when val=0 (val=negative))
    dir = [-1;-1];
        % note: with current settings, an event occurs only if V is
        % decreasing and we cross the minimum voltage, or if V is
        % increasing and we cross the maximum voltage.
        % If you start charging at V < Vmin, no event will occur if
        % we cross Vmin. If you start discharging at V < Vmin, no event
        % will occur (eventually errors will be produced, e.g. if the
        % concentration becomes negative)
        % Similar for starting at V>Vmax (no event if you cross Vmax on
        % discharging, no event ever if you start charging)
end














