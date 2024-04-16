clear;
polars_cell = cell(1,6);
for i = 1:6
    load(sprintf("MatData/polars_average_channel_%d.mat", i));
    polars_cell{i} = polars;
end

f_sound = polar_freq; % sound frequency for beamforming
% f_sound = [2000]; % sound frequency for beamforming
r = 0.057; % coordinate unit length
c = 343.3; % speed of sound

% microphone system parameters definition
m_pos = [r 0 0; r/2 -sqrt(3)/2*r 0; -r/2 -sqrt(3)/2*r 0; -r 0 0; -r/2 sqrt(3)/2*r 0; r/2 sqrt(3)/2*r 0]';
%m_pos = [0 0 0; r 0 0; 0 -r 0; -r 0 0; 0 r 0; 2*r 0 0; 2*r/2 -sqrt(3)/2*2*r 0; -2*r/2 -sqrt(3)/2*2*r 0; -2*r 0 0; -2*r/2 sqrt(3)/2*2*r 0; 2*r/2 sqrt(3)/2*2*r 0; 3*r 0 0; 0 -3*r 0; -3*r 0 0; 0 3*r 0; 2*3*r 0 0; 2*3*r/2 -sqrt(3)/2*2*3*r 0; -2*3*r/2 -sqrt(3)/2*2*3*r 0; -2*3*r 0 0; -2*3*r/2 sqrt(3)/2*2*3*r 0; 2*3*r/2 sqrt(3)/2*2*3*r 0]';
%m_pos = [-5*r 0 0; -4*r 0 0; -3*r 0 0; -2*r 0 0; -r 0 0; 0 0 0; r 0 0; 2*r 0 0; 3*r 0 0; 4*r 0 0; 5*r 0 0]';

% sound source parameters definition
azimuth_deg = 0;
elevation_deg = -90;
azimuth = deg2rad(azimuth_deg);
elevation = deg2rad(elevation_deg);
s_pos = [cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

% kaiser window beamformer algorithm
% define target center beamwidth
theta_CBW = 90; % degree
f_min = c/(numel(m_pos(1,:))*r*sin(deg2rad(theta_CBW/2)));
f_max = c/r;
%check if it is applicable to certain freq
f_sound_index = (f_sound < f_max) & (f_sound > f_min);
f_sound = f_sound(f_sound_index);

% directivity patterns of target frequencies
numFreq = numel(f_sound);
sys = cell(1,numFreq);

for i = 1:numFreq
    sys{i} = db2mag([polars_cell{1}(i, :); polars_cell{2}(i, :); polars_cell{3}(i, :); polars_cell{4}(i, :); polars_cell{5}(i, :); polars_cell{6}(i, :)]');
end


%???
A_SL =(26*numel(m_pos(1,:))*f_sound*r/c)*sin(deg2rad(theta_CBW/2))-12; % sidelobe amplitude
beta = 0.76608*(A_SL-13.26).^0.4 + 0.09834*(A_SL-13.26);

kw_delay = s_pos*m_pos/norm(s_pos)/c;
max_delay = max(kw_delay);  % delay is symetric since the reference point of mic array is in the center
kw_index = kw_delay/max_delay;  % (-1,1)
d_kw = exp(-1j*2*pi*f_sound'*kw_delay);

% weights
w_kw = besseli(0,beta'*sqrt(1-kw_index.^2))./repmat(besseli(0,beta)',1,numel(kw_index));
w_kw = w_kw/sum(sum(abs(w_kw)));
% simulation parameters
sound_delay_angles = deg2rad(0:5:360);
sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
sound_delays = -sound_delay_positions*m_pos/c;
simulations = zeros(numFreq,numel(sound_delay_angles));
for i = 1:numFreq
    simulations(i,:) = (w_kw(i,:).*d_kw(i,:))*(exp(-1j*2*pi*f_sound(i)*sound_delays).*sys{i}).';
end

% simulation polar plot
threshold = -40;
simulations_dB = 20 * log10(abs(simulations));
under_threshold_indices = simulations_dB < threshold;
simulations_dB(under_threshold_indices) = threshold;

tbl = array2table([sound_delay_angles' simulations_dB']);
tbl = renamevars(tbl, "Var1", "Angles (rad)");
polar_freq_labels = string(f_sound);
for i = 1:numFreq
    polar_freq_labels(i) = sprintf("%dHz",f_sound(i));
    tbl = renamevars(tbl, sprintf("Var%d",i+1), sprintf("%dHz",f_sound(i)));
end
figure("Name","Beam Patterns of Kaiser Window(polar plot)");
polarplot(tbl, "Angles (rad)", polar_freq_labels, 'Linewidth', 1);
legend;

thetalim([0 360]);
thetaticks(0:45:315);
rlim([threshold 0]);
rticks(threshold:5:0);
title(sprintf("Kaiser Window with 6 channel avg directivity patterns \n azimuth %d°, elevation %d°", azimuth_deg, elevation_deg));

%
tbl.("Angles (rad)") = tbl.("Angles (rad)")/pi*180;
figure("Name","Beam Patterns of Kaiser Window");
plot(tbl, "Angles (rad)", polar_freq_labels, 'Linewidth', 1);
xlim([0 360]);
title(sprintf("Kaiser Window with 6 channel avg directivity patterns \n azimuth %d°, elevation %d°", azimuth_deg, elevation_deg));
legend;
%%
clear;
polars_cell = cell(1,6);
for i = 1:6
    load(sprintf("MatData/polars_average_channel_%d.mat", i));
    polars_cell{i} = polars;
end

f_sound = polar_freq; % sound frequency for beamforming
r = 0.057; % coordinate unit length
c = 343.3; % speed of sound

% microphone system parameters definition
m_pos = [r 0 0; r/2 -sqrt(3)/2*r 0; -r/2 -sqrt(3)/2*r 0; -r 0 0; -r/2 sqrt(3)/2*r 0; r/2 sqrt(3)/2*r 0]';
%m_pos = [0 0 0; r 0 0; 0 -r 0; -r 0 0; 0 r 0; 2*r 0 0; 2*r/2 -sqrt(3)/2*2*r 0; -2*r/2 -sqrt(3)/2*2*r 0; -2*r 0 0; -2*r/2 sqrt(3)/2*2*r 0; 2*r/2 sqrt(3)/2*2*r 0; 3*r 0 0; 0 -3*r 0; -3*r 0 0; 0 3*r 0; 2*3*r 0 0; 2*3*r/2 -sqrt(3)/2*2*3*r 0; -2*3*r/2 -sqrt(3)/2*2*3*r 0; -2*3*r 0 0; -2*3*r/2 sqrt(3)/2*2*3*r 0; 2*3*r/2 sqrt(3)/2*2*3*r 0]';
%m_pos = [-5*r 0 0; -4*r 0 0; -3*r 0 0; -2*r 0 0; -r 0 0; 0 0 0; r 0 0; 2*r 0 0; 3*r 0 0; 4*r 0 0; 5*r 0 0]';

% sound source parameters definition
elevation_range = -180:1:180; % Range of elevation angles

% Animation loop
for elevation_deg = elevation_range
    azimuth_deg = 0;
    azimuth = deg2rad(azimuth_deg);
    elevation = deg2rad(elevation_deg);
    s_pos = [cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];
    
    % kaiser window beamformer algorithm
    % define target center beamwidth
    theta_CBW = 90; % degree
    f_min = c/(numel(m_pos(1,:))*r*sin(deg2rad(theta_CBW/2)));
    f_max = c/r;
    %check if it is applicable to certain freq
    f_sound_index = (f_sound < f_max) & (f_sound > f_min);
    f_sound = f_sound(f_sound_index);
    
    % directivity patterns of target frequencies
    numFreq = numel(f_sound);
    sys = cell(1,numFreq);
    
    for i = 1:numFreq
        sys{i} = db2mag([polars_cell{1}(i, :); polars_cell{2}(i, :); polars_cell{3}(i, :); polars_cell{4}(i, :); polars_cell{5}(i, :); polars_cell{6}(i, :)]');
    end
    
    
    %???
    A_SL =(26*numel(m_pos(1,:))*f_sound*r/c)*sin(deg2rad(theta_CBW/2))-12; % sidelobe amplitude
    beta = 0.76608*(A_SL-13.26).^0.4 + 0.09834*(A_SL-13.26);
    
    kw_delay = s_pos*m_pos/norm(s_pos)/c;
    max_delay = max(kw_delay);  % delay is symetric since the reference point of mic array is in the center
    kw_index = kw_delay/max_delay;  % (-1,1)
    d_kw = exp(-1j*2*pi*f_sound'*kw_delay);
    
    % weights
    w_kw = besseli(0,beta'*sqrt(1-kw_index.^2))./repmat(besseli(0,beta)',1,numel(kw_index));
    w_kw = w_kw/sum(sum(abs(w_kw)));
    % simulation parameters
    sound_delay_angles = deg2rad(0:5:360);
    sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
    sound_delays = -sound_delay_positions*m_pos/c;
    simulations = zeros(numFreq,numel(sound_delay_angles));
    for i = 1:numFreq
        simulations(i,:) = (w_kw(i,:).*d_kw(i,:))*(exp(-1j*2*pi*f_sound(i)*sound_delays).*sys{i}).';
    end
    
    % simulation polar plot
    threshold = -60;
    simulations_dB = 20 * log10(abs(simulations));
    under_threshold_indices = simulations_dB < threshold;
    simulations_dB(under_threshold_indices) = threshold;
    
    tbl = array2table([sound_delay_angles' simulations_dB']);
    tbl = renamevars(tbl, "Var1", "Angles (rad)");
    polar_freq_labels = string(f_sound);
    for i = 1:numFreq
        polar_freq_labels(i) = sprintf("%dHz",f_sound(i));
        tbl = renamevars(tbl, sprintf("Var%d",i+1), sprintf("%dHz",f_sound(i)));
    end
    polarplot(tbl, "Angles (rad)", polar_freq_labels, 'Linewidth', 1);
    legend;
    
    thetalim([0 360]);
    thetaticks(0:45:315);
    rlim([threshold 0]);
    rticks(threshold:5:0);
    title(sprintf("Kaiser Window with 6 channel avg directivity patterns \n azimuth %d°, elevation %d°", azimuth_deg, elevation_deg));
    drawnow;
end;
