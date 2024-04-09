clear;

step = 5;
angles = 0:step:355;
num_angles = length(angles);
ir_data = cell(1, num_angles);
fft_data = cell(1, num_angles);

polar_freq = [500 1000 2000 4000 8000];
num_polar_freq = length(polar_freq);
polars = zeros(num_polar_freq, num_angles);

monName = ["monCQGLL74L" "monG7SMCCBW" "monKD8N255G" "monGHHX" "mon99PJ"];
channel = 1;
for k = 1:5
    for i = 1:num_angles
        %filename = sprintf("Microphone_Impulse_Responses/ShureSM58_125cm_Normalised_IRs/IRs/ShureSM58_125cm_%dDeg.wav", angles(i));
        filename = sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Hann/IR_%s_%d_Channel_%d.wav", channel, monName(k), angles(i), channel);
        [ir_data{i}, fs] = audioread(filename);
        N = length(ir_data{i});
        fft_data{i} = fft(ir_data{i});
        magnitude = abs(fft_data{1, i}(1:N/2));
        for j = 1:num_polar_freq
            index = 1+mod(i-1+90/step,360/step);
            polars(j, index) = polars(j, index)+magnitude(polar_freq(j)*N/fs);
        end
    end
end
polars = polars/5;
polars = mag2db(polars);

figure;
angles = [angles 360];
polars = [polars polars(:,1)];
tbl = array2table([deg2rad(angles') polars']);
tbl = renamevars(tbl, "Var1", "Angles (rad)");
polar_freq_labels = string(polar_freq);
for i = 1:num_polar_freq
    polar_freq_labels(i) = sprintf("%dHz",polar_freq(i));
    tbl = renamevars(tbl, sprintf("Var%d",i+1), sprintf("%dHz",polar_freq(i)));
end

polarplot(tbl, "Angles (rad)", polar_freq_labels, 'Linewidth', 1);
legend;
thetalim([0 360]);
thetaticks(0:45:315);
rlim([-40 20]);
rticks(-40:5:20);
title(sprintf("Average Directivity Pattern for Channel %d", channel));

save(sprintf("polars_average_channel_%d", channel), "polars", "polar_freq");