c = 343.3;
delay_diff = (sqrt(5)-1)/c;
fix = [3 3 3 -3 -3 -3];

for i = 1:6
    for j = 0:5:355
        fileName = sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Original/IR_mon99PJ_%d_Channel_%d.wav", i, j, i);
        [ir_data, fs] = audioread(fileName);
        [M, I] = max(ir_data);
        N_hann_half = round(fs*delay_diff);
        N_hann = 2*N_hann_half;
        w = hann(N_hann);
        N = length(ir_data);
        w_padd = [zeros(I-N_hann/2, 1);w;zeros(N-I-N_hann/2, 1)];
        ir_data_hann = ir_data.*w_padd;
        index = i + fix(i);
        audiowrite(sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Hann/IR_mon99PJ_%d_Channel_%d.wav", index, j, index), ir_data_hann, fs);
    end
end