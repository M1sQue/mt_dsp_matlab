inverseSineSweep = audioread("SineSweep/inverse_22050.wav");
for i = 1:6
    for j = 0:5:355
        fileName = sprintf("AnechoicRoomMeasurements/Temp/monG7SMCCBW_%d_%d_seg.wav", j, i);
        [ss_seg_data, fs] = audioread(fileName);
        IR_temp = conv(ss_seg_data, inverseSineSweep);
        maxNum = max(IR_temp);
        IR = IR_temp/maxNum;
        audiowrite(sprintf("AnechoicRoomMeasurements/IRs_Channel_%d/IR_monG7SMCCBW_%d_Channel_%d.wav", i, j, i), IR(15*fs:1:length(IR)), fs);
    end
end