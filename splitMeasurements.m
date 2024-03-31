for i = 1:6
    for j = 0:5:115
        fileName = sprintf("AnechoicRoomMeasurements/RecycleBin2/monG7SMCCBW_%d_%d.wav", j, i);
        [ss_data, fs] = audioread(fileName);
        ss_data_1 = ss_data((0*fs+1):1:18*fs);
        audiowrite(sprintf("AnechoicRoomMeasurements/Temp/monG7SMCCBW_%d_%d_seg.wav", j, i),ss_data_1, fs);
        ss_data_2 = ss_data((18*fs+1):1:36*fs);
        audiowrite(sprintf("AnechoicRoomMeasurements/Temp/monG7SMCCBW_%d_%d_seg.wav", (j+120), i),ss_data_2, fs);
        ss_data_3 = ss_data((36*fs+1):1:52*fs);
        audiowrite(sprintf("AnechoicRoomMeasurements/Temp/monG7SMCCBW_%d_%d_seg.wav", (j+240), i),ss_data_3, fs);
    end
end