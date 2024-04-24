function CPSD = calculateCPSD(signal_1, signal_2, n_mics, fs, f_sound)
    CPSD = gpuArray.ones(n_mics, n_mics);
    for i = 1:n_mics
        for j = i:n_mics
            [cpsd_ij, f] = cpsd(signal_1(:,i), signal_2(:,j), [], [], [], fs);
            f_index = compareIndex(f_sound,fs,f);
            CPSD(i,j) = cpsd_ij(f_index);
            CPSD(j,i) = conj(CPSD(i,j));
        end
    end
end

