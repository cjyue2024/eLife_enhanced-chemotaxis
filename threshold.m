% Thefilesrshold: Trajectories length and bacteria swimming speed
% get vars uint the same as imageJ
clear all
files = dir('Data04*.mat');
Magnification = 20;
fps = 20;
for II = 1:length(files)
    load(files(II).name);
    j = 1;
    k = [];
    for i= 1 : length(tracks)
        t = tracks{i}(:, 1);
        x = tracks{i}(:, 2)/Magnification;
        y = tracks{i}(:, 3)/Magnification;
        dx = diff(x);
        dy = diff(y);
        vtemp = (sqrt(dx.^2+dy.^2))*fps;
        Dt = t(end)-t(1);
        if (Dt>=1 && mean(vtemp)>10)
            bacpos{j,:} = tracks{i};
            j = j+1;
            k = [k;i];
        end
    end
    bacpos = Norm_to_um(bacpos,Magnification);
    save([files(II).name(1:end-4),'.mat'],'bacpos','-append');
    clearvars -except files II Magnification fps
end
