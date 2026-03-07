% 计算每条轨迹y位置的平均值，将该条轨迹划分至某一区域（LSW/MS/RSW）
clear all
Files = dir('*um.mat');
for II = 1:length(Files)
    load(Files(II).name,'bacpos','c');
    n_top = 1;
    n_normal = 1;
    n_bottom = 1;
    
    for i = 1:length(bacpos)
        y = mean(bacpos{i}(:,3));
        if abs(y-c(1))<=3
            pos_top{n_top} = bacpos{i};
            n_top = n_top +1;
        end
        
        if abs(y-c(2))<=3
            pos_bottom{n_bottom} = bacpos{i};
            n_bottom = n_bottom +1;
        end
        
        if abs(y-c(1))>3 && abs(y-c(2))>3
            pos_normal{n_normal} = bacpos{i};
            n_normal = n_normal +1;
        end
    end
    pos_top = pos_top';
    if exist('pos_normal','var')
        pos_normal = pos_normal';
    end
    pos_bottom = pos_bottom';
    pos_side = [pos_top;pos_bottom];
    save(Files(II).name,'bacpos','c','pos_*');
    clearvars -except Files II
    
end




