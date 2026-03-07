function MSOD = Cal_MSOD(bacpos)  % MSOD: mean-squared orientational displacement
for i = 1:length(bacpos)
    x = bacpos{i}(:,2);
    y = bacpos{i}(:,3);
    for j = 3:length(x)-2
        vx(j) = (8*(x(j+1)-x(j-1))-(x(j+2)-x(j-2)))/12/0.05;     % 用四阶中心微分方程计算细菌在i点的速度
        vy(j) = (8*(y(j+1)-y(j-1))-(y(j+2)-y(j-2)))/12/0.05;
        v(j) =  sqrt(vx(j).^2 + vy(j).^2);
    end
    vx(1:2) = [];
    vy(1:2) = [];
    v(1:2) = [];
    vx(end-2:end) = [];
    vy(end-2:end) = [];
    v(end-2:end) = [];
    MSD=[];
    for df = 1:20
        cos_theta = (vx(1+df:end).*vx(1:end-df) + vy(1+df:end).*vy(1:end-df))...
            ./sqrt(vx(1:end-df).^2 + vy(1:end-df).^2 )...
            ./sqrt(vx(1+df:end).^2 + vy(1+df:end).^2 );
%         MSD(df) = mean(acos(cos_theta).^2);
        MSD(df) = mean(real(acos(cos_theta)).^2);
    end
    MSOD(i,:) =  MSD;
    clearvars -except bacpos i MSOD
end
end