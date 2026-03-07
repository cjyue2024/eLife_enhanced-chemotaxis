function tumbles = FindTumble_mode(bacpos,id,filterSelect,CalvSelect,a,f,judgeNum)
t = bacpos{id}(:,1);
x = bacpos{id}(:,2);
y = bacpos{id}(:,3);

if filterSelect == 1
    x = smooth(x);
    y = smooth(y);
end

if CalvSelect == 1
    vx = zeros(1,length(x));
    vy = zeros(1,length(y));
    for j = 3:length(x)-2
        vx(j) = (8*(x(j+1)-x(j-1))-(x(j+2)-x(j-2)))/12/0.05;     % 用中心微分方程计算细菌在j点的速度
        vy(j) = (8*(y(j+1)-y(j-1))-(y(j+2)-y(j-2)))/12/0.05;
    end
    v = sqrt(vx.^2+vy.^2);
    dv = [0,diff(v)];
end

if CalvSelect == 2
    vx = [0; diff(x)./0.05];
    vy = [0; diff(y)./0.05];
    v = sqrt(vx.^2+vy.^2);
    dv = [0;diff(v)];
end

speed_drop_threshold = mean(v)*a;
if judgeNum == 2
    tumble_indices = find(dv <0 & abs(dv) > speed_drop_threshold);
else if judgeNum == 3
        tumble_indices = find(dv <0 & abs(dv) > speed_drop_threshold & v<mean(v));
    end
end


% 提取tumble事件
if length(tumble_indices)~=0
    tumble_indices = Seg(tumble_indices);
    for i = 1:length(tumble_indices)
        if length(tumble_indices{i})<f || min(tumble_indices{i})==1 % 第二个判断：要识别前一个run
            tumble_indices{i} = [];
        end
    end
    % 删除空矩阵
    tumble_indices = tumble_indices(~cellfun('isempty', tumble_indices));
    % tumble_indices = cell2mat(tumble_indices);
end

tumbles = struct();
for kk = 1:length(tumble_indices)
    idx = tumble_indices{kk};
    t_start = t(min(idx));
    t_end = t(max(idx));
    
    y_start = y(min(idx));
    
    tumbles(kk).start_time = t_start;
    tumbles(kk).end_time = t_end;
    tumbles(kk).duration = t_end - t_start;
    tumbles(kk).start_speed = v(min(idx));
    tumbles(kk).end_speed = v(max(idx));
    tumbles(kk).y = y_start;
    
    % 计算run angle
    last_run_end = [vx(min(idx)-1),vy(min(idx)-1)];
    if max(idx) == length(vx)
        next_run_start = [vx(max(idx)),vy(max(idx))]; % 轨迹的末尾
    else
        next_run_start = [vx(max(idx)+1),vy(max(idx)+1)]; % 轨迹的末尾
    end
    dot_product = dot( last_run_end, next_run_start);
   
    
    magnitude_a = norm(last_run_end);   % 计算向量的模
    magnitude_b = norm(next_run_start);
    theta_radians = acos(dot_product / (magnitude_a * magnitude_b));    % 计算夹角（以弧度表示）
    tumbles(kk).change_in_direction = theta_radians;
end
end