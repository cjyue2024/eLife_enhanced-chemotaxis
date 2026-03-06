% 坐标变换
function Bacpos = CorrdinateTransformations(Bacpos,k,rotation)
beta = atan(k);  % k为斜率

if (strcmp(class(Bacpos),'double') & strcmp(rotation,'CW'))
    Bacpos(:,3)=Bacpos(:,3)*cos(beta)+Bacpos(:,4)*sin(beta);
    Bacpos(:,4)=Bacpos(:,4)*cos(beta)-Bacpos(:,3)*sin(beta);
    Bacpos(:,3)=Bacpos(:,3);
    Bacpos(:,4)=Bacpos(:,4);
end

if (strcmp(class(Bacpos),'double') & strcmp(rotation,'CCW'))
    Bacpos(:,3)=Bacpos(:,3)*cos(beta)-Bacpos(:,4)*sin(beta);
    Bacpos(:,4)=Bacpos(:,4)*cos(beta)+Bacpos(:,3)*sin(beta);
    Bacpos(:,3)=Bacpos(:,3);
    Bacpos(:,4)=Bacpos(:,4);
end

if (strcmp(class(Bacpos),'cell') & strcmp(rotation,'CW'))
    for i = 1:length(Bacpos)
        Bacpos{i}(:,2)=Bacpos{i}(:,2)*cos(beta)+Bacpos{i}(:,3)*sin(beta);
        Bacpos{i}(:,3)=Bacpos{i}(:,3)*cos(beta)-Bacpos{i}(:,2)*sin(beta);
        Bacpos{i}(:,2)=Bacpos{i}(:,2);
        Bacpos{i}(:,3)=Bacpos{i}(:,3);
    end
end

if (strcmp(class(Bacpos),'cell') & strcmp(rotation,'CCW'))
    for i = 1:length(Bacpos)
        Bacpos{i}(:,2) = Bacpos{i}(:,2)*cos(beta)-Bacpos{i}(:,3)*sin(beta);
        Bacpos{i}(:,3) = Bacpos{i}(:,3)*cos(beta)+Bacpos{i}(:,2)*sin(beta);
        Bacpos{i}(:,2)=Bacpos{i}(:,2);
        Bacpos{i}(:,3)=Bacpos{i}(:,3);
    end
end

end