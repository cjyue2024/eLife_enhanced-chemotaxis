% 通过直线拟合算drift velocity
function [mdt,mdx,stdmdx] = Cal_dt_dx(bacpos)
Dx = zeros(100,length(bacpos));
for i=1:length(bacpos)
    pos = bacpos{i};
    for k = 1:2:length(pos)-1
        Dx(k,i) = mean(pos(k+1:2:end,2)-pos(1:2:end-k,2));
    end
end

Dx(all(Dx==0,2),:)=[];
%
Dx(21:end,:) = [];         % 20/20*2 = 2s ， 30/20*2 = 3s
mdx = []; stdmdx =[];  mdt=[];
for i=1:size(Dx,1)
    mdx(i,1)=mean(nonzeros(Dx(i,:)));
    %     stdmdx(i,1) = std(nonzeros(Dx(i,:)),1);
    stdmdx(i,1) = std(nonzeros(Dx(i,:)),1)/sqrt(length(nonzeros(Dx(i,:))));
    mdt(i,1)=(2*i-1)*0.05;
end
% errorbar(mdt,mdx,stdmdx,'color','c','Linestyle','None');
% hold on;errorbar(mdt,mdx,stdmdx,'color',c);
end