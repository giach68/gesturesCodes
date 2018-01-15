function [ dq ] = gestureDist2( gesture1, gesture2, met)

if met==2
gesture1=gesture1-mean(gesture1);
gesture2=gesture2-mean(gesture2);
end

if met<3
dq=dtw( gesture1(:,1:end),gesture2(:,1:end));
elseif met==3
dq=dtw( gesture1(:,1:end),gesture2(:,1:end),'absolute');  
elseif met==4
dq=dtw( gesture1(:,1:end),gesture2(:,1:end),'squared');  
elseif met==5
dq=dtw( gesture1(1:3,1:end),gesture2(1:3,1:end));
dq=dq+dtw( gesture1(4:6,1:end),gesture2(4:6,1:end));
dq=dq+dtw( gesture1(7:9,1:end),gesture2(7:9,1:end));
dq=dq+dtw( gesture1(10:12,1:end),gesture2(10:12,1:end));
end
end


