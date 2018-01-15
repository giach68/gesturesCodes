function [ dq ] = gestureDistance( gesture1, gesture2, mode, rot )

% per ora solo somma distanze punto punto

if rot==0
 if mode==1
dq=sum(sum((gesture1(:,1:3)-gesture2(:,1:3)).^2,2));
 elseif mode==2
dq=sum(sum(abs(gesture1(:,1:3)-gesture2(:,1:3)),2));
 elseif mode==3
     dq=dtw( gesture1(:,1:3),gesture2(:,1:3));
 elseif mode==4
     for i = 1:size(gesture1,1)
          dq = dq + min(sum(sum((gesture1(:,1:3)-gesture2(i,1:3)).^2,2)));
     end
 end
else
[dq,Z] = procrustes(gesture1(:,1:3),gesture2(:,1:3), 'Scaling',false);
dq=sum(sum((gesture1(:,1:3)-Z).^2,2));
if mode==2
    dq=sum(sum(abs(gesture1(:,1:3)-Z(:,1:3)),2));
elseif mode==3
     dq=dtw( gesture1(:,1:3),Z(:,1:3));
 elseif mode==4
     for i = 1:size(gesture1,1)
          dq = dq+ min(sum(sum((gesture1(:,1:3)-Z(i,1:3)).^2,2)));
     end
end

end


