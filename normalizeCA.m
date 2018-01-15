function [tra len vel] = normalizeCA(tab,nParti,resampling,centering,scaling)
%NORMALIZZACAGLIARITANO Summary of this function goes here
%   Detailed explanation goes here

%estraggo coordinate gesto
coords=tab(:,1:3);
time=tab(:,4);

dist = sqrt((coords(2:end,1)-coords(1:end-1,1)).^2 + (coords(2:end,2)-coords(1:end-1,2)).^2  + (coords(2:end,3)-coords(1:end-1,3)).^2)';
s=cumsum([0 dist]);
vel = [0 dist]';
len = s(end);




if resampling==1
h=s(end)/nParti;
x=spline(s,tab',0:h:s(end));
else

h=time(end)/nParti; %passo tra le distanze
x=spline(time,tab',0:h:time(end)); %vettore che id
end

% plot3(tab(:,1)',tab(:,2)',tab(:,3)')
% hold on 
% plot3(x(1,:),x(2,:),x(3,:))


if centering==1
%centro su centroide
x(1:3,:) = x(1:3,:) - mean(x(1:3,:),2);
else
% centro orig
x(1:3,:) = x(1:3,:) - x(1:3,1);%
end

if scaling ==1
%scalo su lunghezza
tra(:,1:3) = x(1:3,:)'/s(end);

elseif scaling==2
% %scalatura BB
ma = max(x(1:3,:)');
mi = min(x(1:3,:)');
tra(:,1:3) =(x(1:3,1:1:end)'-mi)./(ma-mi);
elseif scaling==3
% %scalo dist origine
y(1:3,:) = x(1:3,:) - x(1:3,1);%
d= sqrt(y(1,:).*y(1,:)+y(2,:).*y(2,:)+y(3,:).*y(3,:));
md=max(d);
%scalo
tra(:,1:3) = x(1:3,:)'/md;
end
    

end


