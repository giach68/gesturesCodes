load gesturenamesF

nParti=39;

smoothing = 2
resampling=2

centering=1; % 1 centroide -   2 origine

scaling=1; % 1 su lunghezza 2 bounding box(protractor 3D), 3 max distanza origine

comptype=1; % main trajectory comparison:  1 squared Euclid 2-absolute(cityblock)  3-DTW  4- min point to point
rotation=0;



cofi = [ 1 0 1 1 1 1 0 0 0 0 0 0 0 0];

train = readtable('datasetFrancese/train_gestures.txt','Delimiter',' ','ReadVariableNames',false,'Format','%f %f %f %f %f %f %f ');
trainSet=table2array(train(:,1:end));
trainGest=cell(1,size(trainSet,1));

%Coordinate da estrarre dai dati 
%legenda
%palm 4:6 wrist  1:3  thumbtip 16:18  indextip 28:30 middletip 40:42 ringtip 52:54


trainGest=cell(1,size(trainSet,1));
k=1;

for i=1:size(trainSet,1)
    idGesture=trainSet(i,1);
    idFinger=trainSet(i,2);
    idSubject=trainSet(i,3);
    idEssay=trainSet(i,4);
    label14(i)=trainSet(i,5);
    label28(i)=trainSet(i,6);
   
    i
    labelTrain(k)=trainSet(i,5);%idGesture;
    gestotxt = sprintf('datasetFrancese/gesture_%i/finger_%i/subject_%i/essai_%i/skeletons_world.txt',idGesture,idFinger,idSubject,idEssay);
    gesto=readtable(gestotxt);
    gesto=table2array(gesto);
    
    if(smoothing==2)
        for(u=1:66)
        gesto(:,u)=smooth( gesto(:,u),'lowess');
    end
    end

% analisi feature dita
    ti = 0:1:(size(gesto)-1);
    dist = sqrt((gesto(2:end,1)-gesto(1:end-1,1)).^2 + (gesto(2:end,2)-gesto(1:end-1,2)).^2  + (gesto(2:end,3)-gesto(1:end-1,3)).^2)';
    s=cumsum([0 dist]);
    vel = [0 dist]'/s(end);

    h=ti(end)/nParti; %passo tra le distanze
    x=spline(ti,gesto',0:h:ti(end))'; %vettore che id

    po{1}=gesto(:,[1:3]);% polso
    po{2}=gesto(:,[4:6]);% palmo
    po{3}=gesto(:,[7:9]);% base pollice
    po{4}=gesto(:,[16:18]);%punta pollice
    po{5}=gesto(:,[19:21]);% base indice
    po{6}=gesto(:,[28:30]);%punta indice 
    po{7}=gesto(:,[55:57]);%base mignolo  
    trainF{k} =  [];
    ii=1;
    for kk=1:7
        for ll=1:kk-1
            dirV{ii} = dirVec(po{kk},po{ll});
             diV{ii}=spline(ti,dirV{ii}',0:h:ti(end));
             trainF{k} =  [trainF{k}(:); diV{ii}(:)];
            ii=ii+1;
          
        end
    end
    v4=dirV{8};
    v5=dirV{20};
    v6 = [v4(:,2).*v5(:,3)-v4(:,3).*v5(:,2), v4(:,3).*v5(:,1)-v4(:,1).*v5(:,3), v4(:,1).*v5(:,2)-v4(:,2).*v5(:,1)];
    norm = sqrt(v6(:,1).*v6(:,1)+v6(:,2).*v6(:,2)+v6(:,3).*v6(:,3));
    dirV{ii}=v6./norm; % normale al palmo
    diV{ii}=spline(ti,dirV{ii}',0:h:ti(end));
    trainF{k} =  [trainF{k}; diV{ii}(:)];
    
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[16:18]);
    v1 = dirVec(m1,m2); % direz poll polso
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[28:30]);
    v2 = dirVec(m1,m2); % direz indice polso
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[4:6]);
    v3 = dirVec(m1,m2); % direz palmo polso
%     a1 = acos( sum(v1.*v2,2));    
%     his1 = histcounts(a1,[0:0.2:1])/numel(a1);
%     
%     a2 = acos( sum(v1.*v3,2));
%     his2 = histcounts(a2,[0:0.2:1])/numel(a2);
    
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[19:21]);
    v4 = dirVec(m1,m2); % direz base indice.polso 2
            
    m1 = gesto(:,[19:21]);
    m2 = gesto(:,[55:57]);
    v5 = dirVec(m1,m2); % direz base indice base mignolo
    
    m1 = gesto(:,[4:6]);
    m2 = gesto(:,[16:18]);
    v7 = dirVec(m1,m2); % direz palmo  poll
       
    m1 = gesto(:,[4:6]);
    m2 = gesto(:,[28:30]);
    v8 = dirVec(m1,m2); % direz palmo  indice
       
    
    v6= [v4(:,2).*v5(:,3)-v4(:,3).*v5(:,2), v4(:,3).*v5(:,1)-v4(:,1).*v5(:,3), v4(:,1).*v5(:,2)-v4(:,2).*v5(:,1)];
    norm = sqrt(v6(:,1).*v6(:,1)+v6(:,2).*v6(:,2)+v6(:,3).*v6(:,3));
    v6=v6./norm; % normale al palmo
    
    m1 = gesto(:,[19:21]);
    m2 = gesto(:,[28:30]);
    v9 = dirVec(m1,m2); % direz base indice punta indice
       
    m1 = gesto(:,[7:9]);
    m2 = gesto(:,[16:18]);
    v10 = dirVec(m1,m2); % direz base pollice punta pollice
    
    w1=spline(ti,v1',0:h:ti(end)); %vettore direz poll polso
    w2=spline(ti,v2',0:h:ti(end)); %vettore drez indice polso
    w3=spline(ti,v3',0:h:ti(end)); %vettore  palmo polso   
    w6=spline(ti,v6',0:h:ti(end)); %normale palmo
    w5=spline(ti,v5',0:h:ti(end)); % base indice-mignolo   
    w7=spline(ti,v7',0:h:ti(end));
    w8=spline(ti,v8',0:h:ti(end)); 
    w9=spline(ti,v9',0:h:ti(end));
    w10=spline(ti,v10',0:h:ti(end));   
        
    speed = spline(ti,vel,0:h:ti(end));

%     a1 = sum(w1.*w2);  % angolo pollicepolso-indicepolso
%     a2 = sum(w1.*w3);  % angolo pollicepolso-palmopolso
%     a3 = sum(w1.*w6);  % angolo pollicepolso-perpendicolare palmo
%     a4 = sum(w1.*w7);  % angolo pollicepolso-palmopollice
%     a5 = sum(w2.*w3);  % angolo indicepolso-palmopolso    
%     a6 = sum(w2.*w6);  % angolo indicepolso-perpendicolare palmo
%     a7 = sum(w2.*w8);   % angolo indicepolso-palmoindice
%     a8 = sum(w6.*w7); % angolo palmo pollice-perpendicolare palmo
%     a9 = sum(w6.*w8);  % angolo palmo indice-perpendicolare palmo
%     
    a1 = real(acos( sum(w1.*w2)));  % angolo pollicepolso-indicepolso
    a2 = real(acos( sum(w1.*w3)));  % angolo pollicepolso-palmopolso
    a3 = real(acos( sum(w1.*w6)));  % angolo pollicepolso-perpendicolare palmo
    a4 = real(acos( sum(w1.*w7)));  % angolo pollicepolso-palmopollice
    a5 = real(acos( sum(w2.*w3)));   % angolo indicepolso-palmopolso    
    a6 = real(acos( sum(w2.*w6)));  % angolo indicepolso-perpendicolare palmo
    a7 = real(acos( sum(w2.*w8)));   % angolo indicepolso-palmoindice
    a8 = real(acos( sum(w6.*w7))); % angolo palmo pollice-perpendicolare palmo
    a9 = real(acos( sum(w6.*w8)));  % angolo palmo indice-perpendicolare palmo
%     a10 = real(acos( sum(w6.*w9))); % angolo pollice -perpendicolare palmo
%     a11 = real(acos( sum(w6.*w10)));  % angolo  indice-perpendicolare palmo
    a10 = real(acos( sum(w3.*w9))); % angolo pollice -perpendicolare palmo
    a11 = real(acos( sum(w3.*w10)));  % angolo  indice-perpendicolare palmo
    a12 = real(acos( sum(w9.*w10)));  % angolo  indice-pollice
    
    trainFinger{k} =  [w6; w5; w1; w7; w8; w2];
    trainFinger2{k} = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; speed];
    trainFall(k,:)= [ trainFinger2{k}(:)'];    

% processa feature traiettoria (28-30 punta indice)

    gesto=gesto(:,[4:6]);%[28:30 4:6 1:3 16:18  40:42 52:54]);
          
    
    [ges ia ic]=unique(gesto(:,1:3),'rows','stable');
    gesto=gesto(ia,:);

    time=0:size(gesto,1)-1;
    time=time';
    gesto=[gesto time]; 
 
   [ trainGest{k} trainL{k} trainV{k} ]= normalizeP(gesto,nParti,resampling,centering,scaling);
 
%     trcf(k,1)=mean(trainV{k}(1:floor(size(trainV{k},1)/2)));
%     trcf(k,2)=mean(trainV{k}(ceil(size(trainV{k},1)/2):end));
%  

      k=k+1;
end



for i=1:size(trainGest,2)
    for j=1:size(trainGest,2)
        DT(i,j)=gestureDistance(trainGest{i}(1:end,:),trainGest{j}(1:end,:),comptype,rotation);
      
    end
end

for i=1:size(trainFinger,2)
 aa{i}=reshape(trainF{i},[66 (nParti+1)]);
end

for i=1:size(trainFinger,2)
    for j=1:size(trainFinger,2)   
        DF(i,j)=gestureDist2(trainFinger{i}(:,:),trainFinger{j}(:,:),1);
        
    end
end


for i=1:size(trainFinger,2)
    for j=1:size(trainFinger,2) 
       
        DFX(i,j)=gestureDist2(aa{i},aa{j},1);
    end
end

for esc=1:size(trainFinger2{1},1)
range=[esc];%[1:esc-1 esc+1:size(trainFinger2{1},1)];



for i=1:size(trainFinger,2)
    for j=1:size(trainFinger,2)   
        DF2(i,j)=gestureDist2(trainFinger2{i}(range,:),trainFinger2{j}(range,:),1);
    end
end


D=DFX.*sqrt(DT);
D(eye(size(D)>0))=inf;

[val ind] = min(D);
est=label14(ind);

NN =  sum(est==label14)/size(label14,2)
corrl=est==label14;


[Y,I] = sort(D);
first=label14(I(1:50,:));

FT =  mean(sum(first==squeeze(label14))/50)


end
for i=1:14
    err(i)=sum(corrl((i-1)*60+1:(i)*60))/60;
end

% DFF = pdist2(testF,trainF)';
% D = DFF.*sqrt(DT);

[Y,I] = sort(D);
k=3;
est2 = mode(label14(I(1:k,:)));
KNN3(smoothing,resampling) =  sum(est2==label14Te)/size(label14Te,2)
k=5;
est2 = mode(label14(I(1:k,:)));
KNN5(smoothing,resampling) =  sum(est2==label14Te)/size(label14Te,2)

corr=(est2==label14Te);
for label=1:14
   pc(label) = sum(corr(label14Te==label))/sum(label14Te==label);
end

pclass(smoothing,resampling,label)=pc(label);

