load gesturenamesF

nParti=39;

for smoothing = 1:2
    
for resampling=1:2; % 1-a passi di arc length regolari o 2-di tempo regolari

centering=1; % 1 centroide -   2 origine

scaling=1; % 1 su lunghezza 2 bounding box(protractor 3D), 3 max distanza origine

comptype=1; % main trajectory comparison:  1 squared Euclid 2-absolute(cityblock)  3-DTW  4- min point to point
rotation=0;



cofi = [ 1 0 1 1 1 1 0 0 0 0 0 0 0 0];

train = readtable('datasetShrec/train_gestures.txt','Delimiter',' ','ReadVariableNames',false,'Format','%f %f %f %f %f %f %f ');
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
   
   
    labelTrain(k)=trainSet(i,5);%idGesture;
    gestotxt = sprintf('datasetShrec/gesture_%i/finger_%i/subject_%i/essai_%i/skeletons_world.txt',idGesture,idFinger,idSubject,idEssay);
    gesto=readtable(gestotxt);
    gesto=table2array(gesto);
    
    if(smoothing==2)
        for(u=1:66)
        gesto(:,u)=smooth( gesto(:,u),'lowess');
    end
    end
      
    gesto=gesto(1:ceil(size(gesto)*0.8),:);

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
    
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[7:9]);
    v11 = dirVec(m1,m2); % direz base indice.polso 2
    
    
    w1=spline(ti,v1',0:h:ti(end)); %vettore direz poll polso
    w2=spline(ti,v2',0:h:ti(end)); %vettore drez indice polso
    w3=spline(ti,v3',0:h:ti(end)); %vettore  palmo polso   
    w6=spline(ti,v6',0:h:ti(end)); %normale palmo
    w5=spline(ti,v5',0:h:ti(end)); % base indice-mignolo  
    w4=spline(ti,v4',0:h:ti(end));
    w7=spline(ti,v7',0:h:ti(end));
    w8=spline(ti,v8',0:h:ti(end)); 
    w9=spline(ti,v9',0:h:ti(end));
    w10=spline(ti,v10',0:h:ti(end));   
    w11=spline(ti,v11',0:h:ti(end));    
  %  speed = spline(ti,vel,0:h:ti(end));

    
    a1 = real(acos( sum(w1.*w2)));  % angolo pollicepolso-indicepolso
    a2 = real(acos( sum(w1.*w3)));  % angolo pollicepolso-palmopolso
    a3 = real(acos( sum(w1.*w6)));  % angolo pollicepolso-perpendicolare palmo
    a4 = real(acos( sum(w1.*w7)));  % angolo pollicepolso-palmopollice
    a5 = real(acos( sum(w2.*w3)));   % angolo indicepolso-palmopolso    
    a6 = real(acos( sum(w2.*w6)));  % angolo indicepolso-perpendicolare palmo
    a7 = real(acos( sum(w2.*w8)));   % angolo indicepolso-palmoindice
    a8 = real(acos( sum(w6.*w7))); % angolo palmo pollice-perpendicolare palmo
    a9 = real(acos( sum(w6.*w8)));  % angolo palmo indice-perpendicolare palmo
    a10 = real(acos( sum(w3.*w9))); % angolo pollice -perpendicolare palmo
    a11 = real(acos( sum(w3.*w10)));  % angolo  indice-perpendicolare palmo
    a12 = real(acos( sum(w9.*w10)));  % angolo  indice-pollice
   
    
    ww{1}=w1; ww{2}=w2; ww{3}=w3; ww{4}=w6; ww{5}=w7; ww{6}=w8; ww{7}=w9; ww{8}=w10;
     trainF2{k} =  [ ];
    count=1;
    for kk=1:8
        for ll=1:kk;
             trainF2{k} =  [trainF2{k}; real(acos( sum(ww{kk}.*ww{ll})));];
        count=count+1;
    end
    end

       
    trainFinger{k} =  [w6; w5; w1; w7; w8; w2];
    trainFinger2{k} = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10;a11; a12];
    trainFall(k,:)= [ trainFinger2{k}(:)'];    

% processa feature traiettoria (28-30 punta indice)

    gesto=gesto(:,[4:6]);%[28:30 4:6 1:3 16:18  40:42 52:54]);
          
    
    [ges ia ic]=unique(gesto(:,1:3),'rows','stable');
    gesto=gesto(ia,:);

    time=0:size(gesto,1)-1;
    time=time';
    gesto=[gesto time]; 
 
   [ trainGest{k} trainL{k} trainV{k} ]= normalizeP(gesto,nParti,resampling,centering,scaling);
 

      k=k+1;
end

test = readtable('datasetShrec/test_gestures.txt','Delimiter',' ','ReadVariableNames',false,'Format','%f %f %f %f %f %f %f ');
testSet=table2array(test(:,1:end));
testGest=cell(1,size(testSet,1));

k=1;
for i=1:size(testSet,1)
    
    idGesture=testSet(i,1);
    idFinger=testSet(i,2);
    idSubject=testSet(i,3);
    idEssay=testSet(i,4);
    label14Te(i)=testSet(i,5);
    label28Te(i)=testSet(i,6);
    labelTest(k)=testSet(i,5);%idGesture;
    gestotxt = sprintf('datasetShrec/gesture_%i/finger_%i/subject_%i/essai_%i/skeletons_world.txt',idGesture,idFinger,idSubject,idEssay);
    gesto=readtable(gestotxt);
    gesto=table2array(gesto);
    
    dist = sqrt((gesto(2:end,1)-gesto(1:end-1,1)).^2 + (gesto(2:end,2)-gesto(1:end-1,2)).^2  + (gesto(2:end,3)-gesto(1:end-1,3)).^2)';
    s=cumsum([0 dist]);
    vel = [0 dist]'/s(end);
    
    if(smoothing==2)
        for(u=1:66)
        gesto(:,u)=smooth( gesto(:,u),'lowess');
    end
    end
    gesto=gesto(1:ceil(size(gesto)*0.2),:);
    
    ti = 0:1:(size(gesto)-1);

    h=ti(end)/nParti; %passo tra le distanze

    po{1}=gesto(:,[1:3]);% polso
    po{2}=gesto(:,[4:6]);% palmo
    po{3}=gesto(:,[7:9]);% base pollice
    po{4}=gesto(:,[16:18]);%punta pollice
    po{5}=gesto(:,[19:21]);% base indice
    po{6}=gesto(:,[28:30]);%punta indice 
    po{7}=gesto(:,[55:57]);%base mignolo  
    
    testF{k} = [];
    ii=1;
    
    for kk=1:7
        for ll=1:kk-1
            dirV{ii} = dirVec(po{kk},po{ll});
             diV{ii}=spline(ti,dirV{ii}',0:h:ti(end));
             testF{k} =  [testF{k}; diV{ii}(:)];
            ii=ii+1;      
        end
    end
    
    v4=dirV{8};
    v5=dirV{20};
    v6 = [v4(:,2).*v5(:,3)-v4(:,3).*v5(:,2), v4(:,3).*v5(:,1)-v4(:,1).*v5(:,3), v4(:,1).*v5(:,2)-v4(:,2).*v5(:,1)];
    norm = sqrt(v6(:,1).*v6(:,1)+v6(:,2).*v6(:,2)+v6(:,3).*v6(:,3));
    dirV{ii}=v6./norm; % normale al palmo
    diV{ii}=spline(ti,dirV{ii}',0:h:ti(end));
    
    testF{k} =  [testF{k}; diV{ii}(:)];
    
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[16:18]);
    v1 = dirVec(m1,m2); % direz poll polso
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[28:30]);
    v2 = dirVec(m1,m2); % direz indice polso
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[4:6]);
    v3 = dirVec(m1,m2); % direz palmo polso
    
    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[19:21]);
    v4 = dirVec(m1,m2); % direz indice.polso 2
        
    m1 = gesto(:,[19:21]);
    m2 = gesto(:,[55:57]);
    v5 = dirVec(m1,m2); % direz base indice base mignolo
    
    v6= [v4(:,2).*v5(:,3)-v4(:,3).*v5(:,2), v4(:,3).*v5(:,1)-v4(:,1).*v5(:,3), v4(:,1).*v5(:,2)-v4(:,2).*v5(:,1)];
    norm = sqrt(v6(:,1).*v6(:,1)+v6(:,2).*v6(:,2)+v6(:,3).*v6(:,3));
    v6=v6./norm;
    
    m1 = gesto(:,[4:6]);
    m2 = gesto(:,[16:18]);
    v7 = dirVec(m1,m2); % direz palmo  poll
        
    m1 = gesto(:,[4:6]);
    m2 = gesto(:,[28:30]);
    v8 = dirVec(m1,m2); % direz palmo  indice
    
    m1 = gesto(:,[19:21]);
    m2 = gesto(:,[28:30]);
    v9 = dirVec(m1,m2); % direz base indice punta indice
       
    m1 = gesto(:,[7:9]);
    m2 = gesto(:,[16:18]);
    v10 = dirVec(m1,m2); % direz base pollice punta pollice

    m1 = gesto(:,[1:3]);
    m2 = gesto(:,[7:9]);
    v11 = dirVec(m1,m2); % direz base indice.polso 2
    
    w1=spline(ti,v1',0:h:ti(end)); %vettore che id
    w2=spline(ti,v2',0:h:ti(end)); %vettore che id
    w3=spline(ti,v3',0:h:ti(end)); %vettore che id
    
    w6=spline(ti,v6',0:h:ti(end));
    w5=spline(ti,v5',0:h:ti(end));   
    
    w7=spline(ti,v7',0:h:ti(end));
    w8=spline(ti,v8',0:h:ti(end));  
    w4=spline(ti,v4',0:h:ti(end));
    w9=spline(ti,v9',0:h:ti(end));
    w10=spline(ti,v10',0:h:ti(end));   
    w11=spline(ti,v11',0:h:ti(end));   
   % speed = spline(ti,vel,0:h:ti(end));   

    a1 = real(acos( sum(w1.*w2)));  % angolo pollicepolso-indicepolso
    a2 = real(acos(sum(w1.*w3)));  % angolo pollicepolso-palmopolso
    a3 = real(acos( sum(w1.*w6)));  % angolo pollicepolso-perpendicolare palmo
    a4 = real(acos( sum(w1.*w7)));  % angolo pollicepolso-palmopollice
    a5 = real(acos( sum(w2.*w3)));   % angolo indicepolso-palmopolso    
    a6 = real(acos( sum(w2.*w6)));  % angolo indicepolso-perpendicolare palmo
    a7 = real(acos( sum(w2.*w8)));   % angolo indicepolso-palmoindice
    a8 = real(acos( sum(w6.*w7))); % angolo palmo pollice-perpendicolare palmo
    a9 = real(acos( sum(w6.*w8)));  % angolo palmo indice-perpendicolare palmo
    a10 = real(acos( sum(w3.*w9))); % angolo pollice -perpendicolare palmo
    a11 = real(acos( sum(w3.*w10)));  % angolo  indice-perpendicolare palmo
    a12 = real(acos( sum(w9.*w10)));  % angolo  indice-pollice
    
    ww{1}=w1;ww{2}=w2;ww{3}=w3; ww{4}=w6; ww{5}=w7; ww{6}=w8; ww{7}=w9; ww{8}=w10;
     testF2{k} =  [];
    count=1;
    for kk=1:8
        for ll=1:kk;
             testF2{k} =  [testF2{k}; real(acos( sum(ww{kk}.*ww{ll})));];
        count=count+1;
    end
    end
 
       
   
    testFinger{k} = [w6; w5; w1; w7; w8; w2];
    testFinger2{k} = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12];
    testFall(k,:)= [ testFinger2{k}(:)'];

%palm wrist thumbtip indextip middletip ringtip
     gesto=gesto(:,[4:6]);%(:,[28:30 4:6 1:3 16:18  40:42 52:54]);
     
    [ges ia ic]=unique(gesto(:,1:3),'rows','stable');
    gesto=gesto(ia,:);
    
    
    [testGest{k} testL{k} testV{k}] = normalizeP(gesto,nParti,resampling,centering,scaling);
    k=k+1;
end



tic
for i=1:size(trainGest,2)
    for j=1:size(testGest,2)
        DT(i,j)=gestureDistance(trainGest{i}(1:end,:),testGest{j}(1:end,:),comptype,rotation);
    end
end
time1=toc

for i=1:size(trainFinger,2)
 aa{i}=reshape(trainF{i},[66 (nParti+1)]);
 aall=[aall aa{i}];
end
%[trainPCA, mapping] = compute_mapping(aall','PCA',20);
bball=[];
for i=1:size(testFinger,2)
 bb{i}=reshape(testF{i},[66 (nParti+1)]);
 bball=[bball bb{i}];
end

tic
for i=1:size(trainFinger,2)
    for j=1:size(testFinger,2) -
        DF(i,j)=gestureDist2(aa{i},bb{j},1);
     end
end

time2=toc


clear aa bb
for i=1:size(trainFinger,2)
 aa{i}=reshape(trainF2{i},[36 (nParti+1)]);
end
for i=1:size(testFinger,2)
 bb{i}=reshape(testF2{i},[36 (nParti+1)]);
end

tic
for i=1:size(trainFinger,2)
    for j=1:size(testFinger,2)      
        DF2(i,j)=gestureDist2(aa{i},bb{j},1);
     end
end
time3=toc

D = DF.*DF2.*sqrt(DT);

[val ind] = min(D);
est=label14(ind);

NN(smoothing, resampling) =  sum(est==label14Te)/size(label14Te,2)
corrl=est==label14Te;
for i=1:14
    err(i)=sum(corrl((i-1)*60+1:(i)*60))/60;
end

[Y,I] = sort(D);
k=3;
est2 = mode(label14(I(1:k,:)));
KNN3(smoothing, resampling) =  sum(est2==label14Te)/size(label14Te,2)
k=5;
est2 = mode(label14(I(1:k,:)));
KNN5(smoothing, resampling) =  sum(est2==label14Te)/size(label14Te,2)

corr=(est2==label14Te);
for label=1:14
   pc(label) = sum(corr(label14Te==label))/sum(label14Te==label);
end

conma=zeros(14,14);

for ii=1:size(label14Te,2)
   
   conma(est(ii),label14Te(ii))= conma(est(ii),label14Te(ii))+1;

end

pclass(smoothing,resampling,label)=pc(label);

end
end


