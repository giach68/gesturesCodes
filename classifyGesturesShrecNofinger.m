load gesturenamesF
dp=0;
nParti=20;
smoothing=0

resampling=1; % 1-a passi di arc length regolari o 2-di tempo regolari

centering=1; % 1 centroide -   2 origine

scaling=1; % 1 su lunghezza 2 bounding box(protractor 3D), 3 max distanza origine

comptype=1; % main trajectory comparison:  1 squared Euclid 2-absolute(cityblock)  3-DTW  4- min point to point
rotation=0;


tic
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
      
% uncomment (in both training and test if smoothed)  
if smoothing
    for(u=1:66)
        gesto(:,u)=smooth( gesto(:,u),'lowess');
    end
end

%qui usa la punta indice
%indextip palm wrist thumbtip  middletip ringtip

     gesto=gesto(:,[4:6]);% 4:6 1:3 16:18  40:42 52:54]);
% if dp
%      [gesto,ix] = dpsimplify(gesto,0.05);
% end
         [ges ia ic]=unique(gesto(:,1:3),'rows','stable');
    gesto=gesto(ia,:);
    
    %opzioni 
    % resampling 1- su lunghezza  2- su tempo
    % centering 1- centroide  2-origine
    % scaling  1-lunghezza 2-bounding box(protractor) - 3-max dist orig
    
    

   [ trainGest{k} trainL{k} trainV{k} ]= normalizeP(gesto,nParti,resampling,centering,scaling);

   
    
%     if(i==1 && idGesture==1) 
%         plot3(trainGest{k}(:,1),trainGest{k}(:,2),trainGest{k}(:,3))
%         hold on
%     end
%      if(i>1 && idGesture==1) 
%         plot3(trainGest{k}(:,1),trainGest{k}(:,2),trainGest{k}(:,3))
%       waitforbuttonpress
%      end
    
      k=k+1;
end
normtime=toc

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
    
%        
if smoothing
    for(u=1:66)
        gesto(:,1)=smooth( gesto(:,u),'lowess');
    end
end   
 

%indextip palm wrist thumbtip middletip ringtip
    gesto=gesto(:,[4:6]);% 4:6 1:3 16:18  40:42 52:54]);
%   if dp
%      [gesto,ix] = dpsimplify(gesto,0.05);
% end   
    [ges ia ic]=unique(gesto(:,1:3),'rows','stable');
    gesto=gesto(ia,:);
    
   
    [testGest{k} testL{k} testV{k}] = normalizeP(gesto,nParti,resampling,centering,scaling);

    k=k+1;
end
time(k)=toc;


% opzioni 
% mode 1- euclidean (squared) 2-absolute  3-DTW
% rot 0 non rotated - 1 rotated
% 

for i=1:size(trainGest,2)
    for j=1:size(testGest,2)
        D(i,j)=gestureDistance(trainGest{i},testGest{j},comptype,rotation);
    end
end
timedist=toc


[val ind] = min(D);


est=label14(ind);

NN =  sum(est==label14Te)/size(label14Te,2)
corrl=est==label14Te;
for i=1:14
    err(i)=sum(corrl((i-1)*60+1:(i)*60))/60;
end


[Y,I] = sort(D);
k=3;
est2 = mode(label14(I(1:k,:)));
KNN =  sum(est2==label14Te)/size(label14Te,2)
k=5;
est2 = mode(label14(I(1:k,:)));
KNN =  sum(est2==label14Te)/size(label14Te,2)

corr=(est2==label14Te);
for label=1:14
   pc(label) = sum(corr(label14Te==label))/sum(label14Te==label);
end

pc


% 
% %etichetto tutti i fine uguali a 3
% 
labelRid=label14;
labelRidTe=label14Te;
labelRid(labelRid==1)=3;
labelRid(labelRid==4)=3;
labelRid(labelRid==5)=3;
labelRid(labelRid==6)=3;
labelRidTe(labelRidTe==1)=3;
labelRidTe(labelRidTe==4)=3;
labelRidTe(labelRidTe==5)=3;
labelRidTe(labelRidTe==6)=3;

labRTr=labelRid(labelRid~=3);
labRTe=labelRidTe(labelRidTe~=3);

ii=labelRid ~=3;
jj=labelRidTe ~=3;
DR=D(ii,jj);

[valr indr] = min(DR);

estr=labRTr(indr);

NNR =  sum(estr==labRTe)/size(labRTe,2)

[Y,I] = sort(DR);
k=3;
estr = mode(labRTr(I(1:k,:)));
KNNR =  sum(estr==labRTe)/size(labRTe,2)
k=5;
estr = mode(labRTr(I(1:k,:)));
KNNR =  sum(estr==labRTe)/size(labRTe,2)
toc
