
load gesturenamesC
train = readtable('datasetNew/train_gestures.txt','ReadVariableNames',false,'Format','%s');
nParti=40;
dp = 0;
trainSet=table2array(train(:,1:end));

%creo un vettore di celle che contiene i gesti normalizzati di training traingest. C'è anche quello coi gesti non normalizzati. Più un vettore
%di etichette (labelTrain) associate
    
k=1;
for i=1:size(trainSet,1)
    gesti=dir(['datasetCagliaritano/' trainSet{i}]);
    
    for j=1:26
        gesto=sprintf('datasetCagliaritano/%s/%s',trainSet{i},gesturenames{j});
        
        tab = table2array(readtable(gesto));
        tab(:,4)=tab(:,4)-tab(1,4);
        if smoothing
            for(u=1:3)
                tab(:,u)=smooth( tab(:,u),'lowess');
            end
        end
        labelTrain(k)=j;
        traingesture{k} = tab;
        traingest{k}=normalizeCA(tab,nParti,resampling,centering,scaling); %vettore di celle con le matrici coi gesti di training
        k=k+1;
    end
    
    
end


% stesso per il test

test=readtable('datasetCagliaritano/test_gestures.txt','ReadVariableNames',false,'Format','%s');

testSet=table2array(test(:,1:end));

k=1;
for i=1:size(testSet,1)
    gesti=dir(['datasetCagliaritano/' testSet{i}]);
    
    for j=1:26
        gesto=sprintf('datasetCagliaritano/%s/%s',testSet{i},gesturenames{j});
        labelTest(k)=j;
        tab = table2array(readtable(gesto));
        tab(:,4)=tab(:,4)-tab(1,4); % sottraggo dal timestamp il tempo iniziale
        if smoothing
            for(u=1:3)
                tab(:,u)=smooth( tab(:,u),'lowess');
            end
        end
        testgesture{k} = tab;
        testgest{k}=normalizeCA(tab,nParti,resampling,centering,scaling);%vettore di celle con le matrici coi gesti di test
        k=k+1;
    end
    
end


% calcolo matrice distanze training/test
tic
for i=1:size(traingest,2)
    for j=1:size(testgest,2)
        D(i,j)=gestureDistance(traingest{i},testgest{j},comptype,rotation);
        
    end
end
iterTime = toc/(size(testgest,2));

[val ind] = min(D); % prendo gli indici della minima distanza training-test

est=labelTrain(ind); % queste sono le etichette stimate con 1-Nearest Neighbor

KNN1 =  sum(est==labelTest)/size(labelTest,2) % questa è la percentuale di attribuzioni corrette (97%)

K=3;

[Y,I] = sort(D);
est2=mode(labelTrain(I(1:K,:)));
KNN3 =  sum(est2==labelTest)/size(labelTest,2)

K=5;

[Y,I] = sort(D);
est2=mode(labelTrain(I(1:K,:)));
KNN5 =  sum(est2==labelTest)/size(labelTest,2)





