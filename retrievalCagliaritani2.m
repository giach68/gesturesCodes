
load gesturenames
train = readtable('datasetNew/all_gestures.txt','ReadVariableNames',false,'Format','%s');

trainSet=table2array(train(:,1:end));

%creo un vettore di celle che contiene i gesti normalizzati di training traingest. C'è anche quello coi gesti non normalizzati. Più un vettore
%di etichette (labelTrain) associate


cond=1;

nParti=40;

for smoothing = 1:2
    
for resampling=1:2; % 1-a passi di arc length regolari o 2-di tempo regolari

centering=1; % 1 centroide -   2 origine

scaling=1; % 1 su lunghezza 2 bounding box(protractor 3D), 3 max distanza origine

comptype=1; % main trajectory comparison:  1 squared Euclid 2-absolute(cityblock)  3-DTW  4- min point to point
rotation=0;


k=1;
for i=1:size(trainSet,1)
      %  trainFeatures{i,1}=trainSet(i,1);
     %   target=trainFeatures{i,1}; ??
             
       gesti=dir(['datasetNew/' trainSet{i}]);
       
        for j=1:26
            gesto=sprintf('datasetNew/%s/%s',trainSet{i},gesturenames{j});
           
            tab = table2array(readtable(gesto));
   if(smoothing==2)
        for(u=1:3)
        tab(:,u)=smooth( tab(:,u),'lowess');
    end
   end       
            
            tab(:,4) = tab(:,4)-tab(1,4);
            labelTrain(k)=j;
            traingesture{k} = tab;
            traingest{k}=normalizeP(tab,nParti,resampling,centering,scaling);


            k=k+1;          
        end
  
end



for i=1:size(traingest,2)
    for j=1:size(traingest,2)
        D(i,j)=gestureDistance(traingest{i},traingest{j},comptype,rotation);
            
    end
end


fp = fopen('distance.matrix','w');

for i = 1:size(traingest,2)
    for j=1:size(traingest,2)
      
        fprintf(fp,'%f ',D(i,j));
    end
    fprintf(fp,'\n');
end

EvaluationSC
NN{cond}=avgFirst_NN;
FT{cond}=avgFirst_Tier;
ST{cond}=avgSecond_Tier;
E{cond}=avgE_Measure;
Dcg{cond}=avgDCG;
MAP{cond}=mAp;
prec{cond}=precision;
reca{cond}=recall;
cond=cond+1;
end
end
