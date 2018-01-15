
load gesturenames
train = readtable('datasetCagliaritano/all_gestures.txt','ReadVariableNames',false,'Format','%s');

trainSet=table2array(train(:,1:end));

%creo un vettore di celle che contiene i gesti. 

k=1;
for i=1:size(trainSet,1)
 
             
       gesti=dir(['datasetCagliaritano/' trainSet{i}]);
       
        for j=1:26
            gesto=sprintf('datasetCagliaritano/%s/%s',trainSet{i},gesturenames{j});
           
            tab = table2array(readtable(gesto));
            tab(:,1)=smooth( tab(:,1));
            tab(:,2)=smooth( tab(:,2));
            tab(:,3)=smooth( tab(:,3));
            
            tab(:,4) = tab(:,4)-tab(1,4); % time stamp relativo ai campioni delle traiettorie
            labelTrain(k)=j;
            traingesture{k} = tab;
            traingest{k}=normalize(tab,100); % gesti ricampionati di lunghezza uniforme, numero di punti uniforme

            k=k+1;          
        end
  
end



for i=1:size(traingest,2)
    for j=1:size(traingest,2)
        D(i,j)=gestureDist(traingest{i},traingest{j});
            %sss(k,:)=traingest{k}(:);
    end
end


fp = fopen('distance.matrix','w');

for i = 1:size(traingest,2)
    for j=1:size(traingest,2)
      
        fprintf(fp,'%f ',D(i,j));
    end
    fprintf(fp,'\n');
end


