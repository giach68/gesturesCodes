%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     The evaluation code for 
%%% SHREC'11 -- Shape Retrieval Contest of Non-rigid 3D Watertight Meshes  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input:
%       classification file  -- "test.cla"
%       distance matrix      -- "result.matrix"

%Output:
%       distance matrix for PSB evaluation code -- "PSBresult.matrix"
%       evaluation file      -- "result.txt"

%Evaluation measures:
%       NN, 1-Tier, 2-Tier, e_Measure, DCG

%Author: Zhouhui Lian
%Email: lianzhouhui@yahoo.com.cn
%Date: February 09, 2011
%@NIST, Gaithersburg, US

%Please cite:
% SHREC'11 Track: Shape Retrieval Contest of Non-rigid 3D Watertight Meshes, 3DOR'11, 2011

%clear;
%The folder that contains the distance matrix and classification file,
%Please change the name of the folder if necessary!!
%Please change the name of the distance matrix file to "result.matrix"!!
%Please change the name of the classification file to "test.cla"!!
filePath = './';
%filePath = '/home/giach/Dropbox/Devel/testAPT/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization
avgFirst_NN = 0;
first_NN = 0;
avgFirst_Tier = 0;
first_Tier = 0;
avgSecond_Tier = 0;
second_Tier = 0;
avgE_Measure = 0;
e_Measure = 0;
idealDCG = 0;
DCG = 0;
avgDCG = 0;
K1 = 0;
K2 = 0;

testCategoryList.categories(1).name(1) = 0;
testCategoryList.categories(1).numModels = 0;

testCategoryList.numCategories = 0;
testCategoryList.numTotalModels = 0;
testCategoryList.modelsNo(1) = 0;
testCategoryList.classNo(1) = 0;

%%%%%%Read the classification file
disFileName = 'distance.matrix'

claFileName = 'datasetCagliaritano/test.cla'

%laFileName = 'DATA_TEX1_FILES/testtest.cla';
%claFileName = 'funzioni_per_test/testReal240.cla'
%claFileName = sprintf('%stestReal200.cla',filePath);
%claFileName = '../MPI-FAUST/test/test.cla'
fp = fopen(claFileName,'r');
%fp = fopen('test.cla','r');
%Check file header
strTemp = fscanf(fp,'%s',1);
if ~strcmp(strTemp,'PSB')
    display('The format of your classification file is incorrect!');
    return;
end
strTemp = fscanf(fp,'%s',1);
if ~strcmp(strTemp,'1')
    display('The format of your classification file is incorrect!');
    return;
end

numCategories = fscanf(fp,'%d',1)
numTotalModels = fscanf(fp,'%d',1)

testCategoryList.numCategories = numCategories;
testCategoryList.numTotalModels = numTotalModels;

currNumCategories = 0;
currNumTotalModels = 0;

for i=1:numCategories
    currNumCategories = i
    testCategoryList.categories(currNumCategories).name = fscanf(fp,'%s',1);
    fscanf(fp,'%d',1);
    numModels = fscanf(fp,'%d',1);
    testCategoryList.categories(currNumCategories).numModels = numModels;
    for j=1:numModels
        currNumTotalModels = currNumTotalModels+1;
      %  testCategoryList.modelsNo(currNumTotalModels) = fscanf(fp,'%d',1)+1;
        testCategoryList.modelsNo(currNumTotalModels) = fscanf(fp,'%d',1);
        testCategoryList.classNo(currNumTotalModels) = currNumCategories;
    end
end

if (currNumTotalModels~=numTotalModels)
    display('The format of your classification file is incorrect!');
    return;
else
    display('The format of your classification file is correct!');
end
fclose(fp);

%%%%%%Read the distance matrix
%disFileName = sprintf('%sresult.matrix',filePath);
fp = fopen(disFileName,'r');
matrixInput = fscanf(fp,'%f');
numElement = size(matrixInput,1);
if (numElement~=(numTotalModels*numTotalModels))
    display('The format of your distance file is incorrect!');
    return;
else
    display('The format of your distance file is correct!');    
end
fclose(fp);


%%%%%%%Output the new distance matrix that can be used by the PSB evaluation code
disNewFileName = sprintf('%sPSBresult.matrix',filePath);
fp = fopen(disNewFileName,'w');
matrixDis(numTotalModels,numTotalModels) = 0;
for i=1:numTotalModels
    for j=1:numTotalModels
        iNew = testCategoryList.modelsNo(i);
        jNew = testCategoryList.modelsNo(j);
        matrixDis(i,j) = matrixInput((iNew-1)*numTotalModels+jNew);
        fprintf(fp,'%.6f ',matrixDis(i,j));             
    end
    fprintf(fp,'\n');
end

fclose(fp);

%%%%%%%Evaluation
matrixNo(1:numTotalModels) = 0;
modelNo(1:numTotalModels) = 0;
tempDis(1:numTotalModels) = 0;

 
 for i = 1:numTotalModels
     precision(i) = 0;
     recall(i) = 0; 
     precS(i,1:testCategoryList.numCategories ) = 0;
     recaS(i,1:testCategoryList.numCategories ) = 0;
 end
 
for i = 1:numTotalModels
    matrixDis(i,i) = -Inf;
    [tempDis, modelNo] = sort(matrixDis(i,:));
    for k = 1:numTotalModels
        matrixNo(k) = testCategoryList.classNo(modelNo(k));
    end
	
	count = 0;
  
	K1 = testCategoryList.categories(matrixNo(1)).numModels-1;
	K2 = 2*K1;
	DCG = 0;
	idealDCG = 1;
	for j = 2:K1
		idealDCG = idealDCG + log(2.0)/log(j);
    end
   
	for j = 1:numTotalModels	
        
		if (matrixNo(j) == testCategoryList.classNo(i))
			count = count+1;
			if (j ~= 1)
				if (j == 2)
					first_NN = first_NN+1;
					DCG = 1;
				else
					DCG = DCG + log(2.0)/log(j-1);
				end
            end
            
		end
		if (j == K1+1)
			first_Tier = (count-1)*1.0/K1;
			avgFirst_Tier = avgFirst_Tier + first_Tier;
		end
		if (j == K2+1)
			second_Tier = (count-1)*1.0/K1;
			avgSecond_Tier = avgSecond_Tier + second_Tier;
        end
      
            recall(j) = (recall(j) + (count-1)*1.0/(testCategoryList.categories(testCategoryList.classNo(j)).numModels-1));
            precision(j) = precision(j) + (count)*1.0/(j);
            
            recaS(j,testCategoryList.classNo(i)) = recaS(j,testCategoryList.classNo(i)) + (count-1)*1.0/(testCategoryList.categories(testCategoryList.classNo(j)).numModels-1);
            precS(j,testCategoryList.classNo(i)) = precS(j,testCategoryList.classNo(i)) + (count)*1.0/(j);
        
             
            
		if (j == 33)
			e_Measure = (count-1)*2.0/(K1+32);
			avgE_Measure = avgE_Measure + e_Measure;
		end
	end
	DCG = DCG/idealDCG;
	avgDCG = avgDCG + DCG;
	
end
precision = precision/testCategoryList.numTotalModels ;
recall = recall/testCategoryList.numTotalModels ;

for i=1:numTotalModels
    for j=1:numCategories
precS(i,j) = precS(i,j)/testCategoryList.categories(j).numModels;
recaS(i,j) = recaS(i,j)/testCategoryList.categories(j).numModels;
    end
end
 for j=1:numCategories
precS(i+1,j) = 0;
recaS(i+1,j) = 1;
clear pcon rcon
[pcon,rcon] = consolidator(recaS(:,j),precS(:,j));
prec2=spline(pcon,rcon,0:0.05:1);
avPre(j)=sum(prec2)/numel(prec2);
    end
mAp=mean(avPre)

avgFirst_Tier = avgFirst_Tier/numTotalModels;
avgSecond_Tier = avgSecond_Tier/numTotalModels;
avgE_Measure = avgE_Measure/numTotalModels;
avgDCG = avgDCG/numTotalModels;
avgFirst_NN = first_NN/numTotalModels;


evalFileName = sprintf('%sresult.txt',filePath);
fp = fopen(evalFileName,'w');

fprintf(fp,'NN:\n ');
fprintf(fp,'%.4f\n ',avgFirst_NN);	
fprintf(fp,'1_Tier:\n ');
fprintf(fp,'%.4f\n ',avgFirst_Tier);
fprintf(fp,'2_Tier:\n ');
fprintf(fp,'%.4f\n ',avgSecond_Tier);
fprintf(fp,'e_Measure:\n ');
fprintf(fp,'%.4f\n ',avgE_Measure);
fprintf(fp,'DCG:\n ');
fprintf(fp,'%.4f\n ',avgDCG);

disp(' ')
disp('------------------------ RESULTS -------------------------');
strTemp = sprintf('NN      1-Tier      2-Tier     e-Measure     DCG\n');
disp(strTemp);
strTemp = sprintf('%.4f  %.4f      %.4f     %.4f        %.4f',avgFirst_NN,avgFirst_Tier,avgSecond_Tier,avgE_Measure,avgDCG);
%strTemp = sprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f',avgFirst_NN,avgFirst_Tier,avgSecond_Tier,avgE_Measure,avgDCG);
disp(strTemp);
disp('----------------------------------------------------------');
disp(' ')

fclose(fp);
% 
% count = 0;
% [recall2 ii] = unique(recall);
% precision2=precision(ii)
% 
% yy=spline(recall2,precision2,0:0.05:1)'
% 
% for k=1:15
% % recall=recaS([1:150 160:20:720],k)';
% % precision=precS([1:150 160:20:720],k)';
% [recall ii] = unique(recaS(:,k))
% precision=precS(ii,k)'
% yy((k-1)*21+1:k*21)=spline(recall,precision,0:0.05:1)'
% end
