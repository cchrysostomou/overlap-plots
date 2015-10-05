%if type = 1, use the ranking method
%if type = 2, use the pdf method
%if type = 3, use the cdf method

%use log data will convert the pertinent numbers into logs

%cell data is the complete cell array which contains transcript levels of
%each cluster in each tissue (i.e. MatrixData97_Mouse5)

%colUse defines which columns contain the numerical data
function results = venRankPDFCDFDiagram3(cellData,colUse,contains_header, type,useLogData,names,varargin)

%list of indexes that store the colored values (R,G,B) colors
colVal = [1,64];

%limits of the possible transparency values
trpVal = [0.15,1];

%how far to separate cirlces
enlargeProp = 1.5;

%this is the linear transformation to change the ranked values of unique
%sequences.  the first row corresponds tot he best rank whil ethe 2nd
%column refers to the worst length.
rankStretch = [30,0];

%load the variable storing the possible colors
load colValVar;
col = colors;
minR = 0;

%these will change any settings if we want to on the go (the transparaence
%and color settings)
if nargin>5
    findStr = cellfun(@ischar,varargin);
    tempVargin = varargin(findStr);
    
    c = find(~cellfun('isempty',regexp(tempVargin,'colThr','once')));
    g = find(~cellfun('isempty',regexp(tempVargin,'trpThr','once')));
    h = find(~cellfun('isempty',regexp(tempVargin,'enlarge','once')));
    z = find(~cellfun('isempty',regexp(tempVargin,'colors', 'once')));
    y = find(~cellfun('isempty',regexp(tempVargin,'ColCutoff','once')));
    
    if ~isempty(c)
        newPos = c(1)*2;
        colVal = varargin{newPos};
    end
    
    if ~isempty(y)
        newPos = y(1)*2;
        minR = varargin{newPos};
    end
    
    if ~isempty(g)
        newPos = g(1)*2;
        trpVal = varargin{newPos};
    end
    
    if ~isempty(h)
        newPos = h(1)*2;
        enlargeProp = varargin{newPos};
    end
    
    if ~isempty(z)
        newPos = z(1)*2;
        col = varargin{newPos};
    end
    
end

if iscell(cellData)
    if contains_header > 0
        subData = cell2mat(cellData(2:end,colUse));
    else
        subData = cell2mat(cellData(:,colUse));
    end
else    
    subData = cellData(:,colUse);
end

totalCounts = sum(subData,1);

legendmax = [0,1];
legendmin = [1,1];


for i = 1:size(subData,2)
    
    tempR1 = [];
    temp = [];
    
    
    %temp stores a list of the data in its sorted positions
    temp = sort(subData(:,i),'descend');
    
    
    temp2 = temp(:,1);
    
    %the second column of temp2 corresponds to frequencies of the total
    %counts
    temp2(:,2) = temp(:,1)/sum(temp(:,1));
    
    
    %temp R1 selects only the unique rows from temp 2 (ignores repetitve
    %sequences)
    
    %rankPos will return the relative rank of these sequence counts
    [tempR1,rankPos] = unique(temp2,'rows');
    tempR1=flipud(tempR1);
    
    %this function will organize the sequences into tied rankings
    %where the largest count has a rank of 1
    tiedR = tiedrank(-temp2(:,1));
    
    

    if useLogData ==1               
        %tempR1(:,3) = (tempR1(:,1)+1)/(totalCounts(i)+legendCount(i,1));        
        tempR1(:,3) = log(tempR1(:,2));        
        
    elseif useLogData == 2        
        tempR1(:,3) = (tempR1(:,2)).^(1/3);
        
    else        
        tempR1(:,3) = tempR1(:,2);
        
    end
    
    
    
    %tempR1(:,6) = flipud(rankPos);
    tempR1(:,6) = unique(tiedR);
    
    
    %this is the minimum possible value for minTr
    minTr = trpVal(1);
    maxTr = trpVal(2);
    
    
    minVC = colVal(1);
    maxVC = colVal(2);
    
    [minVal(i),maxVal(i)] = normalizeColor(tempR1(:,3));
    
    %minRank is just the minimum of the column tempR1.  IN this case
    %they are not ranking data.
    
    
    
    minRank = minVal(i);%tempR1(end-minR,3);
    maxRank = maxVal(i);%;tempR1(1,3);
    
    
    minRuse(i) = minRank;
    maxRuse(i) = maxRank;
    hueTrans_M(i,1) = (maxTr-minTr)/(maxRank-minRank);
    hueTrans_M(i,2) = minTr-hueTrans_M(i,1)*minRank;
    
    hueColor_M(i,1) = (maxVC-minVC)/(maxRank-minRank);
    hueColor_M(i,2) = minVC-hueColor_M(i,1)*minRank;
    
    
    %lets use this system if we just use the indexes of the matirx,
    %that is we just care about the numberr of possible unique
    %frequency values
    %   rankHeightMax = 1;
    %   rankHeightMin = length(tempR1)-1;
    
    %use these two values if you want to use their true ranks as a
    %plot
    %rankHeightMax = tempR1(1,6);
    %rankHeightMin = tempR1(end-1,6);
    
    %we will use these rankings if we just want to only plot sequences
    %in the z axis that are in to top X.  That is, any sequence above a
    %top X rank will be plotted in the Z axis while all other rankings
    %will be a height of 0
    posToFind = find(tempR1(:,6)<=rankStretch(1));    
    if isempty(posToFind)
        posToFind = 1;
    elseif posToFind(end)+1 > size(tempR1,1)
        posToFind = find(tempR1(:,6)<=max(tempR1(:,6)-1));       
    end
    
    
    rankHeightMax = tempR1(1,6);
    rankHeightMin = tempR1(posToFind(end)+1,6);
    
    rankHeight_M(i,1) = (rankStretch(1)-rankStretch(2))/(rankHeightMax-rankHeightMin);
    rankHeight_M(i,2) = rankStretch(2)-rankHeight_M(i,1)*rankHeightMin;
    
    
    
    tempR1(:,6) = tempR1(:,6)*rankHeight_M(i,1)+rankHeight_M(i,2);
    tempR1((tempR1(:,6)<0),6)=0;
    
    tempR1(:,4) = tempR1(:,3)*hueTrans_M(i,1)+hueTrans_M(i,2);
    tempR1(:,5) = tempR1(:,3)*hueColor_M(i,1)+hueColor_M(i,2);     
    tempR1(:,5) = round(tempR1(:,5));
    
    %because we linearize the values, some will be less than 0 and others
    %treater than 1 becuase of the cutoffs., so any transparency less than
    %0 we make 0 and any greater than 1 we make 1
    tempR1(tempR1(:,4)<minTr,4) = minTr;
    tempR1(tempR1(:,4)>maxTr,4) = maxTr;
    
    %same argument for "color values"
    tempR1(tempR1(:,5)<minVC,5) = minVC;
    tempR1(tempR1(:,5)>maxVC,5) = maxVC;
    
    
    
    rankData{i} = tempR1;
    [a1,b1] = min(tempR1(tempR1(:,2)>0,2));
    
    if legendmin(1)>a1
        legendmin(1) = a1;
        legendmin(2) = tempR1(b1,3);
    end
    
    [a2,b2] = max(tempR1(:,2));
    if legendmax(1)<a2
        legendmax(1) = a2;
        legendmax(2) = tempR1(b2,3);
    end
        
    
        
end



numPoints = 50;
labelPoints = 5;
stepF = (legendmax(2)-legendmin(2))/numPoints;

stepLabel = (legendmax(2)-legendmin(2))/labelPoints;

stepLegend{1} = [legendmin(2):stepF:legendmax(2)]';
stepLegend{2} = stepLegend{1};

stepLabel = num2cell([legendmin(2):stepLabel:legendmax(2)]');


for j = 1:size(stepLabel,1)
    if useLogData ==1
        stepLabel{j,1} = exp(stepLabel{j,1});
    elseif useLogData == 2
        stepLabel{j,2} = stepLabel{j,1}^3;
    else
        stepLabel{j,2} = stepLabel{j,1};
        
    end
end


    
for j = 1:size(stepLegend{1},1)
    
    
    %if we are using the LogData parameter, then we we set the 3rd
    %column of tempR1 to equal the log of unique sequence counts.  If
    %not then we set it to the sqrt of unique sequences.
    
    for i = 1:size(subData,2)
        stepLegend{1}(j,i+1) = stepLegend{1}(j,1)*hueTrans_M(i,1)+hueTrans_M(i,2);
        if stepLegend{1}(j,i+1)>maxTr
            stepLegend{1}(j,i+1)=maxTr;
        elseif stepLegend{1}(j,i+1)<minTr
            stepLegend{1}(j,i+1)=minTr;
        end
        
        
        
        stepLegend{2}(j,i+1) = round(stepLegend{2}(j,1)*hueColor_M(i,1)+hueColor_M(i,2));

        if stepLegend{2}(j,i+1)>maxVC
            stepLegend{2}(j,i+1)=maxVC;
        elseif stepLegend{2}(j,i+1)<minVC
            stepLegend{2}(j,i+1)=minVC;
        end

        
        
    end
    
    
end

organizedData = OrganizeOverlappingData(subData,enlargeProp);
f = figure;
for i = 1:length(organizedData)
      makeCube(organizedData{i},rankData(:,organizedData{i}{4}),type,col(organizedData{i}{4}),colVal,trpVal,totalCounts);
end
set(gca,'DataAspectRatio',[1,1,1])

drawLegend(stepLegend{1},stepLegend{2},col,names,stepLabel);


for i = 1:7
 results{1,i} = length(organizedData{i}{1});
 for j = 1:length(organizedData{i}{4})
  pos = organizedData{i}{4}(j);
 results{2,i}(pos) =sum(organizedData{i}{1}(:,j))/totalCounts(pos);
 totalFreq(i,pos) = results{2,i}(pos);
 end

 
 end

    results{2,8} = sum(totalFreq);





%vec contains our normalized data/either the raw data sequences, or the
%"sqrt" frequency, or the "log frequency"

%for the reaminder of the logarithm
   %trim the edges of vec so that we only keep 70% of the sequences
   %fit these sequences to a line,
   %based on this linear fit what are the extremes of the line to use as
   %color thresholds for the data
function [minV,maxV] = normalizeColor(vec)
 cutofC = 0.7;
    mid = (1-cutofC)/2;
    p1 = ceil(mid*size(vec,1))+1;
    p2 = floor((1-mid)*size(vec,1));
    vec(:,2) = 1:length(vec);
    
    dataVal = vec(p1:p2,[1,2]);
    %dataVal(:,2) = 1:length(dataVal);
    [x,y,z] = regression(dataVal(:,2),dataVal(:,1),'one');
    sim(:,1) = vec(:,2)*y+z;
   
    maxV = sim(1,1);
    minV = sim(end,1);

    
    
function drawLegend(thrData,colData,colors_Circle,names,stepLabel)

   

numSamples = size(thrData,2)-1;
legendWid = 0.3;

circleAxes = gca;
pos = get(circleAxes,'position');
set(circleAxes,'Position',[0.05,pos(2),pos(3),pos(4)]);
set(circleAxes,'units','pixels');
pos = get(circleAxes,'position');

newPos_x = pos(1)+pos(3)+80;
newPos_y = pos(2);
newPos_yLen = pos(4);

newPos_xLen = pos(3)*legendWid;

figPos = get(gcf,'position');

if (figPos(3)<=(newPos_x+newPos_xLen+90))
    diffSize = newPos_x+newPos_xLen - figPos(3);
    figPos(3) = diffSize+figPos(3)+90;
end

set(gcf,'position',figPos);

legendA = axes;

set(legendA,'units','pixels');

set(legendA,'position',[newPos_x,newPos_y,newPos_xLen,newPos_yLen]);

set(legendA,'units','normalized');
set(circleAxes,'units','normalized');

posPlot = get(legendA,'xlim');
posPloty= get(legendA,'ylim');

xPos = 0;
yPos = 0;

widSample = 1/numSamples;

ywid = 1;

index = thrData(:,1);
numFaces = size(index,1);
widFace = ywid/(numFaces-1);

vert = [];


for i = 1:numSamples
    
    colVal = colors_Circle{i};
    legendColData(:,1) = thrData(:,1);
    legendColData(:,3) = thrData(:,i+1);
    legendColData(:,2) = colData(:,i+1);
    
    
    startVector = [xPos,yPos];
    for j = 1:numFaces-1
        x(1) = startVector(1);
        y(1) = startVector(2);
        
        try
            c(1,:) = colVal(round(legendColData(j,2)),:);
        catch
            c(1,:) = [1,1,1];
        end
        
        t(1) = legendColData(j,3);
        
        x(2) = startVector(1)+widSample;
        y(2) = startVector(2);
        try
            c(2,:) = colVal(round(legendColData(j,2)),:);
        catch
            c(2,:) = [1,1,1];
        end
        t(2) = legendColData(j,3);
        
        
        x(3) = startVector(1)+widSample;
        y(3) = startVector(2)+widFace;
        
        try
            c(3,:) = colVal(round(legendColData(j+1,2)),:);
        catch
            c(3,:) = [1,1,1];
        end
        
        t(3) = legendColData(j+1,3);
        
        x(4) = startVector(1);
        y(4) = startVector(2)+widFace;
        try
            c(4,:) = colVal(round(legendColData(j+1,2)),:);
        catch
            c(4,:) = [1,1,1];
        end
        t(4) = legendColData(j+1,3);
        
        startVector = [x(4),y(4)];
        
        patch('xdata',x,'ydata',y,'facevertexcdata',c,'facecolor','interp','linestyle','none','cdatamapping','direct','facevertexalphadata',t','facealpha','interp','alphadatamapping','none');
        
    end
    xPos = xPos+widSample;
end

set(legendA,'ylim',[0,1]);

newTickPos = 0:1/(length(stepLabel)-1):1;
set(legendA,'ytick',newTickPos);
stepLabel{1} = 0;%%changes lowest freq to 0 for ease
set(legendA,'yticklabel',stepLabel(:,1));

set(gca,'xtick',[0.33/2,0.33/2+0.33,0.33/2+2*0.33]);
set(gca,'xticklabel',names);
%set(legendA,'yticklabel',yticklab);

function makeCube(seqInfo,rankInfo, plotStyle,colors,colVal,trpVal,totalCounts)

circleArea = seqInfo{2};
rad = sqrt(seqInfo{2}/pi);

numSamples = size(seqInfo{1},2);


myCounts = totalCounts(:,seqInfo{4});

centerCircle = seqInfo{3};

angle = 0:pi/(100*pi):2*pi;

if circleArea == 0
    set(gca,'nextplot','add');
    scatter(centerCircle(1),centerCircle(2),10,'k');
else
    %set(gca,'nextplot','add');
    %scatter(centerCircle(1),centerCircle(2),10,'k');
   
    s2 = zeros(size(seqInfo{1},1),4);%1:size(seqInfo{1},1);
    s2(:,1) = 1:size(s2,1);
    for i = 1:numSamples
        s2(:,i+1) = seqInfo{1}(:,i)/myCounts(i);
    end
    
    
    %%comment out if dont want to use this method
%     for i = 1:size(s2,1)
%         tempAll(i,1) = sqrt(sum(s2(i,2:end).^2));
%     end
%     
%     tempAll(:,2) = 1:length(tempAll);
%         
%    s2 = flipud(sortrows(s2,[2,3,4]));
%    %% 
   
   
    tempAll(:,2) = s2(:,1);
    
    for i = 1:numSamples
    
        setData = [];
        x2 = [];
        x1 = [];
        
        subRank = rankInfo{i};
        temp = seqInfo{1}(:,i);
        
        %    temp = sort(temp,'descend');
        %         try
        %             temp = temp(1:100);
        %         catch
        %             temp = temp;
        %         end
        %
        temp = seqInfo{1}(tempAll(:,2),i);
        %[setData(:,1),setData(:,2)] = unique(temp);
        %setData = flipud(setData);
        setData(:,1) = temp;
        setData(:,2) = 1:length(temp);
        
        rIn = 0;
        fPre = 1;
        
        centerCircle = seqInfo{3};
        
        
        
        for j = 1:size(setData,1)
            val = setData(j,2);
            valFreq = setData(j,1);
            
            %rOut = sqrt(val/pi);
            rOut = val;
            %rIn = sqrt(setData(j+1),2);
            
            if j == 1
                preVal = valFreq;
            else
                preVal = setData(j-1,1);
            end
            
            tempFind = find(subRank(:,1) == valFreq);
            
            if isempty(temp)
                temp1Find = find(subRank(:,1)>valFreq);
                tempFind = temp1Find(1);
            end
            
            colData = tempFind;
            
            if j == 1
                fPre = colData;
            end
            
            
            x2(:,1) = rOut*seqInfo{5}{i}(:,1)+centerCircle(1);
            x2(:,2) = rOut*seqInfo{5}{i}(:,2)+centerCircle(2);
            x2(:,3) = 0;%subRank(colData,6);%valFreq;%seqInfo{5}{i}(:,2);
            
            x1(:,1) = rIn*seqInfo{5}{i}(:,1)+centerCircle(1);
            x1(:,2) = rIn*seqInfo{5}{i}(:,2)+centerCircle(2);
            x1(:,3) = 0;%subRank(fPre,6);%preVal;
            
            vComb = [x1;x2];
            
            
            
            colBeMap = subRank(fPre,5);
            colMap = subRank(colData,5);
            colTrBe = subRank(fPre,4);
            colTr = subRank(colData,4);
            
            
            sumFreq = sum(setData(:,1));
            
            
            
            plotColorMap(colors{i},colMap,colBeMap,colTrBe, colTr, vComb,x1,x2);
            
            
            fPre = colData;
            rIn = rOut;
            
        end
        
%              ang = 0:0.01:2*pi;
%              xp = circleArea*cos(ang)+centerCircle(1);
%              yp = circleArea*sin(ang)+centerCircle(2);
%              set(gca,'nextplot','add');
%              plot(xp,yp);
%             
            
            
            
            
            
         % plotColorMap(colors{i},colMap,colBeMap,colTrBe, colTr, vComb,x1,x2);
    end
end

function plotColorMap(colors,col,colBe,colTrBe,colTr,vert,x1,x2)
if(size(colors,1)==64)
    colormapIndex = colors;
else
    colormapIndex = colormap;
end

colBeColor = colormapIndex(colBe,:);
colColor = colormapIndex(col,:);

c1 = repmat(colBeColor,[size(x1,1),1]);
c2 = repmat(colColor,[size(x2,1),1]);

cT1 = repmat(colTrBe,[size(x1,1),1]);
cT2 = repmat(colTr,[size(x1,1),1]);

c = [c1;c2];

cT = [cT1;cT2];

%v{i}=vComb;
%c{i} = [c1;c2];

u(:,1) = 1:size(x1,1);
u(:,2) = size(x1,1)+1:size(x1,1)+size(x2,1);
for x = 1:length(u)-1
    f1(x,1) = u(x,1);
    f1(x,2) = u(x+1,1);
    f1(x,3) = u(x+1,2);
    f1(x,4) = u(x,2);
end

patch('vertices',vert,'faces',f1,'facevertexcdata',c,'facecolor','flat','linestyle','none','cdatamapping','direct','facevertexalphadata',cT,'facealpha','interp','alphadatamapping','none');


%%This will take the matrix of data and organize it the data into each
%%specific sector,
%data{1} refers to seuqences only foudn in set1
%data{2} refers to sequences in only 2
%data{3} = refers to sequences ionly in 3rd dataset
%data{4} = refers to shared sequences between 1 and 2
%data{5} referes to shared between 1 adn 2
%data{6} refers to shared between 2 and 3
%data{7} refers to seuqences shared in all

%Each data variable has a 3 sub sets o finformation
% first {1} = stores the frequency and rank information for that
% subset of data
% second {2} = refers to the number of unique seuqneces that make
% up that data set
% third {3} = refers to the x and y coordinates of the circle that
% plots the population
function data = OrganizeOverlappingData(subData,enlargeProp)

if nargin==2
    cutoff(1) = 2;
    cutoff(2) = 1;
    cutoff(3) = 30;
else
    cutoff = varargin{1};
end

centerx = 0;
centery = 0;


%%determine how to filter the sequences...that is, it determines at what
%%cutoff to plot sequences
c1 = [];
for i = 1:size(subData,2)
    subSub(:,1) = 1:size(subData,1);
    subSub(:,2) = subData(:,i);
    subSub(:,3) = subSub(:,2)/sum(subSub(:,2));
    
    subData2(:,i) = subSub(:,3);%store the frequency values in subData2
    
    subSub = flipud(sortrows(subSub,3));
    sumV = 0;
    for j = 1:size(subData,1)
        subSub(j,4) = subSub(j,3)+sumV;
        sumV = subSub(j,4);
    end
    
    subSub(:,5) = tiedrank(-subSub(:,2));
    
    %sort sequences in this tissue by rank
    subSub = sortrows(subSub,5);
    
    if (cutoff(1)==1)
        g = find(subSub(:,4)<=cutoff(2)); %chooses sequences below this CDF cutoff for that tissue
    elseif (cutoff(1)==2)
        g = find(subSub(:,5)<=cutoff(3)); %chooses sequences with a ranking below this cutoff (cutoff
    end
    
           
    %%we use this code to create sequences in an order in which we first
    %%sort by highest frequencies in BMPC firs,t then we choose the next
    %%top ranked sequences in LNPC, then we choose the next top ranked
    %%squences in SPPC    
        for z = 1:length(g)
        f = 0;
        for zz = 1:length(c1)
            if subSub(g(z),1) == c1(zz)
                f = 1;
            end
        end
        if f == 0
            c1 = [c1;subSub(g(z),1)];
        end
        end
    
    
end

%c1 = unique(c1);

temp = subData;

subData = subData(c1,:);




totalNumUnique = size(subData,1);
maxHorizDist = enlargeProp*50;%sqrt(totalNumUnique);

%extract sequences found only in dataset1
sub1(:,1) = subData(:,1)>0;
sub1(:,2) = subData(:,2)==0;
sub1(:,3) = subData(:,3)==0;

sumD = sum(sub1,2);


justD1 = find(sumD==3);
data{1}{1} = subData(justD1,1);
data{1}{2} = length(justD1);
data{1}{3} = [centerx-maxHorizDist,centery];
data{1}{4} = 1;
angle = 0:pi/(100*pi):2*pi;

data{1}{5}{1}(:,1) = cos(angle);
data{1}{5}{1}(:,2) = sin(angle);


%extract sequences found only in dataset2
sub1(:,1) = subData(:,1)==0;
sub1(:,2) = subData(:,2)>0;
sub1(:,3) = subData(:,3)==0;
sumD = sum(sub1,2);
angle = 0:pi/(100*pi):2*pi;
justD1 = find(sumD==3);
data{2}{1} = subData(justD1,2);
data{2}{2} = length(justD1);
data{2}{3} = [centerx+maxHorizDist,centery];
data{2}{4} = 2;

data{2}{5}{1}(:,1) = cos(angle);
data{2}{5}{1}(:,2) = sin(angle);



%extract sequences found only in dataset3
sub1(:,1) = subData(:,1)==0;
sub1(:,2) = subData(:,2)==0;
sub1(:,3) = subData(:,3)>0;
sumD = sum(sub1,2);

justD1 = find(sumD==3);
data{3}{1} = subData(justD1,3);
data{3}{2} = length(justD1);
data{3}{3} = [centerx,centery-maxHorizDist*2*cos(30/180*pi)];
data{3}{4} = 3;

angle = 0:pi/(100*pi):2*pi;
data{3}{5}{1}(:,1) = cos(angle);
data{3}{5}{1}(:,2) = sin(angle);


%extract sequences found only in dataset 1 and 2
sub1(:,1) = subData(:,1)>0;
sub1(:,2) = subData(:,2)>0;
sub1(:,3) = subData(:,3)==0;
sumD = sum(sub1,2);

justD1 = find(sumD==3);
data{4}{1} = subData(justD1,1:2);
data{4}{2} = length(justD1);
data{4}{3} = [centerx,centery];
data{4}{4} = [1,2];


subPi = pi/2;
stopPi = 3*pi/2;
step = (stopPi-subPi)/200;
angle = subPi:step:stopPi;


data{4}{5}{1}(:,1) = (cos(angle));
data{4}{5}{1}(:,2) = (sin(angle));


subPi = pi/2;
stopPi = -pi/2;
step = (stopPi-subPi)/200;
angle = subPi:step:stopPi;

data{4}{5}{2}(:,1) =(cos(angle));
data{4}{5}{2}(:,2) = (sin(angle));




%extract sequences found only in dataset 1 and 3
sub1(:,1) = subData(:,1)>0;
sub1(:,2) = subData(:,2)==0;
sub1(:,3) = subData(:,3)>0;
sumD = sum(sub1,2);

justD1 = find(sumD==3);
data{5}{1} = subData(justD1,[1,3]);
data{5}{2} = length(justD1);

data{5}{3} = [centerx-maxHorizDist+maxHorizDist*sin(30/180*pi),centery-maxHorizDist*cos(30/180*pi)];
data{5}{4} = [1,3];

subPi = 30/180*pi;
stopPi = (30+180)/180*pi;
step = (stopPi-subPi)/200;
angle = subPi:step:stopPi;

data{5}{5}{1}(:,1) =(cos(angle));
data{5}{5}{1}(:,2) =(sin(angle));

subPi = (30-180)/180*pi;
stopPi = 30/180*pi;
step = (stopPi-subPi)/200;
angle = subPi:step:stopPi;

data{5}{5}{2}(:,1) = (cos(angle));
data{5}{5}{2}(:,2) = (sin(angle));



%extract sequences found only in dataset 2 and 3
sub1(:,1) = subData(:,1)==0;
sub1(:,2) = subData(:,2)>0;
sub1(:,3) = subData(:,3)>0;
sumD = sum(sub1,2);

justD1 = find(sumD==3);
data{6}{1} = subData(justD1,[2,3]);
data{6}{2} = length(justD1);
data{6}{3} = [ centerx+maxHorizDist-maxHorizDist*sin(30/180*pi),centery-maxHorizDist*cos(30/180*pi)];
data{6}{4} = [2,3];



subPi = (60+90)/180*pi;
stopPi = -30/180*pi;
step = (stopPi-subPi)/200;
angle = subPi:step:stopPi;


data{6}{5}{1}(:,1) =(cos(angle));
data{6}{5}{1}(:,2) = (sin(angle));



subPi = (-30/180)*pi;
stopPi =(-30-180)/180*pi;
step = (stopPi-subPi)/200;
angle = subPi:step:stopPi;

data{6}{5}{2}(:,1) = (cos(angle));
data{6}{5}{2}(:,2) = (sin(angle));




%extract sequences found only in datasets 1, 2, and 3
sub1(:,1) = subData(:,1)>0;
sub1(:,2) = subData(:,2)>0;
sub1(:,3) = subData(:,3)>0;
sumD = sum(sub1,2);

justD1 = find(sumD==3);
data{7}{1} = subData(justD1,:);
data{7}{2} = length(justD1);
data{7}{3} = [centerx,centery-maxHorizDist*tan(30/180*pi)];
data{7}{4} = [1,2,3];

subPi = 90/180*pi;
stopPi = (90+120)/180*pi;
step = (stopPi-subPi)/200;

angle = subPi:step:stopPi;
data{7}{5}{1}(:,1) =(cos(angle));
data{7}{5}{1}(:,2) = (sin(angle));


subPi = (90-120)/180*pi;
stopPi = (90)/180*pi;
step = (stopPi-subPi)/200;

angle = subPi:step:stopPi;

data{7}{5}{2}(:,1) = (cos(angle));
data{7}{5}{2}(:,2) = (sin(angle));



subPi = (90+120)/180*pi;
stopPi = (330)/180*pi;
step = (stopPi-subPi)/200;
angle = subPi:step:stopPi;

data{7}{5}{3}(:,1) = (cos(angle));
data{7}{5}{3}(:,2) = (sin(angle));

