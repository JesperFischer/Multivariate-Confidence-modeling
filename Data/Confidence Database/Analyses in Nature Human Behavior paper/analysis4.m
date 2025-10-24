%--------------------------------------------------------------------------
% Analysis 4 from "The Confidence Database" paper by Rahnev, Desender, Lee, et al.
% 
% This analysis explores the files in the Confidence Database that come from perception
% and memory. Specifically, the code compares the type of confidence scale used in the 
% two fields.
%
% To run this analysis, place the file Database_Information.xlsx in the same folder as 
% this file.
%
% Written by Alan Lee. Last update: Sep 16, 2019.
%--------------------------------------------------------------------------

clear
close
clc

% Load spreadsheet
T = readtable('Database_Information.xlsx');
n_datasets = length(T.Name_in_database);
names_in_spreadsheet = sort(T.Name_in_database);

% create flags for each study
% Rule: each study is true (1) in only one of the following vectors
isPerception = strcmp(T.Category,'Perception');
isMemory = strcmp(T.Category,'Memory');
isCognitive = strcmp(T.Category,'Cognitive');
isMotor = strcmp(T.Category,'Motor');
isMixed = strncmp(T.Category,'Mixed',5);

% count number of studies for each confidence type
% confidence types: n-point (discrete) scale vs continuous scale
studyCatNames = {'Perception','Memory','Cognitive','Motor','Mixed'}';
studyCatNum = numel(studyCatNames);
pointScaleNum = NaN(studyCatNum,1);
continScaleNum = NaN(studyCatNum,1);
for cati = 1:studyCatNum
    studyCat = studyCatNames{cati};
    eval(sprintf('flagvec = is%s;',studyCat));
    pointScaleNum(cati) = nnz(endsWith(T.Confidence_scale(flagvec),'point'));
    continScaleNum(cati) = nnz(strncmp(T.Confidence_scale(flagvec),'continuous',10));
end
confTypeT = table(studyCatNames, pointScaleNum, continScaleNum)

% for perc vs memory studies,
%   count number of studies for each number of points on conf scale
pointBaseVec = 2:11;
pointBaseNum = numel(pointBaseVec);
studyCount = NaN(2,pointBaseNum);
for cati = 1:2
    if cati == 1    % perception studies
        percflag = logical(isPerception);
    else            % memory studies
        percflag = logical(isMemory);
    end
    for pbi = 1:pointBaseNum
        pointBase = pointBaseVec(pbi);
        eval(sprintf('confScaleName = ''%d-point'';',pointBase));
        studyCount(cati,pbi) = nnz(strcmp(T.Confidence_scale(percflag),confScaleName));
    end
end

% collapse scales with large number of points (7-11)
plotPointNvec = 2:7;
plotPointNNum = numel(plotPointNvec);
plotStudyCount = studyCount(:,1:(plotPointNNum-1));
plotStudyCount(:,plotPointNNum) = sum(studyCount(:,plotPointNNum:end),2);
plotStudyCountCell = cellstr(num2str(plotPointNvec'));
for ppni = 1:plotPointNNum
    plotStudyCountCell{ppni} = [plotStudyCountCell{ppni} '-point'];
end
plotStudyCountCell{end} = sprintf('%d-to-%d-point',plotPointNvec(end),pointBaseVec(end));

% append continuous scales to the rightmost end
plotStudyCount(:,end+1) = confTypeT.continScaleNum(1:2)';
plotPointNvec(end+1) = plotPointNvec(end)+1;
plotPointNNum = size(plotStudyCount,2);
plotStudyCountCell{end+1} = 'continuous'; 

% perform Z test of two proportions
totalCounts = sum(plotStudyCount,2);
plotStudyP = plotStudyCount./repmat(totalCounts,[1, plotPointNNum]);
for plotPointi = 1:plotPointNNum
    [h, p, z] = twoPropZTest(...
        plotStudyP(1,plotPointi),plotStudyP(2,plotPointi),...
        totalCounts(1),totalCounts(2));
    testRst(plotPointi).h = h;
    testRst(plotPointi).p = p;
    testRst(plotPointi).z = z;
end

%% plot bar graph showing proportion of studies for scale types
ylimvec = [0, 60];
fs = 24;
xticklabels = plotStudyCountCell;

plotx = plotPointNvec';
ploty = 100*(plotStudyCount./repmat(totalCounts,[1 size(plotStudyCount,2)]))';
ispercStr = {...
    sprintf('Perception studies (N = %d)',sum(plotStudyCount(1,:))),...
    sprintf('Memory studies (N = %d)',sum(plotStudyCount(2,:)))};
bh = bar(plotx,ploty);
ylim(ylimvec);
xlabel('Confidence scale type','fontsize',fs);
ylabel('Percent of studies','fontsize',fs);
legend(bh,ispercStr);
set(gca,'xtick',plotPointNvec,'xticklabel',xticklabels,'fontsize',fs);

plotStudyCount

%% Z test for two-sample proportions
function [h,pval,z] = twoPropZTest(p1,p2,n1,n2)
    p = ((n1*p1)+(n2*p2))/(n1+n2);
    z = (p1-p2)./sqrt(p*(1-p)*((1/n1)+(1/n2)));   
    pval = (1-normcdf(abs(z)))*2; % two-tailed
    h = pval<.05;
end