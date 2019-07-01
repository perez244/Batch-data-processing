%**************************************************************************
% This script is meant to pull selected frequencies and selected parameters
% from all the files in the summary file created with processing script. This 
% script will output everything to 'Analysis for + "Summary file name".xlsx'
% Last edited 6/21/19
% By Jesus Perez 
%**************************************************************************
clc
clear
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% This chunk of code will ask you to enter(paste)the folder address & then
% select the Summary file you'd like to pull from. 
prompt = 'Enter folder address: ';
dd = inputdlg(prompt,'Folder Location',[1 70]);
dd = cell2mat(dd(1,1));
cd(dd) 
filename = uigetfile('*.xlsx');
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% display message on the command window
disp(' ')
disp('File selected!')
disp(' ')
disp('Importing data..')
disp(' ')

%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Import data from summary file
A=importdata(filename,' ',9);
data= A.data;
ColumnTitles = A.colheaders.AMP_S21; % Weird way to store the names like 'frequencies' and individual file names for later but it works...
ColumnTitles(:,2) =[];

fileName = filename(1:end-5);% chop of the .xlsx
outputfile= ['Analysis file for ',fileName,'.xlsx']; % Name the output analysis file
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Allocate data memeory without omega, in the matlab workspace 
Amp_S21 = data.AMP_S21;
Amp_S21(:,2) = [];
Amp_S11 = data.AMP_S11;
Amp_S11(:,2) = [];
Up_S21 = data.UP_S21;
Up_S21(:,2) = [];
Up_S11 = data.UP_S11;  
Up_S11(:,2) = [];
GD_S21 = data.GD_S21;
GD_S21(:,2) = [];
GD_S11 = data.GD_S11;
GD_S11(:,2) = [];
Att_Const = data.ATT_CONST;
Att_Const(:,2) = [];
Phase_Const = data.PHASE_CONST;
Phase_Const(:,2) = [];
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Allow user to select which frequencies they'd like to pull
dataout(:,1)=data.AMP_S21(:,1);
Frequencies = dataout(:,1);
[rowLocation,tf] = listdlg('PromptString','Select frequecies','ListString',num2str(Frequencies));
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% display message on the command window
disp(' ')
disp('Thank you')
disp(' ')
disp('Building analysis file... ')
disp(' ')
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% for loop that matches frequency indices to frequency values.
% For example, frequency 3 has a value 409900000(Hz). (Depending on the
% insturment it came from).
chosenFrequencies = [];
idx = 1;
for R = rowLocation
    
chosenFrequencies(1,idx) = Frequencies(R,1);
idx = idx +1;

end
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Allocate memeory for transposed data
trans_Amp_S21 = [];
trans_Amp_S11 = [];
trans_Up_S21 = [];
trans_Up_S11 = [];
trans_GD_S21 = [];
trans_GD_S11 = [];
trans_Att_Const = [];
trans_Phase_Const = [];
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% The line below will create a new struct(Method that matlab uses to store data) without omega
MatlabData = struct('Amp_S21',Amp_S21,'Amp_S11',Amp_S11,'Up_S21',Up_S21,'Up_S11',Up_S11,'GD_S21',GD_S21,'GD_S11',GD_S11,'Att_Const',Att_Const,'Phase_Const',Phase_Const);
% The line below will create a new struct(Method that matlab uses to store data) to fill with the transposed data
TransposedMatlabData = struct('trans_Amp_S21',trans_Amp_S21,'trans_Amp_S11',trans_Amp_S11,'trans_Up_S21',trans_Up_S21,'trans_Up_S11',trans_Up_S11,'trans_GD_S21',trans_GD_S21,'trans_GD_S11',trans_GD_S11,'trans_Att_Const',trans_Att_Const,'trans_Phase_Const',trans_Phase_Const);
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%The field names are the parameters. Like S21, GD, Phase etc
FN = fieldnames(MatlabData);
fn = fieldnames(TransposedMatlabData);
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Ask the user to select parameters they'd like to look at and loop
% loop to match parameters indices to parameter values(for example,S21)
[parametersIndex,tf2] = listdlg('PromptString','Select parameters','ListString',FN);

chosenParameters = cell(1,length(parametersIndex));
z = 1;
for  r = parametersIndex
    
    chosenParameters{1,z} = FN{r,1};
    z = z +1;

end
chosenParameters = chosenParameters';
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% loop to get the transposed parameter names...for later
neededTransParameters = cell(1,length(parametersIndex));
y = 1;
for x = parametersIndex
    
    neededTransParameters{1,y} = fn{x,1};
    y = y + 1;
    
end
neededTransParameters = neededTransParameters';
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% transpose the data!
for U=1:numel(neededTransParameters)
    if( isnumeric(TransposedMatlabData.(neededTransParameters{U})) )
        P = length(ColumnTitles);%length(chosenFrequencies);
        for j = 1:P %loop through column
            c2 = 1;
            for j2 = rowLocation %loop through row
                TransposedMatlabData.(neededTransParameters{U})(j,c2) = cell2mat({MatlabData.(chosenParameters{U})(j2,j)});
                c2 = c2+1;
            end
        end
    end
end
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
warning( 'off', 'MATLAB:xlswrite:AddSheet' ); % remove warning!
disp('')
disp('Creating Analysis file..')
disp(' ')
% Write different sheets/tabs in the Analysis Summary file for the chosen
% parameters
for q =1:numel(neededTransParameters)
    
    xlswrite(outputfile,TransposedMatlabData.(neededTransParameters{q}),chosenParameters{q, 1},'B1')

end
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%for loop for inputting "transposed" labels and frequencies
ColumnTitles = ColumnTitles';
for w=1:numel(neededTransParameters)
        
    xlswrite(outputfile,ColumnTitles,chosenParameters{w,1},'A1')
 
end


disp('Done!')
disp(' ')
