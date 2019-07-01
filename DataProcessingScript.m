% *************************************************************************
% This script will process the s2p files. The issue of a limited number of 
% letter labels has been fixed by using letters function written in the
% letter.m file. Which MUST be in the same folder as this script. There is
% no need to open the letters.m file.
% Last updated 6/21/19
% By Jesus Perez
%**************************************************************************

% define transmission length or length of the sample
length=.031; % length of wave guide in millimeters

% define folder with the raw data

prompt = 'Enter folder address: ';
Address = inputdlg(prompt,'Folder Location',[1 70]);
folder = cell2mat(Address(1,1));
%folder = 'C:\Users\jap11\Desktop\Jesus Perez Summer 2019\Data Analysis Testing';
%oldoutputfile='Name Summary file here.xlsx';

outputfilename = inputdlg('Name your summary file: ','Summary file name',[1 50]);
outputfilename = cell2mat(outputfilename(1,1));
outputfile = [outputfilename,'.xlsx'];

% read the data in the folder
cd(folder)
fnames=dir([folder,'/*.s2p']);
[m,p]=size(fnames);
trdata=cell(1,m);

% write frequency and omega from the first file to the summary file ----->[]

K=1;
    fname=fnames(K).name;
    A=importdata(fname,' ',9);
    data= A.data;
    
dataout(:,1)=data(:,1);
dataout(:,2)=data(:,1)*2*pi;

warning( 'off', 'MATLAB:xlswrite:AddSheet' ); % remove warning!
disp('Beginning to creat the summary file...')
disp(' ')

xlswrite(outputfile,{'frequency'},'AMP_S21','a1');
xlswrite(outputfile,dataout(:,1),'AMP_S21','a2');
xlswrite(outputfile,{'omega'},'AMP_S21','b1');
xlswrite(outputfile,dataout(:,2),'AMP_S21','b2');

xlswrite(outputfile,{'frequency'},'AMP_S11','a1');
xlswrite(outputfile,dataout(:,1),'AMP_S11','a2');
xlswrite(outputfile,{'omega'},'AMP_S11','b1');
xlswrite(outputfile,dataout(:,2),'AMP_S11','b2');

xlswrite(outputfile,{'frequency'},'UP_S21','a1');
xlswrite(outputfile,dataout(:,1),'UP_S21','a2');
xlswrite(outputfile,{'omega'},'UP_S21','b1');
xlswrite(outputfile,dataout(:,2),'UP_S21','b2');

xlswrite(outputfile,{'frequency'},'UP_S11','a1');
xlswrite(outputfile,dataout(:,1),'UP_S11','a2');
xlswrite(outputfile,{'omega'},'UP_S11','b1');
xlswrite(outputfile,dataout(:,2),'UP_S11','b2');

xlswrite(outputfile,{'frequency'},'GD_S21','a1');
xlswrite(outputfile,dataout(:,1),'GD_S21','a2');
xlswrite(outputfile,{'omega'},'GD_S21','b1');
xlswrite(outputfile,dataout(:,2),'GD_S21','b2');

xlswrite(outputfile,{'frequency'},'GD_S11','a1');
xlswrite(outputfile,dataout(:,1),'GD_S11','a2');
xlswrite(outputfile,{'omega'},'GD_S11','b1');
xlswrite(outputfile,dataout(:,2),'GD_S11','b2');
xlswrite(outputfile,dataout(:,2),'GD_S11','b2');

xlswrite(outputfile,{'frequency'},'ATT_CONST','a1');
xlswrite(outputfile,dataout(:,1),'ATT_CONST','a2');
xlswrite(outputfile,{'omega'},'ATT_CONST','b1');
xlswrite(outputfile,dataout(:,2),'ATT_CONST','b2');

xlswrite(outputfile,{'frequency'},'PHASE_CONST','a1');
xlswrite(outputfile,dataout(:,1),'PHASE_CONST','a2');
xlswrite(outputfile,{'omega'},'PHASE_CONST','b1');
xlswrite(outputfile,dataout(:,2),'PHASE_CONST','b2'); 

Array1 = 1:4000;
newLabels = letters(Array1);

disp('Labels have been created. About to make the individual files and add to the Summary file.')
disp(' ')
disp('Patience...')
disp(' ')
fprintf('There are %d files to process... \n',round(m))
disp(' ')
%write frequency and omega from the first file to the summary file []-----> 

for k=1:m

    fname=fnames(k).name;
    A=importdata(fname,' ',9);
    data=A.data;

dataout(:,1)=data(:,1);
dataout(:,2)=data(:,1)*2*pi;
dataout(:,3:10)=data(:,2:9);
% Phase unwrap
S11p=dataout(:,4);
S21p=dataout(:,6);
omega=dataout(:,2);

S11up=unwrap(S11p*pi/180);
S21up=unwrap(S21p*pi/180);
% group delay
S21gd=-diff(S21up)./diff(omega).*1000000;
S11gd=-diff(S11up)./diff(omega).*1000000;
%attenuation constants ------>[]
% Combines magnitude and phase values into complex forms
S11=(dataout(:,3).*(exp(1i.*dataout(:,4)/180.*pi)));
S21=(dataout(:,5).*(exp(1i.*dataout(:,6)/180.*pi)));
S12=(dataout(:,7).*(exp(1i.*dataout(:,8)/180.*pi)));
S22=(dataout(:,9).*(exp(1i.*dataout(:,10)/180.*pi)));
% calculate Delta S and define z0
Delta_S= ((S11.*S22)-(S21.*S12));
z0 = 50;
% convert to abcd form
A= (((1)+ S11-S22-Delta_S)./(2.*S21));
B= (((1+ S11+S22+Delta_S).*z0)./(2.*S21));
C= (1-S11-S22+Delta_S)./(2.*S21.*z0);
D= (1-S11+S22-Delta_S)./(2.*S21);

Propagation_constant = asinh((B.*C).^(1./2))./length;
att_constant=real(Propagation_constant);
phase_constant=imag(Propagation_constant);

% look-up table for output
title={'frequency (Hz)' 'Omega (w)' 'S11'   'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'};

xlswrite([fname,'.xlsx'],title,'Sheet1','a1:r1');
xlswrite([fname,'.xlsx'],dataout,'Sheet1','a2:j802');
xlswrite([fname,'.xlsx'],S11up,'Sheet1','p2:p802');
xlswrite([fname,'.xlsx'],S21up,'Sheet1','l2:l802');
xlswrite([fname,'.xlsx'],S11gd,'Sheet1','r2:r802');
xlswrite([fname,'.xlsx'],S21gd,'Sheet1','n2:n802');
%writting S21 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,dataout(:,5),'AMP_S21',[newLabels{k+2},'2']);
%writting S21 group delay to summmary table
xlswrite(outputfile,{fname},'GD_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21gd,'GD_S21',[newLabels{k+2},'2']);
%writting S21 unwraped phase to summmary table
xlswrite(outputfile,{fname},'UP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21up,'UP_S21',[newLabels{k+2},'2']);
%writting S11 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S11',[newLabels{k+2},'1'])
xlswrite(outputfile,dataout(:,3),'AMP_S11',[newLabels{k+2},'2']);
%writting S11 group delay to summmary table
xlswrite(outputfile,{fname},'GD_S11',[newLabels{k+2},'1']);
xlswrite(outputfile,S11gd,'GD_S11',[newLabels{k+2},'2']);
%writting S11 unwraped phase to summmary table
xlswrite(outputfile,{fname},'UP_S11',[newLabels{k+2},'1']);
xlswrite(outputfile,S11up,'UP_S11',[newLabels{k+2},'2']);
%writting ATT_CONST unwraped phase to summmary table
xlswrite(outputfile,{fname},'ATT_CONST',[newLabels{k+2},'1']);
xlswrite(outputfile,att_constant,'ATT_CONST',[newLabels{k+2},'2']);
%writting PHASE_CONST unwraped phase to summmary table
xlswrite(outputfile,{fname},'PHASE_CONST',[newLabels{k+2},'1']);
xlswrite(outputfile,phase_constant,'PHASE_CONST',[newLabels{k+2},'2']);

c = k;
fprintf(' %d sheets have been added to the Summary Folder \n',round(c))

end

disp(' ')
disp('Done!')
disp(' ')
