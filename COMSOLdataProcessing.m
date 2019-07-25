% *************************************************************************
% This script takes the scattering parameters from the .xlsx sheets
% produced by the "importingCOMOSOLdata" script and process them to be
% placed in a summary file.( Just like s2p proccessing script)
% Last updated 7/9/19
% By Jesus Perez
%**************************************************************************
% define transmission length or length of the sample
clear
clc

length=.036; % length of wave guide in millimeters

Array1 = 1:4000;
newLabels = letters(Array1);
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

% write frequency and omega from the first file to the summary file ----->[]
K=1;
    fname=fnames(K).name;
    A=importdata(fname,' ',1);
    data= A.data;
    
dataout(:,1)=(data(:,1)).*(1e9);
dataout(:,2)=(data(:,1)*2*pi).*(1e9);

warning( 'off', 'MATLAB:xlswrite:AddSheet' ); % remove warning!
disp('Beginning to create the summary file...')
disp(' ')
% 
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

%create RLGC sheets in the summary folder 
xlswrite(outputfile,{'frequency'},'R','a1');
xlswrite(outputfile,dataout(:,1),'R','a2');
xlswrite(outputfile,{'omega'},'R','b1');
xlswrite(outputfile,dataout(:,2),'R','b2'); 

xlswrite(outputfile,{'frequency'},'L','a1');
xlswrite(outputfile,dataout(:,1),'L','a2');
xlswrite(outputfile,{'omega'},'L','b1');
xlswrite(outputfile,dataout(:,2),'L','b2'); 

xlswrite(outputfile,{'frequency'},'G','a1');
xlswrite(outputfile,dataout(:,1),'G','a2');
xlswrite(outputfile,{'omega'},'G','b1');
xlswrite(outputfile,dataout(:,2),'G','b2'); 

xlswrite(outputfile,{'frequency'},'C','a1');
xlswrite(outputfile,dataout(:,1),'C','a2');
xlswrite(outputfile,{'omega'},'C','b1');
xlswrite(outputfile,dataout(:,2),'C','b2');

disp('Labels have been created. About to make the individual files and add to the Summary file.')
disp(' ')
disp('Patience...')
disp(' ')
fprintf('There are %d files to process... \n',round(m))
disp(' ')

for k=1:m

    fname=fnames(k).name;
    A=importdata(fname,' ',1);
    data=A.data;

DataOut(:,1)=(data(:,1)).*(1e9);
DataOut(:,2)=(data(:,1)*2*pi).*(1e9);
w = data(:,1)*2*pi;
DataOut(:,3:10)=data(:,2:9);

% Phase unwrap
S11p=data(:,3);
S21p=data(:,7);
omega=w;
S11up=unwrap(S11p*pi/180);
S21up=unwrap(S21p*pi/180);

% group delay
S21gd=-diff(S21up)./diff(omega).*1000000;
S11gd=-diff(S11up)./diff(omega).*1000000;

%attenuation constants ------>[]
% Combines magnitude and phase values into complex forms
S11=(data(:,2).*(exp(1i.*data(:,3)/180.*pi)));
S21=(data(:,4).*(exp(1i.*data(:,5)/180.*pi)));
S12=(data(:,6).*(exp(1i.*data(:,7)/180.*pi)));
S22=(data(:,8).*(exp(1i.*data(:,9)/180.*pi)));

% calculate Delta S and define z0
Delta_S= ((S11.*S22)-(S21.*S12));
z0 = 50;

% % convert to abcd form
A= (((1)+ S11-S22-Delta_S)./(2.*S21));
B= (((1+ S11+S22+Delta_S).*z0)./(2.*S21));
C= (1-S11-S22+Delta_S)./(2.*S21.*z0);
D= (1-S11+S22-Delta_S)./(2.*S21);
% 
Propagation_constant = asinh((B.*C).^(1./2))./length;
att_constant=real(Propagation_constant);
phase_constant=imag(Propagation_constant);
% 
% %******************************The RLGC values*****************************
res_top = ((1+S11).*(1+S22)) - ((S12).*(S21));
res_bottom = ((2).*(S21));
res = (res_top)./(res_bottom);
resistance  = (z0).*(res);
R = real(resistance);

induct_top = ((1+S11).*(1+S22)) - ((S12).*(S21));
induct_bottom = (2).*(S21);
induct = (induct_top)./(induct_bottom);
inductance = (z0).*(induct);
L = (imag(inductance))./(w);

G_top = (1+S11).*(1-S22) +  ( ((S12).*(S21)) - ((2).*(S21)) );
G_bottom  =  (z0).*(((1+S11).*(1+S22)) - ((S12).*(S21)));
Gtivity = (G_top)./(G_bottom);
G = real(Gtivity);

C_top = ((1+S11).*(1-S22) +  ( ((S12).*(S21)) - ((2).*(S21)) ));
C_bottom  = ((z0).*(((1+S11).*(1+S22)) - ((S12).*(S21))));
Capacity  = (C_top)./(C_bottom);
Cap = (imag(Capacity))./(w);
% %**************************************************************************

% look-up table for output
title={'frequency (Hz)' 'Omega (w)' 'S11'   'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'  ' '  'R'  'L'  'G'  'C'};

lim =  numel(S12)+1;
lim = num2str(lim);

xlswrite([fname,'.xlsx'],title,'Sheet1',['a1:w' lim]);
xlswrite([fname,'.xlsx'],DataOut,'Sheet1',['a2:j' lim]);
xlswrite([fname,'.xlsx'],S11up,'Sheet1',['p2:p' lim]);
xlswrite([fname,'.xlsx'],S21up,'Sheet1',['l2:l' lim]);
xlswrite([fname,'.xlsx'],S11gd,'Sheet1',['r2:r' lim]);
xlswrite([fname,'.xlsx'],S21gd,'Sheet1',['n2:n' lim]);
% add RLGC values to individual files
xlswrite([fname,'.xlsx'],R,'Sheet1',['t2:t' lim]);
xlswrite([fname,'.xlsx'],L,'Sheet1',['u2:u' lim]);
xlswrite([fname,'.xlsx'],G,'Sheet1',['v2:v' lim]);
xlswrite([fname,'.xlsx'],Cap,'Sheet1',['w2:w' lim]);

%writting S21 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,DataOut(:,7),'AMP_S21',[newLabels{k+2},'2']);
%writting S21 group delay to summmary table
xlswrite(outputfile,{fname},'GD_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21gd,'GD_S21',[newLabels{k+2},'2']);
%writting S21 unwraped phase to summmary table
xlswrite(outputfile,{fname},'UP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21up,'UP_S21',[newLabels{k+2},'2']);
%writting S11 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S11',[newLabels{k+2},'1'])
xlswrite(outputfile,DataOut(:,3),'AMP_S11',[newLabels{k+2},'2']);
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

%write the RLGC values
xlswrite(outputfile,{fname},'R',[newLabels{k+2},'1']);
xlswrite(outputfile,R,'R',[newLabels{k+2},'2']);
xlswrite(outputfile,{fname},'L',[newLabels{k+2},'1']);
xlswrite(outputfile,L,'L',[newLabels{k+2},'2']);
xlswrite(outputfile,{fname},'G',[newLabels{k+2},'1']);
xlswrite(outputfile,G,'G',[newLabels{k+2},'2']);
xlswrite(outputfile,{fname},'C',[newLabels{k+2},'1']);
xlswrite(outputfile,Cap,'C',[newLabels{k+2},'2']);

c = k;
fprintf(' %d sheets have been added to the Summary Folder \n',round(c))

 end
% 
datetime.setDefaultFormats('default','hh:mm a MM/dd/yyyy ')
timeStamp = datetime;
disp(' ')
fprintf('Finished at: %s \n',timeStamp)
disp(' ')


