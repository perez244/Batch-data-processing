% *************************************************************************
% This script will process the s2p files. The issue of a limited number of 
% letter labels has been fixed by using letters function written in the
% letter.m file. Which MUST be in the same folder as this script. There is
% no need to open the letters.m file.
% Last updated 6/21/19
% By Jesus Perez
%*****************************************************************************************************
clear
clc
% define transmission length or length of the sample

length=.036; % length of wave guide in meters
uo = pi*4e-7; % permeability of free space (H/m)
ur = 1; %relative permeability of copper conductor
s =  5.99e7; % conductivity of copper (S/m)
g =   0.0150; % loss tangent for Si @ 10GHz: https://www.researchgate.net/publication/315830159_Study_the_loss_of_microstrip_on_silicon
E = 11.7; % dielectic constant for Silicon:  http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
%*****************************************************************************************************
% create the excel column labels 
Array1 = 1:4000;
newLabels = letters(Array1);
%**************************************************************************
% ask the user for the folder address, a name for the summary folder and
% the height of the average surface roughness in microns
x = inputdlg({'Folder Address','Summary file name','Roughness height (microns)'},'User inputs', [1 100; 1 50; 1 30]);
folder =  x{1};
outputfilename = x{2};
outputfile = [outputfilename,'.xlsx'];
rms = str2double(x{3}); % average copper roughness
rms = rms * 1e-6; % convert micron to meters
%-------------------------------------------------------------------------------------------------------------------------------------------
%read and find out how many s2p files are in the folder
cd(folder)
fnames=dir([folder,'/*.s2p']);
[m,p]=size(fnames);
trdata=cell(1,m);
%-------------------------------------------------------------------------------------------------------------------------------------------
% import the data
K=1;
    fname=fnames(K).name;
    A=importdata(fname,' ',9);
    data= A.data;
    
dataout(:,1)=data(:,1);
dataout(:,2)=data(:,1)*2*pi;
%-------------------------------------------------------------------------------------------------------------------------------------------
% Dialog box to ask the user how many files to put in each summary folder.
% or comment out the 3 lines below and enter in {b = enter number here} line;
prompt = {'How many files would you like per Summary folder? '};
b = inputdlg(prompt,'Number of files per summary folder',[1 50]);
b = str2double(b); % number of files per sheet


a = m; % number of .s2p files we will be crunching `
r = rem(a,b);
E = a - r;
N = E/b; % Number of summary files that will be created

%create Summary file names. File 1 through file N
FileIndex = 1:N;
fileNames = cell(1,N);
for idx = 1:N
    fileNames(1,idx)={[outputfilename,' ',int2str(FileIndex(idx)),'.xlsx']};    
end

warning( 'off', 'MATLAB:xlswrite:AddSheet' ); % remove warning!
fprintf(' There will be %d Summary folders. \n',round(N))
disp(' ')
disp('Beginning to create the summary folders...')
disp(' ')
%-------------------------------------------------------------------------------------------------------------------------------------------
% create and set up the summary folders
jdx = 1;
% input the frequencies and omega into all the summary files
for L =  1:N
    
    outputfile = char(fileNames(1,L));

    xlswrite(outputfile,{'frequency'},'AMP_S21','a1');
    xlswrite(outputfile,dataout(:,1),'AMP_S21','a2');
    xlswrite(outputfile,{'omega'},'AMP_S21','b1');
    xlswrite(outputfile,dataout(:,2),'AMP_S21','b2');

                if(L == 1)
                     fprintf(' %d Summary File has been created \n',round(L))
                else
                    fprintf(' %d Summary Files have been created \n',round(L))
                end


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

    xlswrite(outputfile,{'frequency'},'ATT_CONST_dB','a1');
    xlswrite(outputfile,dataout(:,1),'ATT_CONST_dB','a2');
    xlswrite(outputfile,{'omega'},'ATT_CONST_dB','b1');
    xlswrite(outputfile,dataout(:,2),'ATT_CONST_dB','b2');

    xlswrite(outputfile,{'frequency'},'CORRECTED_ATT_CONST_dB','a1');
    xlswrite(outputfile,dataout(:,1),'CORRECTED_ATT_CONST_dB','a2');
    xlswrite(outputfile,{'omega'},'CORRECTED_ATT_CONST_dB','b1');
    xlswrite(outputfile,dataout(:,2),'CORRECTED_ATT_CONST_dB','b2');

    xlswrite(outputfile,{'frequency'},'PHASE_CONST_dB','a1');
    xlswrite(outputfile,dataout(:,1),'PHASE_CONST_dB','a2');
    xlswrite(outputfile,{'omega'},'PHASE_CONST_dB','b1');
    xlswrite(outputfile,dataout(:,2),'PHASE_CONST_dB','b2'); 

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

    xlswrite(outputfile,{'frequency'},'dB_a_diel','a1');
    xlswrite(outputfile,dataout(:,1),'dB_a_diel','a2');
    xlswrite(outputfile,{'omega'},'dB_a_diel','b1');
    xlswrite(outputfile,dataout(:,2),'dB_a_diel','b2'); 

    xlswrite(outputfile,{'frequency'},'dB_corrected_metal_S21','a1');
    xlswrite(outputfile,dataout(:,1),'dB_corrected_metal_S21','a2');
    xlswrite(outputfile,{'omega'},'dB_corrected_metal_S21','b1');
    xlswrite(outputfile,dataout(:,2),'dB_corrected_metal_S21','b2');

    xlswrite(outputfile,{'frequency'},'New_total_corrected_metal_S21','a1');
    xlswrite(outputfile,dataout(:,1),'New_total_corrected_metal_S21','a2');
    xlswrite(outputfile,{'omega'},'New_total_corrected_metal_S21','b1');
    xlswrite(outputfile,dataout(:,2),'New_total_corrected_metal_S21','b2');
   
end


Array1 = 1:4000;
newLabels = letters(Array1);

disp(' ')
disp('Labels have been created. About to make the individual files and add to the Summary files.')
disp(' ')
disp('Patience...')
disp(' ')
fprintf('There are %d s2p files to process... \n',round(m))
disp(' ')
%-------------------------------------------------------------------------------------------------------------------------------------------
% input data into their respective summary folder
index = 1;
Y = 1;
for k=1:m

    if(rem(k,b) == 0)% 
        %check
        if(index == N)% N=10
            outputfile = char(fileNames(1,index));
                fname=fnames(k).name;
                A=importdata(fname,' ',9);
                data=A.data;

                dataout(:,1)=data(:,1);
                freq = dataout(:,1);
                dataout(:,2)=data(:,1)*2*pi;
                w = data(:,1)*2*pi;
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
                A = (((1+S11).*(1-S22)) + ((S12).*(S21)))./(2.*S21);
                B = (z0.*(((1+S11).*(1+S22)) - ((S12).*(S21))))./(1.*S21);
                C = (((1-S11).*(1-S22)) - ((S12).*(S21)))./(z0.*2.*S21);
                D = (((1-S11).*(1+S22)) + ((S12).*(S21)))./(2.*S21);

                Propagation_constant = asinh((B.*C).^(1./2))./length; % progagation constant
                total_att_constant=real(Propagation_constant); % alpha
                phase_constant=imag(Propagation_constant); % beta

                total_att_constant = -20.*log10(total_att_constant); % alpha in dB
                phase_constant_dB = -20.*log10(phase_constant); % beta in dB

                % correct the attentuation constant with the roughness coefficient
                SKD = 1./sqrt(s*pi*uo*ur*freq); % skin depth
                x = (rms./SKD).^2;
                Ksr = 1+(2*atan(1.4*x))./pi; % copper surface roughness correction coefficent

                a_diel = 2.3* freq*tan(g)*sqrt(E);% dielectric attentuation

                a_con = total_att_constant-a_diel; % conductor attentuation without dielectric portion (intel paper)
                a_cr = a_con.*Ksr; %(intel)

                Corrected_ATT = total_att_constant.* Ksr; % correction attentuation

                %**************************************************************************
                % New equations for S21 correction
                % method of dB-ing them individually then adding them together at the end
                % theorized based of Intel paper

                linear_a_diel = 2.3* freq*tan(g)*sqrt(E);% dielectric attentuation
                dB_a_diel = -20.*log10(linear_a_diel);

                linear_corrected_metal_S21 = (a_cr.*8.5850);
                dB_corrected_metal_S21 = -20.*log10(linear_corrected_metal_S21);

                New_total_corrected_metal_S21 = dB_a_diel + dB_corrected_metal_S21;

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
                titles={'frequency (Hz)' 'Omega (w)' 'S11'  'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'  ' '  'R'  'L'  'G'  'C' ' ' 'Fixed_ATT_CONST '};

                lim =  numel(S12)+1;
                lim = num2str(lim);

                xlswrite([fname,'.xlsx'],titles,'Sheet1',['a1:y' lim]);
                xlswrite([fname,'.xlsx'],dataout,'Sheet1',['a2:j' lim]);
                xlswrite([fname,'.xlsx'],S11up,'Sheet1',['p2:p' lim]);
                xlswrite([fname,'.xlsx'],S21up,'Sheet1',['l2:l' lim]);
                xlswrite([fname,'.xlsx'],S11gd,'Sheet1',['r2:r' lim]);
                xlswrite([fname,'.xlsx'],S21gd,'Sheet1',['n2:n' lim]);
                % add RLGC values to individual files
                xlswrite([fname,'.xlsx'],R,'Sheet1',['t2:t' lim]);
                xlswrite([fname,'.xlsx'],L,'Sheet1',['u2:u' lim]);
                xlswrite([fname,'.xlsx'],G,'Sheet1',['v2:v' lim]);
                xlswrite([fname,'.xlsx'],Cap,'Sheet1',['w2:w' lim]);
                xlswrite([fname,'.xlsx'],Corrected_ATT,'Sheet1',['y2:y' lim]);

                %writting S21 amplitude to summmary table
                xlswrite(outputfile,{fname},'AMP_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,dataout(:,7),'AMP_S21',[newLabels{Y+2},'2']);
                %writting S21 group delay to summmary table
                xlswrite(outputfile,{fname},'GD_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S21gd,'GD_S21',[newLabels{Y+2},'2']);
                %writting S21 unwraped phase to summmary table
                xlswrite(outputfile,{fname},'UP_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S21up,'UP_S21',[newLabels{Y+2},'2']);
                %writting S11 amplitude to summmary table
                xlswrite(outputfile,{fname},'AMP_S11',[newLabels{Y+2},'1'])
                xlswrite(outputfile,dataout(:,3),'AMP_S11',[newLabels{Y+2},'2']);
                %writting S11 group delay to summmary table
                xlswrite(outputfile,{fname},'GD_S11',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S11gd,'GD_S11',[newLabels{Y+2},'2']);
                %writting S11 unwraped phase to summmary table
                xlswrite(outputfile,{fname},'UP_S11',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S11up,'UP_S11',[newLabels{Y+2},'2']);
                %writting ATT_CONST unwraped phase to summmary table
                xlswrite(outputfile,{fname},'ATT_CONST_dB',[newLabels{Y+2},'1']);
                xlswrite(outputfile,total_att_constant,'ATT_CONST_dB',[newLabels{Y+2},'2']);
                %writting Corrected atten_CONST unwraped phase to summmary table

                xlswrite(outputfile,{fname},'CORRECTED_ATT_CONST_dB',[newLabels{Y+2},'1']);
                xlswrite(outputfile,Corrected_ATT,'CORRECTED_ATT_CONST_dB',[newLabels{Y+2},'2']);

                %writting PHASE_CONST unwraped phase to summmary table
                xlswrite(outputfile,{fname},'PHASE_CONST_dB',[newLabels{Y+2},'1']);
                xlswrite(outputfile,phase_constant_dB,'PHASE_CONST_dB',[newLabels{Y+2},'2']);
                %write the RLGC values
                xlswrite(outputfile,{fname},'R',[newLabels{Y+2},'1']);
                xlswrite(outputfile,R,'R',[newLabels{Y+2},'2']);
                xlswrite(outputfile,{fname},'L',[newLabels{Y+2},'1']);
                xlswrite(outputfile,L,'L',[newLabels{Y+2},'2']);
                xlswrite(outputfile,{fname},'G',[newLabels{Y+2},'1']);
                xlswrite(outputfile,G,'G',[newLabels{Y+2},'2']);
                xlswrite(outputfile,{fname},'C',[newLabels{Y+2},'1']);
                xlswrite(outputfile,Cap,'C',[newLabels{Y+2},'2']);
                % xlswrite(outputfile,{fname},'a_diel',[newLabels{Y+2},'1']);
                % xlswrite(outputfile,a_diel,'a_diel',[newLabels{Y+2},'2']);
                % xlswrite(outputfile,{fname},'corrected_metal_S21_dB',[newLabels{Y+2},'1']);
                % xlswrite(outputfile,corrected_metal_S21_dB,'corrected_metal_S21_dB',[newLabels{Y+2},'2']);
                % xlswrite(outputfile,{fname},'total_corrected_metal_S21',[newLabels{Y+2},'1']);
                % xlswrite(outputfile,total_corrected_metal_S21,'total_corrected_metal_S21',[newLabels{Y+2},'2']);

                xlswrite(outputfile,{fname},'dB_a_diel',[newLabels{Y+2},'1']);
                xlswrite(outputfile,dB_a_diel,'dB_a_diel',[newLabels{Y+2},'2']);

                xlswrite(outputfile,{fname},'dB_corrected_metal_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,dB_corrected_metal_S21,'dB_corrected_metal_S21',[newLabels{Y+2},'2']);

                xlswrite(outputfile,{fname},'New_total_corrected_metal_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,New_total_corrected_metal_S21,'New_total_corrected_metal_S21',[newLabels{Y+2},'2']);

                    c = k;
                           if (k == 1)
                               fprintf(' %d sheet has been added to its assigned summary folder \n',round(c))
                           else
                               fprintf(' %d sheets have been added to their assigned summary folder \n',round(c))
                           end
                 Y = Y+1;
                           
        else
            outputfile = char(fileNames(1,index));
                fname=fnames(k).name;
                A=importdata(fname,' ',9);
                data=A.data;

                dataout(:,1)=data(:,1);
                freq = dataout(:,1);
                dataout(:,2)=data(:,1)*2*pi;
                w = data(:,1)*2*pi;
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
                A = (((1+S11).*(1-S22)) + ((S12).*(S21)))./(2.*S21);
                B = (z0.*(((1+S11).*(1+S22)) - ((S12).*(S21))))./(1.*S21);
                C = (((1-S11).*(1-S22)) - ((S12).*(S21)))./(z0.*2.*S21);
                D = (((1-S11).*(1+S22)) + ((S12).*(S21)))./(2.*S21);

                Propagation_constant = asinh((B.*C).^(1./2))./length; % progagation constant
                total_att_constant=real(Propagation_constant); % alpha
                phase_constant=imag(Propagation_constant); % beta

                total_att_constant = -20.*log10(total_att_constant); % alpha in dB
                phase_constant_dB = -20.*log10(phase_constant); % beta in dB

                % correct the attentuation constant with the roughness coefficient
                SKD = 1./sqrt(s*pi*uo*ur*freq); % skin depth
                x = (rms./SKD).^2;
                Ksr = 1+(2*atan(1.4*x))./pi; % copper surface roughness correction coefficent

                a_diel = 2.3* freq*tan(g)*sqrt(E);% dielectric attentuation

                a_con = total_att_constant-a_diel; % conductor attentuation without dielectric portion (intel paper)
                a_cr = a_con.*Ksr; %(intel)

                Corrected_ATT = total_att_constant.* Ksr; % correction attentuation

                %**************************************************************************
                % New equations for S21 correction
                % method of dB-ing them individually then adding them together at the end
                % theorized based of Intel paper

                linear_a_diel = 2.3* freq*tan(g)*sqrt(E);% dielectric attentuation
                dB_a_diel = -20.*log10(linear_a_diel);

                linear_corrected_metal_S21 = (a_cr.*8.5850);
                dB_corrected_metal_S21 = -20.*log10(linear_corrected_metal_S21);

                New_total_corrected_metal_S21 = dB_a_diel + dB_corrected_metal_S21;

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
                titles={'frequency (Hz)' 'Omega (w)' 'S11'  'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'  ' '  'R'  'L'  'G'  'C' ' ' 'Fixed_ATT_CONST '};

                lim =  numel(S12)+1;
                lim = num2str(lim);

                xlswrite([fname,'.xlsx'],titles,'Sheet1',['a1:y' lim]);
                xlswrite([fname,'.xlsx'],dataout,'Sheet1',['a2:j' lim]);
                xlswrite([fname,'.xlsx'],S11up,'Sheet1',['p2:p' lim]);
                xlswrite([fname,'.xlsx'],S21up,'Sheet1',['l2:l' lim]);
                xlswrite([fname,'.xlsx'],S11gd,'Sheet1',['r2:r' lim]);
                xlswrite([fname,'.xlsx'],S21gd,'Sheet1',['n2:n' lim]);
                % add RLGC values to individual files
                xlswrite([fname,'.xlsx'],R,'Sheet1',['t2:t' lim]);
                xlswrite([fname,'.xlsx'],L,'Sheet1',['u2:u' lim]);
                xlswrite([fname,'.xlsx'],G,'Sheet1',['v2:v' lim]);
                xlswrite([fname,'.xlsx'],Cap,'Sheet1',['w2:w' lim]);
                xlswrite([fname,'.xlsx'],Corrected_ATT,'Sheet1',['y2:y' lim]);

                %writting S21 amplitude to summmary table
                xlswrite(outputfile,{fname},'AMP_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,dataout(:,7),'AMP_S21',[newLabels{Y+2},'2']);
                %writting S21 group delay to summmary table
                xlswrite(outputfile,{fname},'GD_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S21gd,'GD_S21',[newLabels{Y+2},'2']);
                %writting S21 unwraped phase to summmary table
                xlswrite(outputfile,{fname},'UP_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S21up,'UP_S21',[newLabels{Y+2},'2']);
                %writting S11 amplitude to summmary table
                xlswrite(outputfile,{fname},'AMP_S11',[newLabels{Y+2},'1'])
                xlswrite(outputfile,dataout(:,3),'AMP_S11',[newLabels{Y+2},'2']);
                %writting S11 group delay to summmary table
                xlswrite(outputfile,{fname},'GD_S11',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S11gd,'GD_S11',[newLabels{Y+2},'2']);
                %writting S11 unwraped phase to summmary table
                xlswrite(outputfile,{fname},'UP_S11',[newLabels{Y+2},'1']);
                xlswrite(outputfile,S11up,'UP_S11',[newLabels{Y+2},'2']);
                %writting ATT_CONST unwraped phase to summmary table
                xlswrite(outputfile,{fname},'ATT_CONST_dB',[newLabels{Y+2},'1']);
                xlswrite(outputfile,total_att_constant,'ATT_CONST_dB',[newLabels{Y+2},'2']);
                %writting Corrected atten_CONST unwraped phase to summmary table

                xlswrite(outputfile,{fname},'CORRECTED_ATT_CONST_dB',[newLabels{Y+2},'1']);
                xlswrite(outputfile,Corrected_ATT,'CORRECTED_ATT_CONST_dB',[newLabels{Y+2},'2']);

                %writting PHASE_CONST unwraped phase to summmary table
                xlswrite(outputfile,{fname},'PHASE_CONST_dB',[newLabels{Y+2},'1']);
                xlswrite(outputfile,phase_constant_dB,'PHASE_CONST_dB',[newLabels{Y+2},'2']);
                %write the RLGC values
                xlswrite(outputfile,{fname},'R',[newLabels{Y+2},'1']);
                xlswrite(outputfile,R,'R',[newLabels{Y+2},'2']);
                xlswrite(outputfile,{fname},'L',[newLabels{Y+2},'1']);
                xlswrite(outputfile,L,'L',[newLabels{Y+2},'2']);
                xlswrite(outputfile,{fname},'G',[newLabels{Y+2},'1']);
                xlswrite(outputfile,G,'G',[newLabels{Y+2},'2']);
                xlswrite(outputfile,{fname},'C',[newLabels{Y+2},'1']);
                xlswrite(outputfile,Cap,'C',[newLabels{Y+2},'2']);
                % xlswrite(outputfile,{fname},'a_diel',[newLabels{Y+2},'1']);
                % xlswrite(outputfile,a_diel,'a_diel',[newLabels{Y+2},'2']);
                % xlswrite(outputfile,{fname},'corrected_metal_S21_dB',[newLabels{Y+2},'1']);
                % xlswrite(outputfile,corrected_metal_S21_dB,'corrected_metal_S21_dB',[newLabels{Y+2},'2']);
                % xlswrite(outputfile,{fname},'total_corrected_metal_S21',[newLabels{Y+2},'1']);
                % xlswrite(outputfile,total_corrected_metal_S21,'total_corrected_metal_S21',[newLabels{Y+2},'2']);

                xlswrite(outputfile,{fname},'dB_a_diel',[newLabels{Y+2},'1']);
                xlswrite(outputfile,dB_a_diel,'dB_a_diel',[newLabels{Y+2},'2']);

                xlswrite(outputfile,{fname},'dB_corrected_metal_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,dB_corrected_metal_S21,'dB_corrected_metal_S21',[newLabels{Y+2},'2']);

                xlswrite(outputfile,{fname},'New_total_corrected_metal_S21',[newLabels{Y+2},'1']);
                xlswrite(outputfile,New_total_corrected_metal_S21,'New_total_corrected_metal_S21',[newLabels{Y+2},'2']);

                    c = k;
                           if (k == 1)
                               fprintf(' %d sheet has been added to its assigned summary folder \n',round(c))
                           else
                               fprintf(' %d sheets have been added to their assigned summary folder \n',round(c))
                           end
                    index = index +1;
                    Y = 1;
        end
        
    else
        outputfile = char(fileNames(1,index));
        
        fname=fnames(k).name;
                A=importdata(fname,' ',9);
                data=A.data;

                dataout(:,1)=data(:,1);
                freq = dataout(:,1);
                dataout(:,2)=data(:,1)*2*pi;
                w = data(:,1)*2*pi;
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
                    A = (((1+S11).*(1-S22)) + ((S12).*(S21)))./(2.*S21);
                    B = (z0.*(((1+S11).*(1+S22)) - ((S12).*(S21))))./(1.*S21);
                    C = (((1-S11).*(1-S22)) - ((S12).*(S21)))./(z0.*2.*S21);
                    D = (((1-S11).*(1+S22)) + ((S12).*(S21)))./(2.*S21);

                    Propagation_constant = asinh((B.*C).^(1./2))./length; % progagation constant
                    total_att_constant=real(Propagation_constant); % alpha
                    phase_constant=imag(Propagation_constant); % beta

                    total_att_constant = -20.*log10(total_att_constant); % alpha in dB
                    phase_constant_dB = -20.*log10(phase_constant); % beta in dB

                    % correct the attentuation constant with the roughness coefficient
                    SKD = 1./sqrt(s*pi*uo*ur*freq); % skin depth
                    x = (rms./SKD).^2;
                    Ksr = 1+(2*atan(1.4*x))./pi; % copper surface roughness correction coefficent

                    a_diel = 2.3* freq*tan(g)*sqrt(E);% dielectric attentuation

                    a_con = total_att_constant-a_diel; % conductor attentuation without dielectric portion (intel paper)
                    a_cr = a_con.*Ksr; %(intel)

                    Corrected_ATT = total_att_constant.* Ksr; % correction attentuation

                    %**************************************************************************
                    % New equations for S21 correction
                    % method of dB-ing them individually then adding them together at the end
                    % theorized based of Intel paper

                    linear_a_diel = 2.3* freq*tan(g)*sqrt(E);% dielectric attentuation
                    dB_a_diel = -20.*log10(linear_a_diel);

                    linear_corrected_metal_S21 = (a_cr.*8.5850);
                    dB_corrected_metal_S21 = -20.*log10(linear_corrected_metal_S21);

                    New_total_corrected_metal_S21 = dB_a_diel + dB_corrected_metal_S21;

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
                    titles={'frequency (Hz)' 'Omega (w)' 'S11'  'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'  ' '  'R'  'L'  'G'  'C' ' ' 'Fixed_ATT_CONST '};

                    lim =  numel(S12)+1;
                    lim = num2str(lim);

                    xlswrite([fname,'.xlsx'],titles,'Sheet1',['a1:y' lim]);
                    xlswrite([fname,'.xlsx'],dataout,'Sheet1',['a2:j' lim]);
                    xlswrite([fname,'.xlsx'],S11up,'Sheet1',['p2:p' lim]);
                    xlswrite([fname,'.xlsx'],S21up,'Sheet1',['l2:l' lim]);
                    xlswrite([fname,'.xlsx'],S11gd,'Sheet1',['r2:r' lim]);
                    xlswrite([fname,'.xlsx'],S21gd,'Sheet1',['n2:n' lim]);
                    % add RLGC values to individual files
                    xlswrite([fname,'.xlsx'],R,'Sheet1',['t2:t' lim]);
                    xlswrite([fname,'.xlsx'],L,'Sheet1',['u2:u' lim]);
                    xlswrite([fname,'.xlsx'],G,'Sheet1',['v2:v' lim]);
                    xlswrite([fname,'.xlsx'],Cap,'Sheet1',['w2:w' lim]);
                    xlswrite([fname,'.xlsx'],Corrected_ATT,'Sheet1',['y2:y' lim]);

                    %writting S21 amplitude to summmary table
                    xlswrite(outputfile,{fname},'AMP_S21',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,dataout(:,7),'AMP_S21',[newLabels{Y+2},'2']);
                    %writting S21 group delay to summmary table
                    xlswrite(outputfile,{fname},'GD_S21',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,S21gd,'GD_S21',[newLabels{Y+2},'2']);
                    %writting S21 unwraped phase to summmary table
                    xlswrite(outputfile,{fname},'UP_S21',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,S21up,'UP_S21',[newLabels{Y+2},'2']);
                    %writting S11 amplitude to summmary table
                    xlswrite(outputfile,{fname},'AMP_S11',[newLabels{Y+2},'1'])
                    xlswrite(outputfile,dataout(:,3),'AMP_S11',[newLabels{Y+2},'2']);
                    %writting S11 group delay to summmary table
                    xlswrite(outputfile,{fname},'GD_S11',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,S11gd,'GD_S11',[newLabels{Y+2},'2']);
                    %writting S11 unwraped phase to summmary table
                    xlswrite(outputfile,{fname},'UP_S11',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,S11up,'UP_S11',[newLabels{Y+2},'2']);
                    %writting ATT_CONST unwraped phase to summmary table
                    xlswrite(outputfile,{fname},'ATT_CONST_dB',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,total_att_constant,'ATT_CONST_dB',[newLabels{Y+2},'2']);
                    %writting Corrected atten_CONST unwraped phase to summmary table

                    xlswrite(outputfile,{fname},'CORRECTED_ATT_CONST_dB',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,Corrected_ATT,'CORRECTED_ATT_CONST_dB',[newLabels{Y+2},'2']);

                    %writting PHASE_CONST unwraped phase to summmary table
                    xlswrite(outputfile,{fname},'PHASE_CONST_dB',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,phase_constant_dB,'PHASE_CONST_dB',[newLabels{Y+2},'2']);
                    %write the RLGC values
                    xlswrite(outputfile,{fname},'R',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,R,'R',[newLabels{Y+2},'2']);
                    xlswrite(outputfile,{fname},'L',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,L,'L',[newLabels{Y+2},'2']);
                    xlswrite(outputfile,{fname},'G',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,G,'G',[newLabels{Y+2},'2']);
                    xlswrite(outputfile,{fname},'C',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,Cap,'C',[newLabels{Y+2},'2']);
                    % xlswrite(outputfile,{fname},'a_diel',[newLabels{Y+2},'1']);
                    % xlswrite(outputfile,a_diel,'a_diel',[newLabels{Y+2},'2']);
                    % xlswrite(outputfile,{fname},'corrected_metal_S21_dB',[newLabels{Y+2},'1']);
                    % xlswrite(outputfile,corrected_metal_S21_dB,'corrected_metal_S21_dB',[newLabels{Y+2},'2']);
                    % xlswrite(outputfile,{fname},'total_corrected_metal_S21',[newLabels{Y+2},'1']);
                    % xlswrite(outputfile,total_corrected_metal_S21,'total_corrected_metal_S21',[newLabels{Y+2},'2']);

                    xlswrite(outputfile,{fname},'dB_a_diel',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,dB_a_diel,'dB_a_diel',[newLabels{Y+2},'2']);

                    xlswrite(outputfile,{fname},'dB_corrected_metal_S21',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,dB_corrected_metal_S21,'dB_corrected_metal_S21',[newLabels{Y+2},'2']);

                    xlswrite(outputfile,{fname},'New_total_corrected_metal_S21',[newLabels{Y+2},'1']);
                    xlswrite(outputfile,New_total_corrected_metal_S21,'New_total_corrected_metal_S21',[newLabels{Y+2},'2']);

            Y = Y+1;

            c = k;
                if (k == 1)
                    fprintf(' %d sheet has been added to its assigned summary folder \n',round(c))
                else
                    fprintf(' %d sheets have been added to their assigned summary folder \n',round(c))
                end        
    end
    
end

datetime.setDefaultFormats('default','hh:mm a MM/dd/yyyy ')
timeStamp = datetime;
disp(' ')
fprintf('Finished at: %s \n',timeStamp)
disp(' ')