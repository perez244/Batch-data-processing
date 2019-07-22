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
%-------------------------------------------------------------------------------------------------------------------------------------------
% define folder with the raw data. Either have a dialog box pop and
% enter(paste) address, or comment out the 3 lines below and manually input
% paste the address in folder = 'inside here'
prompt = 'Enter folder address: ';
Address = inputdlg(prompt,'Folder Location',[1 70]);
folder = cell2mat(Address(1,1));
%folder = 'C:\Users\jap11\Desktop\Jesus Perez Summer 2019\Data Process Sampling';
%-------------------------------------------------------------------------------------------------------------------------------------------
% name your summary file
outputfilename = inputdlg('Name your summary file: ','Summary file name',[1 50]);
outputfilename = cell2mat(outputfilename(1,1));
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
b = str2double(b);
%b = enter # of files here;

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
fprintf(' There will be %d Summary files. \n',round(N))
disp(' ')
disp('Beginning to create the summary files...')
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

    fprintf(' %d Summary Files have been created \n',round(L))
    
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

    if(rem(k,b) == 0)% b = 5
        %check
        if(index == N)% N=3
            %do nothing
        else
            OutPutFile = char(fileNames(1,index));
            fname=fnames(k).name;
            A=importdata(fname,' ',9);
            data=A.data;
            
                    dataout(:,1)=data(:,1);
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
                    A= (((1)+ S11-S22-Delta_S)./(2.*S21));
                    B= (((1+ S11+S22+Delta_S).*z0)./(2.*S21));
                    C= (1-S11-S22+Delta_S)./(2.*S21.*z0);
                    D= (1-S11+S22-Delta_S)./(2.*S21);

                    Propagation_constant = asinh((B.*C).^(1./2))./length;
                    att_constant=real(Propagation_constant);
                    phase_constant=imag(Propagation_constant);

                    %******************************The RLGC values*****************************
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
                    %**************************************************************************

                    % look-up table for output
                    title={'frequency (Hz)' 'Omega (w)' 'S11'   'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'  ' '  'R'  'L'  'G'  'C'};

                    xlswrite([fname,'.xlsx'],title,'Sheet1','a1:w1');
                    xlswrite([fname,'.xlsx'],dataout,'Sheet1','a2:j802');
                    xlswrite([fname,'.xlsx'],S11up,'Sheet1','p2:p802');
                    xlswrite([fname,'.xlsx'],S21up,'Sheet1','l2:l802');
                    xlswrite([fname,'.xlsx'],S11gd,'Sheet1','r2:r802');
                    xlswrite([fname,'.xlsx'],S21gd,'Sheet1','n2:n802');
                    %writting S21 amplitude to summmary table
                    xlswrite(OutPutFile,{fname},'AMP_S21',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,dataout(:,5),'AMP_S21',[newLabels{Y+2},'2']);
                    %writting S21 group delay to summmary table
                    xlswrite(OutPutFile,{fname},'GD_S21',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,S21gd,'GD_S21',[newLabels{Y+2},'2']);
                    %writting S21 unwraped phase to summmary table
                    xlswrite(OutPutFile,{fname},'UP_S21',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,S21up,'UP_S21',[newLabels{Y+2},'2']);
                    %writting S11 amplitude to summmary table
                    xlswrite(OutPutFile,{fname},'AMP_S11',[newLabels{Y+2},'1'])
                    xlswrite(OutPutFile,dataout(:,3),'AMP_S11',[newLabels{Y+2},'2']);
                    %writting S11 group delay to summmary table
                    xlswrite(OutPutFile,{fname},'GD_S11',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,S11gd,'GD_S11',[newLabels{Y+2},'2']);
                    %writting S11 unwraped phase to summmary table
                    xlswrite(OutPutFile,{fname},'UP_S11',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,S11up,'UP_S11',[newLabels{Y+2},'2']);
                    %writting ATT_CONST unwraped phase to summmary table
                    xlswrite(OutPutFile,{fname},'ATT_CONST',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,att_constant,'ATT_CONST',[newLabels{Y+2},'2']);
                    %writting PHASE_CONST unwraped phase to summmary table
                    xlswrite(OutPutFile,{fname},'PHASE_CONST',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,phase_constant,'PHASE_CONST',[newLabels{Y+2},'2']);
                    %write the RLGC values
                    xlswrite(OutPutFile,{fname},'R',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,R,'R',[newLabels{Y+2},'2']);
                    xlswrite(OutPutFile,{fname},'L',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,L,'L',[newLabels{Y+2},'2']);
                    xlswrite(OutPutFile,{fname},'G',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,G,'G',[newLabels{Y+2},'2']);
                    xlswrite(OutPutFile,{fname},'C',[newLabels{Y+2},'1']);
                    xlswrite(OutPutFile,Cap,'C',[newLabels{Y+2},'2']);
            
                    c = k;
                    fprintf(' %d sheets have been added to their assigned summary folder \n',round(c))
                    index = index +1;
                    Y = 1;
        end
        
    else
        OutPutFile = char(fileNames(1,index));
        
        fname=fnames(k).name;
        A=importdata(fname,' ',9);
        data=A.data;
        
            dataout(:,1)=data(:,1);
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
            A= (((1)+ S11-S22-Delta_S)./(2.*S21));
            B= (((1+ S11+S22+Delta_S).*z0)./(2.*S21));
            C= (1-S11-S22+Delta_S)./(2.*S21.*z0);
            D= (1-S11+S22-Delta_S)./(2.*S21);

            Propagation_constant = asinh((B.*C).^(1./2))./length;
            att_constant=real(Propagation_constant);
            phase_constant=imag(Propagation_constant);

            %******************************The RLGC values*****************************
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
            %**************************************************************************

            % look-up table for output
            title={'frequency (Hz)' 'Omega (w)' 'S11'   'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'  ' '  'R'  'L'  'G'  'C'};

            xlswrite([fname,'.xlsx'],title,'Sheet1','a1:w1');
            xlswrite([fname,'.xlsx'],dataout,'Sheet1','a2:j802');
            xlswrite([fname,'.xlsx'],S11up,'Sheet1','p2:p802');
            xlswrite([fname,'.xlsx'],S21up,'Sheet1','l2:l802');
            xlswrite([fname,'.xlsx'],S11gd,'Sheet1','r2:r802');
            xlswrite([fname,'.xlsx'],S21gd,'Sheet1','n2:n802');
            %writting S21 amplitude to summmary table
            xlswrite(OutPutFile,{fname},'AMP_S21',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,dataout(:,5),'AMP_S21',[newLabels{Y+2},'2']);
            %writting S21 group delay to summmary table
            xlswrite(OutPutFile,{fname},'GD_S21',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,S21gd,'GD_S21',[newLabels{Y+2},'2']);
            %writting S21 unwraped phase to summmary table
            xlswrite(OutPutFile,{fname},'UP_S21',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,S21up,'UP_S21',[newLabels{Y+2},'2']);
            %writting S11 amplitude to summmary table
            xlswrite(OutPutFile,{fname},'AMP_S11',[newLabels{Y+2},'1'])
            xlswrite(OutPutFile,dataout(:,3),'AMP_S11',[newLabels{Y+2},'2']);
            %writting S11 group delay to summmary table
            xlswrite(OutPutFile,{fname},'GD_S11',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,S11gd,'GD_S11',[newLabels{Y+2},'2']);
            %writting S11 unwraped phase to summmary table
            xlswrite(OutPutFile,{fname},'UP_S11',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,S11up,'UP_S11',[newLabels{Y+2},'2']);
            %writting ATT_CONST unwraped phase to summmary table
            xlswrite(OutPutFile,{fname},'ATT_CONST',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,att_constant,'ATT_CONST',[newLabels{Y+2},'2']);
            %writting PHASE_CONST unwraped phase to summmary table
            xlswrite(OutPutFile,{fname},'PHASE_CONST',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,phase_constant,'PHASE_CONST',[newLabels{Y+2},'2']);
            %write the RLGC values
            xlswrite(OutPutFile,{fname},'R',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,R,'R',[newLabels{Y+2},'2']);
            xlswrite(OutPutFile,{fname},'L',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,L,'L',[newLabels{Y+2},'2']);
            xlswrite(OutPutFile,{fname},'G',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,G,'G',[newLabels{Y+2},'2']);
            xlswrite(OutPutFile,{fname},'C',[newLabels{Y+2},'1']);
            xlswrite(OutPutFile,Cap,'C',[newLabels{Y+2},'2']);

            Y = Y+1;

            c = k;
            fprintf(' %d sheets have been added to their assigned summary folder \n',round(c))

    end
    

end

datetime.setDefaultFormats('default','hh:mm a MM/dd/yyyy ')
timeStamp = datetime;
disp(' ')
fprintf('Finished at: %s \n',timeStamp)
disp(' ')
