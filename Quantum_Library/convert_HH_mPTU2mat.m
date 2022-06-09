function convert_HH_mPTU2mat(filename,mode)

[folder,filen,ext] = fileparts(filename);
% for matlab written HH data files

if nargin <2
    mode = 3;
end

max_file_size = 2*1024;

fid=fopen(filename);

%% constants
chunk_size = 1000000;

%% read data

data = [];
save_count = 0;
cnt_ph = 0;
cnt_ov = 0;
cnt_ma = 0;

% choose right decode function
if mode==3, %T3 mode

        %---------------T3----------------------------------------
        
        OverflowCorrection = 0;
        T3WRAPAROUND = 1024;
      
        while 1            
            
            data_h = zeros(chunk_size,3);
            data_h_count = 1;
            
            T3Record = fread(fid, chunk_size, 'ubit32');     % all 32 bits:
            nsync = bitand(T3Record,1023);       % the lowest 10 bits:
            dtime = bitand(bitshift(T3Record,-10),32767);   % the next 15 bits:
            channel = bitand(bitshift(T3Record,-25),63);   % the next 6 bits:
            special = bitand(bitshift(T3Record,-31),1);   % the last bit:
            
            for recnum=1:length(T3Record)
                
                if special(recnum) == 0   % this means a regular input channel
                    true_nSync = OverflowCorrection + nsync(recnum);
                    %  one nsync time unit equals to "syncperiod" which can be
                    %  calculated from "SyncRate"
                    data_h(data_h_count,:)=[channel(recnum)+1,true_nSync,dtime(recnum)];
                    data_h_count=data_h_count+1;
                    %cnt_ph = cnt_ph + 1;
                else    % this means we have a special record
                    if channel(recnum) == 63  % overflow of nsync occured
                        %  nsync indicates the number of overflows - THIS IS NEW IN FORMAT V2.0
                        OverflowCorrection = OverflowCorrection + T3WRAPAROUND * nsync(recnum);
                    %    cnt_ov = cnt_ov + 1;
                    end;
                end;
            end;
            
            if length(T3Record)~=chunk_size,
                data = [data;data_h(1:data_h_count,:)];
                break,
            else
                data = [data;data_h(1:data_h_count,:)];
            end
            
            s = whos('data');
            if s.bytes/(1024^2)>max_file_size
                save([filen '_' num2str(save_count) '.mat'],'data','-v7.3')
                data =[];
                save_count = save_count+1;
            end
                
                
                
        end;
        
        save([filen '_' num2str(save_count) '.mat'],'data','-v7.3')
        
elseif mode==2,
    %---------------T2----------------------------------------
  %  disp('Not supportet, yet')
    OverflowCorrection = 0;
    %T2WRAPAROUND_V1=33552000;
    T2WRAPAROUND_V2=33554432; % = 2^25  IMPORTANT! THIS IS NEW IN FORMAT V2.0
    
    while 1
        data_h = zeros(chunk_size,2);
        data_h_count = 1;
        
        
        T2Record = fread(fid, chunk_size, 'ubit32');     % all 32 bits:
        
        dtime = bitand(T2Record,33554431);   % the last 25 bits:
        
        channel = bitand(bitshift(T2Record,-25),63);   % the next 6 bits:
        
        special = bitand(bitshift(T2Record,-31),1);   % the last bit:
        
        %             % the resolution in T2 mode is 1 ps  - IMPORTANT! THIS IS NEW IN FORMAT V2.0
        for recnum=1:length(T2Record)
            timetag = OverflowCorrection + dtime(recnum);
            if special(recnum) == 0   % this means a regular photon record
%                 cnt_ph = cnt_ph + 1;
                data_h(data_h_count,:)=[channel(recnum)+1,timetag];
                data_h_count=data_h_count+1;
            else    % this means we have a special record
                if channel(recnum) == 63  % overflow of dtime occured
                    OverflowCorrection = OverflowCorrection + T2WRAPAROUND_V2 * dtime(recnum);
%                     cnt_ov = cnt_ov + dtime;
                elseif channel(recnum) == 0  % Sync event
                    data_h(data_h_count,:)=[channel(recnum),timetag];
                    data_h_count=data_h_count+1;
                    
                end;
            end;
        end;
        
         if length(T2Record)~=chunk_size,
                data = [data;data_h(1:data_h_count,:)];
                break,
            else
                data = [data;data_h(1:data_h_count,:)];
            end
            
            s = whos('data');
            if s.bytes/(1024^2)>max_file_size
                save([filen '_' num2str(save_count) '.mat'],'data','-v7.3')
                data =[];
                save_count = save_count+1;
            end
                
             
                
        end;
        
        save([filen '_' num2str(save_count) '.mat'],'data','-v7.3')  
        
    end;

fclose(fid);
fprintf(1,'Ready!  \n\n');
%fprintf(1,'\nStatistics obtained from the data:\n');
%fprintf(1,'\n%i photons, %i overflows, %i markers.',cnt_ph, cnt_ov, cnt_ma);
fprintf(1,'\n');
end

