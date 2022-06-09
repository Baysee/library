function [hist_ch,binaxis] = coincidence_hist(TimeStamps,Channel,ch1,ch2,bin_average,windowsize_bin)
%% Analysis script for T2 mode 

arraysize = ceil(windowsize_bin/(bin_average));

binaxis = ((0:2*arraysize-1)-arraysize)*bin_average;

hist_ch = [];
hist_ch(1:2*arraysize) = 0;

for n1=1:length(TimeStamps)
    
    if Channel(n1)==ch1
        
        for n2=(n1+1):length(TimeStamps)
           
            if abs(TimeStamps(n1)-TimeStamps(n2))>=windowsize_bin,
                break
            end
            
            
            if Channel(n2)==ch2
                h = -(TimeStamps(n1)-TimeStamps(n2));
                hist_pos = ceil(h/bin_average)+arraysize;
                hist_ch(hist_pos) = hist_ch(hist_pos) + 1;
            end
            
        end
        
        for n2=(n1-1):-1:1
            
             if abs(TimeStamps(n1)-TimeStamps(n2))>=windowsize_bin,
                break
            end
            
            if Channel(n2)==ch2
                h = -(TimeStamps(n1)-TimeStamps(n2));
                hist_pos = ceil(h/bin_average)+arraysize;
                hist_ch(hist_pos) = hist_ch(hist_pos) + 1;
            end
            
        end
       
    end    
end        