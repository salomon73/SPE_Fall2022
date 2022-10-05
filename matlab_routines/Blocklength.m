function blocklengths = Blocklength(Array)
%Blocklength return the length of non-zero blocks in number of elements of
%array Array
lengthindex=1;
blocklengths(lengthindex)=0;
for i=1:length(Array)
    if(Array(i)~=0)
        blocklengths(lengthindex)=blocklengths(lengthindex)+1;
        try
            if(Array(i+1)==0)
                lengthindex=lengthindex+1;
                blocklengths=[blocklengths,0];
            end
        catch
            return
        end
    end
end

