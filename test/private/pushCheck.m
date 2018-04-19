function pushCheck()
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));

paths=genpath(ravenPath);

for p=regexp(paths,':','split')
    try
        tmp=regexp(ls([p{1} '/*.m']),'[\n\t]','split');
        tmp=tmp(~cellfun('isempty',tmp));
        %fn{end+1}=tmp;
        for file=tmp
            disp(file)
            checkcode(file)
            disp('----------------------------------')
        end
    catch Ex
    end
end
end
