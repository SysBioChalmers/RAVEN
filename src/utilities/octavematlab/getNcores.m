function ncores=getNcores()
    % getNcores
    %   Get the number of cores available for computing
    %
    %   Usage: ncores=getNcores()

    if isOctave
        ncores=nproc();
    else
        ncores = evalc('feature(''numcores'')');
        ncores = strsplit(ncores, 'MATLAB was assigned: ');
        ncores = regexp(ncores{2},'^\d*','match');
        ncores = ncores{1};
    end
end
    