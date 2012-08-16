function [ alpha, bsvm, Sigma, posw ] = mwu_native( K, y, C, options, verbose)

    libname = 'libmwu';
    if ~libisloaded(libname)
        place = '../mwu_native_cpp/';
        libfile = [place libname '.so'];
        hdrfile = [place 'mwu_main.h'];  
        
        loadlibrary(libfile, hdrfile);
        assert(libisloaded(libname));
    end
    
    m = size(K.data,3);
    n = length(y);

    if C > 1e6
      normtype = 0;
    else
      normtype = 2;
    end
    
    [ Sigma alpha bsvm posw ] = calllib(libname, 'run_mwu_cpp',...
        zeros(m,1), zeros(n,1), 0, zeros(n,1),...
        K.data, y, n, m, options.mwuEPS, options.mwu_iter_mult, options.mwu_cutoff,...
        C, normtype, verbose);

    posw = (posw == 1);
        
    alpha(~posw) = 0;
    posw = find(posw);
end

