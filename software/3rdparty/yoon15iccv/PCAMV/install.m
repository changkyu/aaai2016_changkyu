mvstr = version;
i = strfind(mvstr,'.');
mv = str2num(mvstr(1:i(1)-1))*100 + str2num(mvstr(i(1)+1:i(2)-1));

clear mexoptions
if mv < 703
    % If Matlab version is older than 7.3
    mexoptions{1} = '-DOLDMATLABAPI';
else
    mexoptions{1} = '-largeArrayDims';
end

mexoptions{2} = '-DNOTHREADS';

mex( mexoptions{:}, 'subtract_mu.cpp' )
mex( mexoptions{:}, 'errpca_pt.cpp' )
mex( mexoptions{:}, 'errpca_diag.cpp' )

