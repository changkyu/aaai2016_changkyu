mvstr = version;
i = strfind(mvstr,'.');
mv = str2num(mvstr(1:i(1)-1))*100 + str2num(mvstr(i(1)+1:i(2)-1));

clear mexoptions libfiles
if mv < 703
    % If Matlab version is older than 7.3
    mexoptions{1} = '-DOLDMATLABAPI';
else
    mexoptions{1} = '-largeArrayDims';
end

if 0
    % The following lines are needed in Windows
    % Path to pthreads.h
    mexoptions{2} = '-I"C:\Unix\pthreads-win32"';
    % Pthreads library (e.g., pthreads for Win32)
    libfiles{1} = 'c:\Unix\pthreads-win32\libpthreadGCE2.a';
    % pthreadGCE2.dll from pthreads for Win32 is used
else
    libfiles = {};
end

mex( mexoptions{:}, 'subtract_mu.cpp', libfiles{:} )
mex( mexoptions{:}, 'errpca_pt.cpp', libfiles{:} )
mex( mexoptions{:}, 'errpca_diag.cpp', libfiles{:} )

