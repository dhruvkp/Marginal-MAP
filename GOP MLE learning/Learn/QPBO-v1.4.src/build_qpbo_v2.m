function build_qpboMex
% build_qpboMex builds package qpboMex
%
% Anton Osokin (firstname.lastname@gmail.com),  24.09.2014

codePath = './';

srcFiles = { 'QPBO_mex.cpp'};
allFiles = '';
for iFile = 1 : length(srcFiles)
    allFiles = [allFiles, ' ', srcFiles{iFile}];
end

cmdLine = ['mex -DCOST_TYPE=double CXXFLAGS="$CXXFLAGS -std=c++11" -largeArrayDims -v',allFiles,' -output QPBO_double_mex.mexw64 ', '-I', codePath];
eval(cmdLine);




