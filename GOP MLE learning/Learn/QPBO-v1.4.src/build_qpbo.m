function build_qpboMex
% build_qpboMex builds package qpboMex
%
% Anton Osokin (firstname.lastname@gmail.com),  24.09.2014

codePath = './';

srcFiles = { 'qpboMex.cpp'};
allFiles = '';
for iFile = 1 : length(srcFiles)
    allFiles = [allFiles, ' ', srcFiles{iFile}];
end

cmdLine = ['mex -v',allFiles,' -output qpboMex -largeArrayDims ', '-I', codePath];
eval(cmdLine);




