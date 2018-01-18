function [x  ] = get_sentence( dictionary,sample,output )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x=blanks(length(output));
for t=1:length(output)
    x(t)=dictionary{sample(output(t))};
end
x=string(x);

end

