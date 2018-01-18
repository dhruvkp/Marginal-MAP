temp1 = theta_c_learned_qpbo(:,:,1,1);
temp2 = theta_c_learned_qpbo(:,:,1,2);
temp3 = theta_c_learned_qpbo(:,:,2,2);%...
%     + theta_c_learned_qpbo(:,:,1,2)...
%     + theta_c_learned_qpbo(:,:,2,1)...
%     +theta_c_learned_qpbo(:,:,2,2);
% % temp = temp ./4;
theta1 = temp1(4:787, 1:5);
theta2 = temp2(4:787, 1: 3);
theta3 = temp3(4:787, 1: 3);
training = images(:,1)-1;
training = training(4:787,:);
h = 1./(1+ exp(-1*training'*theta1))+1./(1+ exp(-1*training'*theta2))+1./(1+ exp(-1*training'*theta3));
%h = h - mean(h);
% h(h< mean(h))=0;
% h(h>= mean(h))=1;
h = h/3;
image = 1./(1+exp(-1*theta1*round(h)'))+1./(1+exp(-1*theta2*round(h)'))+1./(1+exp(-1*theta3*round(h)'));
image(image <  mean(image))= 0;
image(image >=  mean(image))= 1;
image = vec2mat(image,28);
image1 = vec2mat(training,28);
i = mat2gray(image);
i1 = mat2gray(image1');
i2 = mat2gray(image);
imshow(i), figure, 
imshow(i1)
% , figure, 
%imshow(i2)
