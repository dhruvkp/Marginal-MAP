addpath('../../Learn/v4')


images = loadMNISTImages('train-images-idx3-ubyte');
% labels = loadMNISTLabels('train-labels-idx1-ubyte');
 images(images>0)=1;
 images=images+1;
  h=3;

 images = [repmat(ones(1,60000)*-1,h,1);images];
 hidden = [1:h];
output = [2+h:784+h];
input=[1+h];
x_domain = [1 2];
A = zeros(784+h, 784+h);
h_con = zeros(1,784+h);
h_con(h+1:784+h) = ones(1, 784);
o_con = zeros(1, 784+h);
o_con(1:h) = ones(1,h);
A(1:h,:) = repmat(h_con, h, 1);
A(h+1: 784+h,:) = repmat(o_con, 784,1);
theta_c_learned_qpbo=learn(A,input,output,hidden,x_domain,[images(:,2)';images(:,1)'],0.1,'qpbo');

% We are using display_network from the autoencoder code
%display_network(images(:,1:100)); % Show the first 100 images
% disp(labels(1:10));