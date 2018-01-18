function delta_theta=M_step(adjac_mat,q_of_h,input_tree,theta,hidden_indexes,output_indexes,input_indexes,training_data, lambda)


delta_theta=zeros(size(theta));
temp_theta=theta;
for i=1:size(training_data,1)
    temp_example=training_data(i,:);
    temp_example1=training_data(i,:);
    % calculate first part of the drivetives
    for n=1:size(q_of_h,1)

        temp_example(1,hidden_indexes)=q_of_h(n,1);
        for j=1:size(theta,1)
            for k=j+1:size(theta,1)
                if adjac_mat(k,j)>0
                    delta_theta(j,k,temp_example(j),temp_example(k))= ...
                        delta_theta(j,k,temp_example(j),temp_example(k))+q_of_h(n,i+size(hidden_indexes,2));
                    delta_theta(k,j,temp_example(k),temp_example(j))=delta_theta(j,k,temp_example(j),temp_example(k));
                end
            end
        end
    end
        % calculate second part of the drivetives
        for j=1:size(theta,1)
            for k=j+1:size(theta,1)
                if adjac_mat(k,j)>0
                    dif_of_pobserved=0;
                    temp_example1(1,hidden_indexes)=ones(size(hidden_indexes))*-1;
                    temp_example1(1,output_indexes)=ones(size(output_indexes))*-1;
                    p_of_observed=inference(input_tree,theta,temp_example1);
                    for l=1:2
                            for m=1:2
                                if (ismember(j,input_indexes)&&temp_example1(1,j)~=l)...
                                        || (ismember(k,input_indexes)&& temp_example1(1,k)~=m)
                                        dif_of_pobserved=0;
                                else

                                    temp_example1(1,j)=l;
                                    temp_example1(1,k)=m;
                                    dif_of_pobserved=inference(input_tree,theta,temp_example1);
                                    delta_theta(j,k,l,m)= ...
                                    delta_theta(j,k,l,m)-(dif_of_pobserved/p_of_observed);
                                    delta_theta(k,j,m,l)=delta_theta(j,k,l,m);
                                end
                            end

                    end


                end
            end
        end



end
delta_theta =delta_theta - lambda * theta;
end