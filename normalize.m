function [measure,norm_data] = normalize(dataset, cell_line, treatment, hp, hn)
    header = dataset.colheaders;
    data_matrix = dataset.data;
    celline_index = strmatch(['TR:' cell_line],header);
    cellline_filter = find(data_matrix(:,celline_index)==1);
    data_matrix = data_matrix(cellline_filter,:);
    
    control_index = strmatch('TR:DMSO', header);
    control_filter = find(data_matrix(:,control_index)==1);
    data_control = data_matrix(control_filter,:);
    
    treatment_index = strmatch(['TR:' treatment],header);
    treatment_filter = find(data_matrix(:,treatment_index)==1);
    data_treatment = data_matrix(treatment_filter,:);
    
    da_index = strmatch('DA:',header);
    
    measure={};
    norm_data = [];
    for i=1:length(da_index)  % extract median from raw data
        da = header{da_index(i)};
        node_name = da(4:length(da));
       
        dv_index = strmatch(['DV:' node_name],header);
        if isempty(dv_index)
            disp(strcat('error: DV is not found for ',node_name));
        else
            control = data_control(:,dv_index);
            treatment = data_treatment(:,dv_index);
            group = unique(data_control(:,da_index(i)));
            group = group(~isnan(group));
            norm_tmp = [];
            
            da_value = data_control(:,da_index(i));
            
            
            for j=1:length(group)
                ctrl_tmp = median(control(da_value==group(j)));
                trt_tmp = median(treatment(da_value==group(j)));
                %norm_level = (trt_tmp/ctrl_tmp)^hn/(hp^hn+(trt_tmp/ctrl_tmp)^hn);
                norm_tmp = [norm_tmp; trt_tmp];
            end
            if length(norm_tmp)>7
                if isempty(measure)
                    measure = {node_name};
                else
                    measure=[measure node_name];
                end
                norm_data = [norm_data norm_tmp];
            end
        end
    end
    
    % compute z-score with respect to phosphorylation level at time 0
    for i = 1:size(norm_data,2)        
        norm_data(:,i) = (norm_data(:,i)-norm_data(1,i))/std(norm_data(:,i));
    end
    
    for i = 2:size(norm_data,1)
        norm_data(i,:) = (norm_data(i,:)-min(norm_data(i,:)))/(max(norm_data(i,:))-min(norm_data(i,:)));
    end
end













