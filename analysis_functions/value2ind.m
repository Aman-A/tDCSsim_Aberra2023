function val_ind = value2ind(val,vec)
%VALUE2IND Takes an input value val and vector vec, finds closest value in vector to
%input val, returns the index of that value 
%picks out first column if time input is actually a matrix, accomodates
%for input of 'raw' array from event structure
%
% AUTHOR : Aman Aberra
if isrow(vec)
   vec = vec'; 
end
vec = vec(:,1);
diff_vec = abs(vec-val);
val_ind = find(diff_vec == min(diff_vec));
if length(val_ind)> 1
    val_ind = val_ind(2); %ensures only one value is chosen, if in between, round up            
end
end