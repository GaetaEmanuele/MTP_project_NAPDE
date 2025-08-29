function [index] = check_dependece(tent1,tent2)
index = [];
t1 = tent1(:,2);
t2 = tent2(:,2);
tol = 1e-4;
for i=1:length(t1)
    for j=1:length(t2)
        if (abs(t1(i)-t2(j))<tol) && (t1(i)>0)
            index = [index , i];
        end
    end
end


end