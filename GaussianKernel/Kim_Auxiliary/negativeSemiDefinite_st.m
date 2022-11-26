clear; 
close all 


list = importdata("eigenValuesTest.txt"); 


y=zeros(length(list)); 


scatter(list(:,1), y, 'ro'); hold on; 
scatter(list(:,2), y, 'bo')



