function phit = loadPHIT_v1(fName)

[a,b,c] = xlsread(fName);
head = b(1,:);
iPID = strcmp(head, 'PID');
iAge = strcmp(head, 'Age');
iSex = strcmp(head, 'Sex');
iEmail = strcmp(head, 'Phishing_ID');
iClick = strcmp(head, 'Click_Status');
b = b(2:end,:);
c = c(2:end,:);
%% process PHIT
sID = unique(a(:,iPID));
for i = 1:length(sID)
    ind = find(a(:,iPID) == sID(i));
    
    phit(i).sex = a(ind(1),iSex);
    phit(i).age = a(ind(1),iAge);
    phit(i).emails = c(ind,iEmail);
    phit(i).click = a(ind,iClick);
end
