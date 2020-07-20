function [head, data] = readFileSubject_v2(fname, ageGroup)

% from Ziad 4/4/18:
% userId - subject response (1 - safe with high confidence, 2 - safe with low confidence, 3 - scam with low confidence, 4 - scam with high confidence)
% reactTime - reaction time in seconds 
% category - PHIT Email Category (and my custom categories for pooled scam/safe emails)
% type - weapon of influence (for pooled emails no weapon of influence - it will just record ?Type?)
% realID - real email identifier (scam or safe)
% emailCode - unique ID of each email - used to locate specific emails within excel files


fid = fopen(fname);

% a hack for dealing with switch of format in old
% userId	reactTime	category	type	hasAtt	realId	emailCode	score
hdr = textscan(fid, '%s%s%s%s%s%s%s%s', 1, 'delimiter', ',');
if strcmp(hdr{end}, 'score') ~= 1
    % if old format start again with 7 columns of data
    fclose(fid);
    fid = fopen(fname);
    hdr = textscan(fid, '%s%s%s%s%s%s%s', 1, 'delimiter', ',');
    dat = textscan(fid, '%f%f%s%s%f%s%s', 'delimiter', ',');
else
    % read 8 columns of data
    dat = textscan(fid, '%f%f%s%s%f%s%s%s', 'delimiter', ',');
end
% dat = textscan(fid, '%s%s%s%s%s%s%s', 'delimiter', ',');

for i = 1:length(hdr)
    head{i} = hdr{i}{1};
end


data.fname          = fname;
data.ageGroup       = ageGroup;
data.age            = str2num(fname(end-7:end-6));
data.sex            = fname(end-4:end-4);
data.userId         = dat{1};
data.RT             = dat{2};
data.category       = dat{3};
data.type           = dat{4};
data.hasAtt         = dat{5};
data.realId         = dat{6};
data.emailCode      = dat{7};
if isempty(data.age)
    data.age = nan;
end
for i = 1:length(data.emailCode)
    ind = data.emailCode{i} == ' ';
    data.emailCode{i}(ind) = [];
end


fclose(fid);

% userId	reactTime	category	type	hasAtt	realId	emailCode
% 2	15.72611658	Security	Liking	0	Scam	002363c