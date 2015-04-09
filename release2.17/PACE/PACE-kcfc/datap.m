function [ID_subj,ID_class,data] = datap(filename,nrow);

%  Transform the input data set into one with format for kcfc.m.
% Input 
%  filename - file name of data set with the following format:
%             1st column = ID number of subject
%             2nd column = External class label
%             3rd column = Recording time
%             4th column = Oberservations measured at the recording time
%  nrow     - number of rows of the input data set in "filename"
% Output 
%  ID_subj  - n x 1 vector of ID numbers of n subjects
%  ID_class - n x 1 vector of external class labels of n subjects
%  data     - object of input data for kcfc.m, including
%             data.isobs - n x m matrix of indicators for data status.
%                          isobs(i,j)=1: the ith curve is observed at time Tin(i,j).
%                          isobs(i,j)=0: no observation for the ith curve at time Tin(i,j).
%             data.Tin   - n x m matrix of recording time points of n subjects 
%                          m is the maximun number of time points for subjects.
%             data.Yin   - n x m matrix of observation corresponding to time points Tin

  fid = fopen(filename,'rt');
  data_in = reshape(fscanf(fid,'%f'),[4,nrow]);
  data_in = data_in';
  data_subj = data_in(:,1);
  data_class = data_in(:,2);
  data_t = data_in(:,3);
  data_y = data_in(:,4);
  [ID_subj,loc_m] = unique(data_subj);
  ID_subj = data_subj(sort(loc_m));
  n = length(ID_subj);
  ID_class = zeros(n,1);
  for i = 1:n;
      tempid = find(data_subj == ID_subj(i));
      nobs(i) = length(tempid);
      ID_class(i) = data_class(tempid(1));
  end;
  m = max(nobs);
  Tin = zeros(n,m); Yin = Tin; isobs = Tin;
  for i = 1:n;      
      tempid = find(data_subj == ID_subj(i));
      nobs(i) = length(tempid);
      Tin(i,[1:nobs(i)]) = data_t(tempid)'; 
      Yin(i,[1:nobs(i)]) = data_y(tempid)';
      isobs(i,[1:nobs(i)]) = 1;
  end;
  data.isobs = isobs;
  data.Tin = Tin;
  data.Yin = Yin;