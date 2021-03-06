%This function performs functional evaluation for the existing subjects
% based on the FPCA() function call.
% ypred=FPCAeval(yy, subjectID, newt)
%=======
% Input:
%=======
%  yy: an aggregated object, returned from FPCA().
%
%  subjectID: a 1*m (m <= n) row vector, denotes the subject ID
%             to be predicted.
%             if set to [], then subjectID = 1:n
%
%  newt:  i) a 1*mm row vector of new time points for all subjects
%            defined in subject ID,
%           e.g., if subjectID = 1:3, request return of the predicted 
%           measurements for the first three subjects, evaluated
%           at the same time points at newt.
%         ii) a 1*mm cell array of new time points, where newt{1}
%             denotes the time points for the subject whose ID is
%             subjectID(1). This allows different subjects to be
%             evaluated at different new time points. 
%========
% Output:
%======== 
%   
%    ypred: 1*m cell array of the predicted measurements corresponding
%          to new time points defined in newt.
%
%  example:
%
%  yy = FPCA(y,t,p);
%  subjectID = [1 2 6];
%  newt = 0.1:0.1:0.5;
%  ypred = FPCAeval(yy,subjectID, newt);  %evaluate subect 1,2,6 at
%                                         %the same time points
%  or
%  subjectID = [1 2 6];
%  newt = {[0.1 0.2], [0.3 0.4 0.5], [0.6 0.8 1]};
%  ypred = FPCAeval(yy,subjectID,newt)    %evaluate subject 1,2,6 at
%                                         %different time points
%                                          
% See also FPCA, PCA
%   
function ypred=FPCAeval(yy, subjectID, newt)

  if isempty(subjectID)
     subjectID = 1:length(getVal(yy,'y'));
  end
  if length(newt) == 1
      newt = newt{1};
  end 
  if iscell(newt) == 1
     sameTIME = 0;
     if length(subjectID) ~= length(newt)
         fprintf(1,'Error: length of the subjectID and length of the newt must be the same!\n');
         ypred = [];
         return;
     end
  elseif isnumeric(newt) == 1
     sameTIME = 1;
     if size(newt,1) > 1
        newt = newt';
     end
  else
     fprintf(1,'Error: newt must be either a cell array or a numeric vector!\n');
  end

  out1 = getVal(yy,'out1copy');
  mu = getVal(yy,'mucopy');
  phi = getVal(yy,'phicopy');
  xi_est = getVal(yy,'xi_est');
  no_opt = getVal(yy,'no_opt');

  if sameTIME
    newmu = interp1(out1,mu, newt,'spline');
    newphi = interp1(out1, phi,newt, 'spline');
    if no_opt == 1
       if size(newphi,2) > 1
         newphi = newphi'; 
       end
    end
    ypred = repmat(newmu,length(subjectID),1)+xi_est(subjectID,:)*newphi';
    ypred = num2cell(ypred,2);
  else
    ypred = cell(1,length(subjectID));
    for i=1:length(subjectID)
        newmu = interp1(out1,mu,newt{i},'spline');
        newphi =interp1(out1,phi,newt{i},'spline');
        if no_opt == 1
           if size(newphi,2) > 1
              newphi = newphi';
           end 
        end
        ypred{i} = newmu + xi_est(subjectID(i),:)*newphi';
    end
  end

end

