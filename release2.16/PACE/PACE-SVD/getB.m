% function [b] = getB(pc_x, yy, y,t_y, isYFun, K_x, K_y,method)
%=======
%Input:
%=======
%  pc_x:  The first K_x FPC scores of X
%    yy:  The returned values from FPCA. See FPCA() for more details.
% y    :  The original response variable, used when the response is scalar        
% t_y  :  The original input observation time for response
% isYFun:  a logical value, isYFun = 1 (Y is a functional response)
%                           isYFun = 0 (Y is a scalar response)
% family:  a character string naming the distribution of the response Y.
%          Accepted families include: 'normal'(default), 'binomial', 'poisson'.
%          used only when the response Y is a scalar.
% link:    link function used in the model. Valid options:
%              family          link
%             'normal'        'identity'(default)  
%             'binomial'      'logit'(default),'probit','comploglog'
%             'poisson'       'log'(default)
% K_x:   positive integer; number of principal components of predictor x used in regression
%        must be smaller than or equal to the maximum number of principal components give by FPAC
%        default is the number of pc selected by FPCA
% K_y:   positive integer; number of principal components of response y used in regression
%        must be smaller than or equal to the maximum number of principal components give by FPAC
%        default is the number of pc selected by FPCA
% method:  a string taking the following two values:
%            'YPC': fitting the functional linear regression;
%            'YCO': predicting with a functional linear regression. 
%          Refer to the input FIT of FPCreg for more details.
%=======
%Output:
%=======
%  b:   estimated coefficients of beta in the basis represatation 
%       a K_x*K_y matrix if Y is functional; a 1*K_x vector is Y is scalar 

function [b] = getB(pc_x, pc_y,K_x, K_y, isYFun)
        
    if isYFun == 1
           
            pc_y = pc_y(:,1:K_y);     % truncated PC for y
            b = pinv(pc_x' * pc_x) * pc_x' * pc_y;
            
    else
        
        if strcmp(family,'normal') == 1 && strcmp(link,'identity') == 1
            pc_y = y-mean(y);
            b = zeros(1,K_x);
            for k = 1:K_x
                b(k) = pc_y*pc_x(:,k)/(pc_x(:,k)'*pc_x(:,k));
            end
        else
            b = glmfit(pc_x,y',family,'link',link)';
        end

    end

