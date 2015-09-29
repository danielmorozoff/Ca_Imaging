% SPICA Performs spatial ica.
%
% standardtemporalica receivea a volumestack of image data and performs
% temporal ICA. It returns the largest component
% X is an X x Y x T volume stack
%
% type: function
%
% inputs:
%   volumestack:  3-D array of imaging data X x Y x T
%    
% outputs:
%   outputtrace:  1-D array containing trace
%   outputimage:  2-D image
%
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 09:31pm

function [ica_A,ica_S] = spica(volumestack)
    
    % make input for ICA into a time x space matrix
    datasize = size(volumestack);    
    X = reshape(volumestack,(datasize(1)*datasize(2)),datasize(3))';
    
    ica_A_guess = randn(size(X));
    termtol = 1e-6;
    maxrounds = 1000;

    numSamples = size(X,2);

    B = ica_A_guess;
    BOld = zeros(size(B));

    iternum = 0;
    minAbsCos = 0;

    errvec = zeros(maxrounds,1);
    try
        while (iternum < maxrounds) && ((1 - minAbsCos)>termtol)
            iternum = iternum + 1;
            % Symmetric orthogonalization.
            B = (X * ((X' * B) .^ 2)) / numSamples;

            % eran's simple inverse method
            %        B = B * real(inv(B' * B)^(1/2));
            % a numerically more accurate inverse
            B = B * real( ((B'*B)\speye(size(B'*B)))^(1/2));
            % Test for termination condition.
            minAbsCos = min(abs(diag(B' * BOld)));

            BOld = B;
            errvec(iternum) = (1 - minAbsCos);
        end
    
        if iternum<maxrounds
            fprintf('Convergence in %d rounds.\n', iternum)
        else
            fprintf('Failed to converge; terminating after %d rounds, current change in estimate %3.3g.\n', ...
                iternum, 1-minAbsCos)
        end
    catch
        outputtrace = [];
        outputimage = [];
        return
    end
        
    % output of ICA
    ica_A = B;
    ica_W = ica_A';
    ica_S = ica_W * X;      
    
end