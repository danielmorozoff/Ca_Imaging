% STANDARDTEMPORALICA Performs spatial ica.
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

function [outputtrace,outputimage] = standardica(volumestack)
    
    % make input for ICA into a time x space matrix
    datasize = size(volumestack);    
    X = reshape(double(volumestack),(datasize(1)*datasize(2)),datasize(3));
    
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
    % compute temporal signals
    ica_sig = ica_W * X;      

    % find strongest temporal compoenent
    [~, maxloc] = max(sum(ica_sig,2));    
    
    icatrace = ica_sig(maxloc,:);
    icatrace = (icatrace-median(icatrace))/median(icatrace);
    icaimage = reshape(ica_W(maxloc,:),datasize(1),datasize(2));
    minimg = min(icaimage(:));
    maximg = max(icaimage(:));
    normicaimage = (icaimage-minimg)/(maximg - minimg);

    % segment as necessary
    bwicaimage = im2bw(normicaimage,graythresh(normicaimage));
%    bwicaimage = im2bw(icaimage,.5);
    bwicaimage = bwmorph(bwicaimage,'open');
    bwicaimage = bwmorph(bwicaimage,'open');

    cellmap = bwlabel(bwicaimage);
    numcells = max(cellmap(:));
    outputimage = zeros(datasize(1),datasize(2),numcells);
    outputtrace = zeros(datasize(3),numcells);
    for i=1:numcells
        outputimage(:,:,i) = icaimage .* (cellmap==i);
        outputtrace(:,i) = X' * reshape(outputimage(:,:,i),[datasize(1)*datasize(2),1]);
    end
    keyboard
end