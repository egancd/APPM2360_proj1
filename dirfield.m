    % To run the script not as a function, uncomment
    % lines 4-7 and comment lines 9 and 53

%f = @(t,y) y;
%IndVarVals = -10:1:10;
%DepVarVals = -10:1:10;
%SegmentLength = 0.35;

function dirfield(f,IndVarVals,DepVarVals,SegmentLength)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot direction field for the first order ODE y' = f(t,y)
%   IndVarVals - vector of grid points of the independent variable
%   DepVarVals - vector of grid points of the dependent variable
%   f is an anonymous (@) function or name of an m-file in single quotes
%   SegmentLength - length of tick marks
%
% Usage example: y' = -y^2 + t
%   show direction field for t in [-1,3], y in [-2,2]
%   using a spacing of .2 for both t and y:
%
%   f = @(t,y) -y^2+t
%   dirfield(f, -1:.2:3, -2:.2:2, 0.35)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set up the grid on which the direction field is plotted
    % and compute the spacing between grid points.

[tg,yg] = meshgrid(IndVarVals,DepVarVals);
DeltaInd = IndVarVals(2) - IndVarVals(1); 
DeltaDep = DepVarVals(2) - DepVarVals(1);

    % Fix up the function if necessary.

if isa(f,'function_handle')
    f = eval(vectorize(f));
end

    % Evaluate the function on the grid and compute the scaling factor, s,
    % which helps keep the line segments' lengths relatively constant, and
    % create the direction field on the grid.

z = feval(f,tg,yg); 
s = 1./max(1/DeltaInd,abs(z)./DeltaDep)*SegmentLength;
quiver(tg,yg,s,s.*z,0,'.r'); hold on;
quiver(tg,yg,-s,-s.*z,0,'.r');

    % Fix up the axes.

axis([IndVarVals(1)-DeltaInd/2,IndVarVals(end)+DeltaInd/2, ...
      DepVarVals(1)-DeltaDep/2,DepVarVals(end)+DeltaDep/2])

end