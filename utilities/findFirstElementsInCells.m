function data = findFirstElementsInCells(x)

     s = size(x);
     data = zeros(s);
     for i=1:s(1)
         for j=1:s(2)
             if ~isempty(x{i,j})
                 data(i,j) = x{i,j}(1);
             end
         end
     end