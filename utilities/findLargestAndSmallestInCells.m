function [smallest,largest] = findLargestAndSmallestInCells(data)

    N = length(data);
    smallest = min(data{1});
    largest = max(data{1});
    for i=2:N
        a = min(data{i});
        b = max(data{i});
        if a < smallest
            smallest = a;
        end
        if b > largest
            largest = b;
        end
    end
        