function class = knn(test, train_data, train_labels, k)
    distances = pdist2(train_data, test);
    
    [~, ind] = sort(distances);
    
    class = mode(train_labels(ind(1:k)));
end