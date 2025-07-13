function A_inv = computeTruncatedSVDInverse(A, threshold)
    [U, S, V] = svd(A);
    above_indices = diag(S) > threshold;
    S_inv = diag(1 ./ diag(S(above_indices, above_indices)));
    A_inv = V(:, above_indices) * S_inv * U(:, above_indices)';
end