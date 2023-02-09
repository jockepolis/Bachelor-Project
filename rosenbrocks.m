q = integral(rosenbrocks,-5,5)

function y = rosenbrocks(x)
    y = sum(100*(x(:,1:end-1).^2 - x(:,2:end)).^2 + (x(:,1:end-1)-1).^2);
end