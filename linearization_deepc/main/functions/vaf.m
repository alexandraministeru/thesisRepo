function v = vaf(y,yhat)
v = max(0, (1 - (norm(y-yhat)^2)/(norm(y))^2)*100);
end

