function mat = constructHankelMat(data,i,s,N)
    mat = hankel(data(i:i+s-1),data(i+s-1:i+N+s-2));
end