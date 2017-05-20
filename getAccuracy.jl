#Hao 10/20/2015

function mkmat_incidence_factor(b)
    factor=unique(b)
    coMat= spzeros(length(b),length(factor))

    dictFactor = Dict()
    index=1
    for i in factor
        dictFactor[i]=index
        index+=1
    end

    nrow=1
    for i in b
        myindex=dictFactor[i]
        coMat[nrow,myindex]=1
        nrow=nrow+1
    end
    return full(coMat),factor
end


function get_correlation(file) #ID, fiexed, phenotype, ebv
    d=readdlm(file);
    fixed=int(d[:,2]);
    X,effects=mkmat_incidence_factor(fixed);
    X=[ones(size(X,1)) X];
    y=[d[:,3] d[:,4]]
    betaHat=pinv(X'X)*X'y;
    res=y-X*betaHat
    varcov=res'res/(size(res,1)-length(effects))
    cor=varcov[1,2]/sqrt(varcov[1,1]*varcov[2,2])
    return cor,betaHat
end
