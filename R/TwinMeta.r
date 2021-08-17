### Format a number with comma delimited thousands
.s = function(x){
    formatC(x = x, 
            digits = ceiling(log10(max(x)+1)),
            big.mark = ',',
            big.interval = 3)
}

# "pmax" function that can handle empty array
.my.pmax = function(x, val){
    if( length(x) == 0 ){
        return(x)
    } else {
        return(pmax.int(x, val));
    }    
}

# Avoid creating zero p-values
.pv.nz = function(x){ .my.pmax(x, .Machine$double.xmin) }

# Orthogonalize and normalize covariates
.orthonormalizeCovariates = function(cvrt, gene){
    # Combine and add constant element
    if( NROW(cvrt) > 0 ){
        cvrtM = rbind(matrix(1,1,ncol(gene)), cvrt);
    } else {
        cvrtM = matrix(1,1,ncol(gene));
    }
    
    # Standardize and orthogonolize via QR decomposition
    q = qr(t(cvrtM));
    if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * ncol(gene) ){
        stop("Colinear or zero covariates detected");
    }
    zcvrt = t( qr.Q(q) );
    rm(q);
    
    return(zcvrt);
}

# Residualize the gene expression or genotype matrix
.ResidualizeSlice = function(slice, zcvrt){
    rowsq1 = rowSums(slice^2);
    slice = slice - tcrossprod(slice, zcvrt) %*% zcvrt;
    rowsq2 = rowSums(slice^2);
    
    # kill rows colinear with the covariates
    delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps);
    if(any(delete.rows)){
        slice[delete.rows,] = 0;
        rowsq2[delete.rows] = 1;
    }
    slice = slice / sqrt(rowsq2);
    
    rm(rowsq1, rowsq2, delete.rows);
    return(slice);
}

EstimateACE_SqD = function(gene, cvrt, twininfo){
    
    # Checks
    {
        # gene
        stopifnot( is.matrix(gene) );
        stopifnot( is.numeric(gene) );
        stopifnot( all(is.finite(gene)) );
        stopifnot( nrow(gene) > 0 );
        stopifnot( ncol(gene) > 0 );
        
        # cvrt
        if( NROW(cvrt) > 0 ){
            stopifnot( is.matrix(cvrt) );
            stopifnot( is.numeric(cvrt) );
            stopifnot( all(is.finite(cvrt)) );
            
            stopifnot( ncol(gene) == ncol(cvrt) );
        }
    }
    
    # Process twininfo
    {
        message("Parsing twininfo parameter") 
        ti = .ProcessTwininfo(twininfo = twininfo, colnames(gene));
    } # ti with (idsMZ1, idsMZ2, idsDZ1, idsDZ2, idsS, Nm, Nd, N)

    # Residualize and standardize gene expression data
    {
        message("Orthonormalizing covariates");
        zcvrt = .orthonormalizeCovariates(cvrt, gene);
        
        message("Residualizing gene expression data");
        gene = .ResidualizeSlice(gene, zcvrt);
        
        rm(zcvrt);
    } # gene, -zcvrt
    
    message("Estimating ACE model for each gene") 
    
    # Get mean squares
    {
        Rm = rowMeans((gene[,ti$idsMZ1, drop = FALSE] - gene[,ti$idsMZ2, drop = FALSE])^2)/2; # var of MZ pair differences / 2
        Rd = rowMeans((gene[,ti$idsDZ1, drop = FALSE] - gene[,ti$idsDZ2, drop = FALSE])^2)/2; # var of DZ pair differences / 2
        Ra = rowMeans(gene^2);
        
        # head(cbind(Rm, Rd, Ra))
    } # Rm, Rd, Ra
    
    acelist = vector('list', 4);
    
    # Estimate model with A,C,E all included
    {
        # Full model
        #             2*E = c(0,0,2) %*% rez = Rm
        #   A       + 2*E = c(1,0,2) %*% rez = Rd
        # 2*A + 2*C + 2*E = c(2,2,2) %*% rez = Ra
        
        # mat = matrix(nrow = 3, ncol = 3, data = c(
        #   c(0,0,2),
        #   c(1,0,2),
        #   c(2,2,2)))
        
        acelist[[1]] = cbind(
            A = 2 * (Rd - Rm),
            C = Rm + Ra - 2 * Rd,
            E = Rm);
        
        # head(round(acelist[[1]],1))
    } # acelist[[1]]
    
    # Estimate model with C,E included
    {
        #       2*E = c(0,2) %*% rez = 2 * Rm
        #       2*E = c(0,2) %*% rez = 2 * Rd
        # 2*C + 2*E = c(2,2) %*% rez = 2 * Ra
        
        #       2*E = c(0,2) %*% rez = 2 * (Rm*Nmz + Rd*Ndz) / (Nmz + Ndz)
        # 2*C + 2*E = c(2,2) %*% rez = 2 * Ra
        
        R2 = (Rm * ti$Nm + Rd * ti$Nd) / (ti$Nm + ti$Nd);
        acelist[[2]] = cbind(
            A = 0,
            C = Ra - R2,
            E = R2);
        rm(R2);
        
        # head(round(acelist[[2]],1));
    } # acelist[[2]]
    
    # Estimate model with A,E included (approximate, very close)
    {
        #       2*E = c(0,2) %*% rez = 2 * Rm
        #   A + 2*E = c(1,2) %*% rez = 2 * Rd
        # 2*A + 2*E = c(2,2) %*% rez = 2 * Ra
        
        # weighted regression
        tmpA = 2 * ((Ra - Rd) * ti$Nd + (Ra - Rm) * ti$Nm*2) / (ti$Nd + 4 * ti$Nm);
        acelist[[3]] = cbind(
            A = tmpA,
            C = 0,
            E = Ra - tmpA);
        rm(tmpA);
        
        # head(round(acelist[[3]],1));
    } # acelist[[3]]
    
    # Estimate model with E only
    {
        acelist[[4]] = cbind(
            A = 0,
            C = 0,
            E = Ra);
        # head(round(acelist[[4]],1));
    } # acelist[[4]]
    
    # Select best model for each gene
    {
        # weights from ACE parameters to 2*c(Rm, Rd, Ra)
        ace2rrr = matrix(nrow = 3, ncol = 3, data = c(
            c(0,0,2),
            c(1,0,2),
            c(2,2,2)));
        
        # best group averages
        bstmean = 2 * cbind(Rm, Rd, Ra);
        
        
        bestace = matrix(NA_real_, nrow = nrow(gene), ncol = 3);
        bestmft = rep(+Inf, nrow(gene));
        
        for( i in seq_along(acelist) ){ # i = 1;
            
            # Current model estimates
            ace = acelist[[i]];
            
            # calculate the misfit (SSE for the model - SSE of best fit).
            mft = (ace %*% ace2rrr - bstmean)^2 %*% c(ti$Nm, ti$Nd, (ti$N^2 - ti$N)/2 - ti$Nm - ti$Nd);
            mft[rowSums(ace<0) > 0L] = +Inf;
            
            set = which(mft < bestmft);
            if(length(set)>0){
                bestace[set,] = ace[set,];
                bestmft[set ] = mft[set ];
            }
            rm(ace, set, mft);
        }
        rm(i);
        rm(ace2rrr, bstmean, bestmft);
    }
 
    rm(acelist);
    rm(Rm, Rd, Ra);
    rm(ti);
    # rm(idsMZ1, idsMZ2, idsDZ1, idsDZ2, Nm, Nd, N);

    return(bestace);   
}

.ProcessTwininfo = function(twininfo, colnms){
    
    twininfo[[1]] = as.character(twininfo[[1]]);
    twininfo[[2]] = as.character(twininfo[[2]]);
    
    # Check twininfo
    {
        stopifnot( !any(duplicated(c(twininfo[[1]],twininfo[[2]]))) );
        stopifnot( all(c(twininfo[[1]],twininfo[[2]]) %in% colnms) );
        stopifnot( all(twininfo[[3]] %in% c("DZ","MZ")) );
    }

    twininfo1 = match(twininfo[[1]], colnms)
    twininfo2 = match(twininfo[[2]], colnms)
    
    MZset = which(twininfo[[3]] == "MZ");
    idsMZ1 = twininfo1[MZset];
    idsMZ2 = twininfo2[MZset];
    
    DZset = which(twininfo[[3]] == "DZ");
    idsDZ1 = twininfo1[DZset];
    idsDZ2 = twininfo2[DZset];
    
    idsS = seq_along(colnms);
    idsS = idsS[!(idsS %in% c(idsMZ1,idsMZ2,idsDZ1,idsDZ2))];
    
    Nm = length(MZset);
    Nd = length(DZset);
    N = length(colnms);
    rm(twininfo1, twininfo2, MZset, DZset);
    
    return(list(idsMZ1 = idsMZ1, idsMZ2 = idsMZ2, idsDZ1 = idsDZ1, idsDZ2 = idsDZ2, idsS = idsS, Nm = Nm, Nd = Nd, N = N));
}

# Generate artificial data set for testing the package
TwinMeta_simulate = function(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt, ACEparam = NULL, MAF = NULL){

    # Default parameters
    if(FALSE){
        Nm = 1000
        Nd = 2000
        Ns = 3000
        
        Ngene = 10000
        Nsnps = 100000
        Ncvrt = 10
        
        ACEparam = c(3,4,5);
        ACEparam = matrix(data = runif(Ngene*3), nrow = Ngene, ncol = 3);
        
        MAF = 0.2;
        MAF = runif(Nsnps)^2*0.49+0.01;
            
        # sim = TwinMeta_simulate( Nm = 10, Nd = 10, Ns = 10, Ngene = 1000, Nsnps = 1000, Ncvrt = 10)
    }
    
    # Checks and defaults
    {
        stopifnot(length(Nm) == 1);
        stopifnot(is.numeric(Nm));
        stopifnot(is.finite(Nm));
        stopifnot(Nm >= 0);
        stopifnot(Nm == as.integer(Nm))
        
        stopifnot(length(Nd) == 1);
        stopifnot(is.numeric(Nd));
        stopifnot(is.finite(Nd));
        stopifnot(Nd >= 0);
        stopifnot(Nd == as.integer(Nd))

        stopifnot(length(Ns) == 1);
        stopifnot(is.numeric(Ns));
        stopifnot(is.finite(Ns));
        stopifnot(Ns >= 0);
        stopifnot(Ns == as.integer(Ns))
        
        stopifnot(length(Ngene) == 1);
        stopifnot(is.numeric(Ngene));
        stopifnot(is.finite(Ngene));
        stopifnot(Ngene >= 1);
        stopifnot(Ngene == as.integer(Ngene))
        
        stopifnot(length(Nsnps) == 1);
        stopifnot(is.numeric(Nsnps));
        stopifnot(is.finite(Nsnps));
        stopifnot(Nsnps >= 1);
        stopifnot(Nsnps == as.integer(Nsnps))
        
        stopifnot(length(Ncvrt) == 1);
        stopifnot(is.numeric(Ncvrt));
        stopifnot(is.finite(Ncvrt));
        stopifnot(Ncvrt >= 0);
        stopifnot(Ncvrt == as.integer(Ncvrt))
        
        stopifnot(Nm + Nd + Ns >= 1);
        
        if( is.null(ACEparam) ){
            A = 3 + 3 * runif(Ngene);
            C = 4 + 4 * runif(Ngene);
            E = 5 + 4 * runif(Ngene);
        } else {
            stopifnot( is.numeric(ACEparam) );
            stopifnot(all( ACEparam >= 0 ));
            
            if( is.matrix(ACEparam) ){
                
                stopifnot( nrow(ACEparam) == Ngene );
                stopifnot( ncol(ACEparam) == 3 );
                
                A = ACEparam[,1];
                C = ACEparam[,2];
                E = ACEparam[,3];
                
            } else {
                
                stopifnot( length(ACEparam) == 3 )
                
                A = ACEparam[1];
                C = ACEparam[2];
                E = ACEparam[3];
                
            }
        }
        
        if( is.null(MAF) ){
            MAF = runif(Nsnps)^2*0.49+0.01;
        } else {
            stopifnot( is.numeric(MAF) );
            stopifnot(all( MAF >= 0 ));
            stopifnot(all( MAF <= 1 ));
            stopifnot( length(MAF) %in% c(1, Nsnps) );
        }
    } # A, C, E, MAF
    
    # Make integers
    {
        Nm = as.integer(Nm);
        Nd = as.integer(Nd);
        Ns = as.integer(Ns);
        Ngene = as.integer(Ngene);
        Nsnps = as.integer(Nsnps);
        Ncvrt = as.integer(Ncvrt);
    } # Nm, Nd, Ns, Ngene, Nsnps, Ncvrt
    
    # Total number of samples
    N = 2L * Nm + 2L * Nd + Ns;
    
    # Generate set ids
    {
        idsMZ1 = seq_len(Nm)
        idsMZ2 = seq_len(Nm) + Nm;
        idsDZ1 = seq_len(Nd) + 2*Nm;
        idsDZ2 = seq_len(Nd) + 2*Nm + Nd;
        idsS   = seq_len(Ns) + 2*Nm + 2*Nd;
        
        # stopifnot(c(idsMZ1,idsMZ2,idsDZ1,idsDZ2,idsS) == seq_len(N));
        # stopifnot(tail(idsS,1) == N)
    } # idsMZ1, idsMZ2, idsDZ1, idsDZ2, idsS  

    # Generate Gene expression
    {
        message("Gene: Start generating gene expression")
        gene = matrix(NA_real_, Ngene, N);
        
        # cov(tmp1, tmp2) = sqrt(vr) * (cv/vr) * sqrt(vr) = cv;
        # var( (cv/vr) * tmp1 ) = vr * (cv/vr)^2 = cv^2 / vr

        message("Gene: Generating MZ sample pairs")
        vr = A + C + E;
        cv = A + C;
        for( i in seq_along(idsMZ1) ){ # i = 1
            tmp1 = rnorm(Ngene) * sqrt(vr);
            tmp2 = (cv/vr) * tmp1 + sqrt(vr - cv^2 / vr) * rnorm(Ngene) ;
            gene[, idsMZ1[i]] = tmp1;
            gene[, idsMZ2[i]] = tmp2;
            rm(tmp1, tmp2);
        }
        rm(i, cv);
        
        message("Gene: Generating DZ sample pairs")
        cv = A/2 + C;
        for( i in seq_along(idsDZ1) ){ # i = 1
            tmp1 = rnorm(Ngene) * sqrt(vr);
            tmp2 = (cv/vr) * tmp1 + sqrt(vr - cv^2 / vr) * rnorm(Ngene) ;
            gene[, idsDZ1[i]] = tmp1;
            gene[, idsDZ2[i]] = tmp2;
            rm(tmp1, tmp2);
        }
        rm(i, cv);
        
        message("Gene: Generating singleton samples")
        for( i in seq_along(idsS) ){ # i = 1
            gene[, idsS[i]] = rnorm(Ngene) * sqrt(vr);
        }
        rm(i);

        message("Gene: Setting row and column names")
        colnames(gene) = sprintf("Sample_%06d", seq_len(N));
        rownames(gene) = sprintf("Gene_%06d", seq_len(Ngene));
        
        rm(vr);
        rm(A, C, E);

        message("Gene: Done generating gene expression")
    } # gene

    # Generate genotypes
    {
        message("SNPs: Start generating genotypes")
        snps = matrix(NA_integer_, Nsnps, N);

        message("SNPs: Generating DZ sample pairs")
        for( i in seq_along(idsDZ1) ){ # i = 1
            a1 = as.integer(runif(Nsnps) < MAF);
            a2 = as.integer(runif(Nsnps) < MAF);
            a3 = as.integer(runif(Nsnps) < MAF);

            snps[, idsDZ1[i]] = a1 + a2;
            snps[, idsDZ2[i]] = a2 + a3;
            rm(a1, a2, a3);
        }
        rm(i);
        
        message("SNPs: Generating MZ sample pairs")
        # probs: MAF^2, 2*MAF*(1-MAF), (1 - MAF)^2
        # 2: MAF^2
        # 1: MAF^2 + 2*MAF*(1-MAF) = 2 * MAF - MAF^2
        # 0: Otherwise
        thr2 = MAF^2;
        thr1 = 2 * MAF - MAF^2;
        for( i in seq_along(idsMZ1) ){ # i = 1
            
            tmp = runif(Nsnps);
            snp = integer(Nsnps);
            snp[tmp < thr1] = 1L;
            snp[tmp < thr2] = 2L;
                  
            snps[, idsMZ1[i]] = snp;
            snps[, idsMZ2[i]] = snp;
            rm(tmp, snp);
        }
        rm(i);
        
        message("SNPs: Generating singleton samples")
        for( i in seq_along(idsS) ){ # i = 1
            
            tmp = runif(Nsnps);
            snp = integer(Nsnps);
            snp[tmp < thr1] = 1L;
            snp[tmp < thr2] = 2L;
            
            snps[, idsS[i]] = snp;
            rm(tmp, snp);
        }
        rm(i, thr1, thr2);
        
        
        # min(snps)
        # 
        # cor(c(snps[,idsMZ1]),c(snps[,idsMZ2]))
        # cor(c(snps[,idsDZ1]),c(snps[,idsDZ2]))
        # 
        # cor(c(snps[,c(idsMZ1,idsDZ1)]),c(snps[,idsS]))
        
        message("SNPs: Setting row and column names")
        colnames(snps) = colnames(gene);
        rownames(snps) = sprintf("Snp_%06d", seq_len(Nsnps));
        
        rm(MAF);
        
        message("SNPs: Done generating genotypes")
    } # snps
    
    # Covariates
    {
        message("Cvrt: Start generating covariates") 
        cvrt = rnorm(Ncvrt*N);
        dim(cvrt) = c(Ncvrt, N);
        colnames(cvrt) = colnames(gene);
        rownames(cvrt) = sprintf("Covt_%06d", seq_len(Ncvrt));
        message("Cvrt: Done generating covariates") 
    } # cvrt
    
    # Generate twininfo
    {
        twininfo = data.frame(
            first  = colnames(gene)[c(idsMZ1,idsDZ1)],
            second = colnames(gene)[c(idsMZ2,idsDZ2)],
            type = rep(c("MZ","DZ"), c(Nm, Nd)));
    } # twininfo
    
    return(list(gene = gene, snps = snps, cvrt = cvrt, twininfo = twininfo));
}

# Main function for eQTL testing on twins
TwinMeta_testAll = function(gene, snps, cvrt, twininfo, pvthreshold){
    
    # Checks
    {
        message("Checking input parameters") 
        
        # gene
        stopifnot( is.matrix(gene) );
        stopifnot( !is.null(rownames(gene)) );
        stopifnot( !is.null(colnames(gene)) );
        stopifnot( is.numeric(gene) );
        stopifnot( all(is.finite(gene)) );
        stopifnot( nrow(gene) > 0 );
        stopifnot( ncol(gene) > 0 );
        
        # SNPs
        stopifnot( is.matrix(snps) );
        stopifnot( !is.null(rownames(snps)) );
        stopifnot( !is.null(colnames(snps)) );
        stopifnot( is.numeric(snps) );
        stopifnot( all(is.finite(snps)) );
        stopifnot( nrow(snps) > 0 );
        stopifnot( ncol(snps) > 0 );
        
        # cvrt
        if( NROW(cvrt) > 0 ){
            stopifnot( is.matrix(cvrt) );
            stopifnot( !is.null(rownames(cvrt)) );
            stopifnot( !is.null(colnames(cvrt)) );
            stopifnot( is.numeric(cvrt) );
            stopifnot( all(is.finite(cvrt)) );
            
            stopifnot( ncol(gene) == ncol(cvrt) );
            stopifnot( all(colnames(gene) == colnames(cvrt)) );
        }
        
        # columns
        stopifnot( ncol(gene) == ncol(snps) );
        stopifnot( all(colnames(gene) == colnames(snps)) );
        
        # pvthreshold
        stopifnot( length(pvthreshold) == 1 );
        stopifnot( is.numeric(pvthreshold) );
        stopifnot( is.finite(pvthreshold) );
        stopifnot( pvthreshold > 0 );
    }
    
    # Process twininfo
    {
        message("Parsing twininfo parameter") 
        ti = .ProcessTwininfo(twininfo = twininfo, colnames(gene));
    } # ti with (idsMZ1, idsMZ2, idsDZ1, idsDZ2, idsS, Nm, Nd, N)

    # Perform ACE model estimation, using SqD method
    {
        ace = EstimateACE_SqD(gene = gene, cvrt = cvrt, twininfo = twininfo);
        
        ace = ace / rowSums(ace);
        
        # rhoMZ = ace %*% c(1,   1, 0);
        # rhoDZ = ace %*% c(0.5, 1, 0);
        # corrT1T2 = (Nd * rhoDZ + 2 * Nm * rhoMZ) / N;
        corrT1T2 = as.vector(ace %*% ((ti$Nd * c(1, 1, 0) + 2 * ti$Nm * c(0.5, 1, 0)) / ti$N));
        
        # cor(as.vector(tt1), as.vector(tt2))
        # [1] 0.2177539
        # > mean(corrT1T2)
        # [1] 0.2186725
        
        ttmultiplier = 1 / sqrt(2 + 2 * corrT1T2);
        rm(corrT1T2, ace);
    } # ttmultiplier
    
    # Split the samples into 2 groups, each without related individuals
    {
        message("Splitting data into two groups of independent samples") 
        
        # split the data (Avoid R stupidity with ifs
        if( length(ti$idsS) == 0 ){
            idsS1 = c();
            idsS2 = c();
        } else if( length(ti$idsS) == 1 ){
            idsS1 = ti$idsS;
            idsS2 = c();
        } else {
            idsS1 = ti$idsS[seq(from = 1, to = length(ti$idsS), by = 2)];
            idsS2 = ti$idsS[seq(from = 2, to = length(ti$idsS), by = 2)];
        }
        
        ids1 = c(ti$idsMZ1, ti$idsDZ1, idsS1);
        ids2 = c(ti$idsMZ2, ti$idsDZ2, idsS2);
        rm(idsS1, idsS2);
        
        stopifnot( all(sort(c(ids1, ids2)) == 1:ncol(gene)) );
        
        gene1 = gene[, ids1, drop = FALSE];
        gene2 = gene[, ids2, drop = FALSE];

        if( NROW(cvrt) > 0){
            cvrt1 = cvrt[, ids1, drop = FALSE];
            cvrt2 = cvrt[, ids2, drop = FALSE];
        } else {
            cvrt1 = NULL;
            cvrt2 = NULL;
        }
    } # ids1, ids2, gene1, cvrt1, gene2, cvrt2

    # Preprocess set 1 (gene1, cvrt1)
    {
        message("Residualizing gene expression data in set 1") 

        zcvrt1 = .orthonormalizeCovariates(cvrt1, gene1);
        gene1 = .ResidualizeSlice(gene1, zcvrt1);
    } # zcvrt1, gene1
    
    # Preprocess set 2 (gene2, cvrt2)
    {
        message("Residualizing gene expression data in set 2") 
        
        zcvrt2 = .orthonormalizeCovariates(cvrt2, gene2);
        gene2 = .ResidualizeSlice(gene2, zcvrt2);
    } # zcvrt2, gene2
    
    # Loop over SNPs, in blocks
    {
        # Functions consistent through the loop
        {
            # Correlation to t-statistic for slice1
            dfFull1 = ncol(gene1) - 1 - nrow(zcvrt1);
            testfn1 = function(x){ return(x * sqrt(dfFull1 / (1 - pmin(x^2,1)))); }
            
            # Correlation to t-statistic for slice2
            dfFull2 = ncol(gene2) - 1 - nrow(zcvrt2);
            testfn2 = function(x){ return(x * sqrt(dfFull2 / (1 - pmin(x^2,1)))); }
            
            # Threshold for abs(z)
            absthr =  qnorm(pvthreshold/2, lower.tail = FALSE);
            absthr[pvthreshold >= 1] = 0;
            absthr[pvthreshold <= 0] = Inf;

            # Get p-value from abs(z)
            pvfun = function(absz){
                return( .pv.nz(pnorm(absz, lower.tail = FALSE)*2));
            }
            
        } # dfFull1, testfn1, dfFull2, testfn2, absthr, pvfun
        
        blocksize = 1024L;
        Nsnps = nrow(snps);
        nsteps = ceiling(Nsnps/blocksize);
        
        # Create .listBuilder objects to collect results
        {
            collect.geneid = vector('list', nsteps);
            collect.snpsid = vector('list', nsteps);
            collect.zscore = vector('list', nsteps);
            collect.pvalue = vector('list', nsteps);
        } # collect.geneid, collect.snpsid, collect.zscore, collect.pvalue
        
        for( part in seq_len(nsteps) ){ # part = 1L
            fr = (part-1L)*blocksize + 1L;
            to = min(part*blocksize, Nsnps);
            
            message('Testing SNPs ', .s(fr), ' - ', .s(to), ' of ', .s(Nsnps), '. ', round(100*to/Nsnps,2), '%');
            
            
            # Extract SNP slices
            {
                slice1 = snps[fr:to, ids1, drop = FALSE];
                slice2 = snps[fr:to, ids2, drop = FALSE];
            } # slice1, slice2
            
            # Residualize and standardize slice1, slice2
            {
                slice1 = .ResidualizeSlice(slice1, zcvrt1);
                slice2 = .ResidualizeSlice(slice2, zcvrt2);
            } # slice1, slice2
            
            # Getting t-statistics for slice1, slice2, combining them
            {
                tt1 = testfn1(tcrossprod(gene1, slice1));
                tt2 = testfn2(tcrossprod(gene2, slice2));
                rm(slice1, slice2);
                
                zstat = (tt1 + tt2) * ttmultiplier;
                rm(tt1,tt2);
                abszz = abs(zstat);
                # qqnorm(as.vector(zstat))
            } # zstat, abszz
            
            # Select z-stats exceeding threshold, record them in collect*
            {
                selected = which(abszz > absthr, arr.ind = TRUE, useNames = FALSE);
                if( length(selected) > 0 ){
                    collect.geneid[[part]] = selected[,1];
                    collect.snpsid[[part]] = selected[,2] + (fr - 1L);
                    collect.zscore[[part]] = zstat[selected];
                    collect.pvalue[[part]] = pvfun(abszz[selected]);
                }
                rm(selected);
                rm(zstat, abszz)
            } # collect.* updated
        }
        rm(part, blocksize, Nsnps, nsteps, fr, to);
        rm(dfFull1, testfn1, dfFull2, testfn2, absthr, pvfun)
    } # collect*
    
    rm(ttmultiplier);
    rm(gene1, cvrt1, gene2, cvrt2, zcvrt1, zcvrt2); # snps1, snps2
    rm(ids1, ids2);
    
    # Form the resulting data frame
    {
        message('Collecting findings in a data frame')
        
        # gene names factor
        gene.factor = unlist(collect.geneid, recursive = FALSE, use.names = FALSE);
        if( length(gene.factor) == 0 ){
            # message('No significant associations found at pvthreshold = ', pvthreshold);
            warning('No significant associations found at pvthreshold = ', pvthreshold);
            return(data.frame(gene = character(), snps = character(), zscore = numeric(0), pvalue = numeric(0)));
        }
        levels(gene.factor) = rownames(gene);
        class(gene.factor) = 'factor';
        rm(collect.geneid);
        
        # SNPs names factor
        snps.factor = unlist(collect.snpsid, recursive = FALSE, use.names = FALSE);
        levels(snps.factor) = rownames(snps);
        class(snps.factor) = 'factor';
        rm(collect.snpsid);
        
        # Z-scores
        zscore = unlist(collect.zscore, recursive = FALSE, use.names = FALSE);
        rm(collect.zscore);
        
        # P-values
        pvalue = unlist(collect.pvalue, recursive = FALSE, use.names = FALSE);
        rm(collect.pvalue);
        
        results = data.frame(
            gene = gene.factor,
            snps = snps.factor,
            zscore = zscore,
            pvalue = pvalue);
        rm(gene.factor, snps.factor, zscore, pvalue);
    } # results
    
    return(results);
}
