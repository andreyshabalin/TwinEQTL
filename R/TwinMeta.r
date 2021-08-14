library(MatrixEQTL)

### Format a number with comma delimited thousands
.s = function(x)formatC(x=x,digits=ceiling(log10(max(x)+1)),big.mark=',',big.interval=3);

.my.pmax = function(x, val){
    # minimum "pmin" function that can handle empty array
    if(length(x) == 0){
        return(x)
    } else {
        return(pmax.int(x,val));
    }    
}

.pv.nz = function(x){return( .my.pmax(x, .Machine$double.xmin) )}

# Accumulate the results in slices of a data frame
.listBuilder = setRefClass(".listBuilder",
	fields = list(
		dataEnv = "environment",
		n = "integer",
		count = 'numeric'
	),
	methods = list(
		initialize = function(){
			.self$dataEnv = new.env(hash = TRUE);
			.self$n = 0L;
			.self$count = 0;
# 			cumlength <<- 0;
			return(.self);
		},
		add = function(x){
			if(length(x) > 0){
				n <<- n + 1L;
# 				cumlength <<- cumlength + length(x);
				assign(paste(n), x, dataEnv );
				.self$count = .self$count + length(x);
			}
			return(.self);
		},
		set = function(i,x){
			i = as.integer(i);
			if(length(x) > 0){
				if(i>n)
					n <<- i;
				assign(paste(i), x, dataEnv );
			}
			return(.self);
		},
		get = function(i){
			return(base::get(paste(i),dataEnv));
		},
		list = function(){
			if(n==0)	return(list());
			result = vector("list",n);
			for( i in 1:n ){
				result[[i]] = .self$get(i);
			}
			return(result);
		},
		unlist = function(){
			return(base::unlist(.self$list(), recursive=FALSE, use.names = FALSE));
		},
		show = function(){
			cat(".listBuilder object.\nIternal object in TwinMeta package.\n");
			cat("Number of elements:", .self$n, "\n");
		}
	)
)

.SlicedDataOffsets = function(x){
	rez = integer(x$nSlices()+1);
	for( i in 1:x$nSlices()){
		rez[i+1] = rez[i] + NROW( x$getSliceRaw(i) );
	}
	return(rez);
}

TwinMeta_simulate = function(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt){

    # Default parameters
    if(FALSE){
        Nm = 1000
        Nd = 2000
        Ns = 3000
        
        Ngene = 10000
        Nsnps = 100000
        Ncvrt = 10
        
        # sim = TwinMeta_simulate( Nm = 10, Nd = 10, Ns = 10, Ngene = 1000, Nsnps = 1000, Ncvrt = 10)
    }
    
    # Checks
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
    }
    
    # Make integers
    {
        Nm = as.integer(Nm);
        Nd = as.integer(Nd);
        Ns = as.integer(Ns);
        Ngene = as.integer(Ngene);
        Nsnps = as.integer(Nsnps);
        Ncvrt = as.integer(Ncvrt);
    }
    
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
        message("Gene: Generating ACE parameters for each gene");        
        a = 3 + 3 * runif(Ngene);
        c = 4 + 4 * runif(Ngene);
        e = 5 + 4 * runif(Ngene);
        
        # Genotype effect component
        message("Gene: Generating A component")
        A = matrix(0, Ngene, N);
        A[,-idsMZ2] = rnorm(Ngene*(N-Nm))*a;
        
        A[,idsMZ2] = A[,idsMZ1];
        A[,idsDZ2] = A[,idsDZ1] * 0.5 + sqrt(0.75) * A[,idsDZ2];

        # cor(c(A[,idsDZ1]), c(A[,idsDZ2]))
        # cor(c(A[,idsMZ1]), c(A[,idsMZ2]))
        
        message("Gene: Generating C component")
        C = matrix(0, Ngene, N);
        C[,-c(idsMZ2,idsDZ2)] = rnorm(Ngene*(N-Nm-Nd))*c;
        
        C[,idsMZ2] = C[,idsMZ1];
        C[,idsDZ2] = C[,idsDZ1];
        
        # cor(c(C[,idsDZ1]), c(C[,idsDZ2]))
        # cor(c(C[,idsMZ1]), c(C[,idsMZ2]))
        
        message("Gene: Generating E component");
        E = rnorm(Ngene*N)*e;
        dim(E) = c(Ngene, N);
        
        message("Gene: Generating gene expression as A+C+E")
        gene = A + C + E;
        colnames(gene) = sprintf("Sample_%06d", seq_len(N));
        rownames(gene) = sprintf("Gene_%06d", seq_len(Ngene));
        
        rm(A, C, E);
        rm(a, c, e);
        
        message("Gene: Done generating gene expression")
    }
    
    # Generate genotypes
    {
        message("SNPs: Start generating genotypes")
        message("SNPs: Generating MAF for each SNP") 
        maf = runif(Nsnps)^2*0.49+0.01;
        # A1 = rbinom(Nsnps*N, size = 1, prob = maf);
        
        message("SNPs: Generating Allele 1 with given MAF");
        A1 = matrix(NA_integer_, Nsnps, N);
        A1[,-c(idsMZ2,idsDZ2)] = as.integer(runif(Nsnps*(N - Nm - Nd)) < maf);

        message("SNPs: Generating Allele 2 with given MAF");
        A2 = matrix(NA_integer_, Nsnps, N);
        A2[,-idsMZ2] = as.integer(runif(Nsnps*(N - Nm)) < maf);
        
        # system.time({ A1 = rbinom(Nsnps*N, size = 1, prob = maf); dim(A1) = c(Nsnps, N); })
        # system.time({ A1 = as.integer(runif(Nsnps*N) < maf); dim(A1) = c(Nsnps, N); })
        
        # A1[1:10,1:10]
        # as.matrix(maf[1:10])
        
        message("SNPs: Generating SNP matrix");
        snps = matrix(-1L, nrow = Nsnps, ncol = N);
        
        message("SNPs: Filling SNP matrix for MZ twins");
        snps[,idsMZ1] = A1[,idsMZ1] + A2[,idsMZ1];
        snps[,idsMZ2] = snps[,idsMZ1];
        
        message("SNPs: Filling SNP matrix for DZ twins");
        snps[,idsDZ1] = A1[,idsDZ1] + A2[,idsDZ1];
        snps[,idsDZ2] = A1[,idsDZ1] + A2[,idsDZ2];
        
        message("SNPs: Filling SNP matrix for Singletons");
        snps[,idsS  ] = A1[,idsS  ] + A2[,idsS  ];
        
        # min(snps)
        # 
        # cor(c(snps[,idsMZ1]),c(snps[,idsMZ2]))
        # cor(c(snps[,idsDZ1]),c(snps[,idsDZ2]))
        # 
        # cor(c(snps[,c(idsMZ1,idsDZ1)]),c(snps[,idsS]))
        
        colnames(snps) = colnames(gene);
        rownames(snps) = sprintf("Snp_%06d", seq_len(Nsnps));
        
        rm(maf);
        rm(A1, A2);
        
        message("SNPs: Done generating genotypes")
    }
    
    # Covariates
    {
        message("Cvrt: Start generating covariates") 
        cvrt = rnorm(Ncvrt*N);
        dim(cvrt) = c(Ncvrt, N);
        colnames(cvrt) = colnames(gene);
        rownames(cvrt) = sprintf("Covt_%06d", seq_len(Ncvrt));
        message("Cvrt: Done generating covariates") 
    }
    
    # Generate twininfo
    {
        twininfo = data.frame(
            first  = colnames(gene)[c(idsMZ1,idsDZ1)],
            second = colnames(gene)[c(idsMZ2,idsDZ2)],
            type = rep(c("MZ","DZ"), c(Nm, Nd)));
    } # twininfo
    
    return(list(gene = gene, snps = snps, cvrt = cvrt, twininfo = twininfo));
}

if(FALSE){

    devtools::install_github("andreyshabalin/TwinMeta@main")

    library(TwinMeta);
    set.seed(18090212+3)
    sim = TwinMeta_simulate( Nm = 1000, Nd = 2000, Ns = 3000, Ngene = 1000, Nsnps = 1000, Ncvrt = 10);
    gene = sim$gene; snps = sim$snps; cvrt = sim$cvrt; twininfo = sim$twininfo; pvthreshold = 1000 / (nrow(snps)*nrow(gene))
    # rm(sim)
    
}
TwinMeta_testAll = function(gene, snps, cvrt, twininfo, pvthreshold){
    
    # Checks
    {
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
        if( !is.null(cvrt) ){
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
        
        # twininfo
        stopifnot( !any(duplicated(c(twininfo[[1]],twininfo[[2]]))) );
        stopifnot( all(c(twininfo[[1]],twininfo[[2]]) %in% colnames(gene)) );
        
        stopifnot( all(twininfo[[3]] %in% c("DZ","MZ")) );
        
        # pvthreshold
        stopifnot( length(pvthreshold) == 1 );
        stopifnot( is.numeric(pvthreshold) );
        stopifnot( is.finite(pvthreshold) );
        stopifnot( pvthreshold > 0 );
    }
    
    # Process twininfo
    {
        twininfo1 = match(twininfo[[1]], colnames(gene))
        twininfo2 = match(twininfo[[2]], colnames(gene))
        
        MZset = which(twininfo[[3]] == "MZ");
        idsMZ1 = twininfo1[MZset];
        idsMZ2 = twininfo2[MZset];
        
        DZset = which(twininfo[[3]] == "DZ");
        idsDZ1 = twininfo1[DZset];
        idsDZ2 = twininfo2[DZset];
        
        idsS = 1:ncol(gene)
        idsS = idsS[!(idsS %in% c(idsMZ1,idsMZ2,idsDZ1,idsDZ2))];
        
        Nm = length(MZset);
        Nd = length(DZset);
        N = ncol(gene)
        rm(twininfo1, twininfo2, MZset, DZset);
    } # idsMZ1, idsMZ2, idsDZ1, idsDZ2, Nm, Nd, N
    
    # Process covariates
    {
        # Combine and add constant element
        if( NROW(cvrt)>0 ){
            cvrtM = rbind(matrix(1,1,ncol(gene)), cvrt);
        } else {
            cvrtM = matrix(1,1,ncol(gene));
        }
        
        # Standardize and orthogonolize via QR decomposition
        q = qr(t(cvrt));
        if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * ncol(gene) ){
            stop("Colinear or zero covariates detected");
        }
        zcvrt = t( qr.Q(q) );
        rm(q);
    } # cvrtM, zcvrt
    
    # Residualize and standardize gene expression data
    {
        rowsq1 = rowSums(gene^2);
        gene = gene - tcrossprod(gene, zcvrt) %*% zcvrt;
        rowsq2 = rowSums(gene^2);
        
        # kill rows colinear with the covariates
        delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
        if(any(delete.rows)){
            stop('Genes colinear with covariates: ', paste0(rownames(gene)[delete.rows], collapse = ', '))
        }
        gene = gene / sqrt(rowsq2);
        
        rm(rowsq1, rowsq2, delete.rows);
        rm(zcvrt);
    } # gene
    
    # Perform ACE model estimation, using SqD method
    {
        # Get mean squares
        {
            Rm = rowMeans((gene[,idsMZ1] - gene[,idsMZ2])^2)/2; # var of MZ pair differences / 2
            Rd = rowMeans((gene[,idsDZ1] - gene[,idsDZ2])^2)/2; # var of DZ pair differences / 2
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
            
            R2 = (Rm * Nm + Rd * Nd) / (Nm + Nd);
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
            tmpA = 2 * ((Ra - Rd)*Nd + (Ra - Rm)*Nm*2) / (Nd + Nm*4);
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
                mft = (ace %*% ace2rrr - bstmean)^2 %*% c(Nm, Nd, (N^2 - N)/2 - Nm - Nd);
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
        
        # cleanup
        rm(acelist);
        rm(Rm, Rd, Ra);
    } # bestace
    
    # Split the samples into 2 groups, each without related individuals
    {
        # split the data (Avoid R stupidity with ifs
        if( length(idsS) == 0 ){
            idsS1 = c();
            idsS2 = c();
        } else if( length(idsS) == 1 ){
            idsS1 = idsS;
            idsS2 = c();
        } else {
            idsS1 = idsS[seq(from = 1, to = length(idsS), by = 2)];
            idsS2 = idsS[seq(from = 2, to = length(idsS), by = 2)];
        }
        
        ids1 = c(idsMZ1, idsDZ1, idsS1);
        ids2 = c(idsMZ2, idsDZ2, idsS2);
        rm(idsS1, idsS2);
        
        stopifnot( all(sort(c(ids1, ids2)) == 1:ncol(gene)) );
        
        gene1 = gene[, ids1, drop = FALSE];
        snps1 = snps[, ids1, drop = FALSE];
        cvrt1 = cvrtM[,ids1, drop = FALSE];
        
        gene2 = gene[, ids2, drop = FALSE];
        snps2 = snps[, ids2, drop = FALSE];
        cvrt2 = cvrtM[,ids2, drop = FALSE];
    } # gene1, snps1, cvrt1, gene2, snps2, cvrt2

    # Get corr(T1,T2) from bestace
    {
        bestace = bestace / rowSums(bestace);
        
        # rhoMZ = bestace %*% c(1,   1, 0);
        # rhoDZ = bestace %*% c(0.5, 1, 0);
        # corrT1T2 = (Nd * rhoDZ + 2 * Nm * rhoMZ) / N;
        corrT1T2 = as.vector(bestace %*% ((Nd * c(1,   1, 0) + 2 * Nm * c(0.5, 1, 0)) / N));
        
        # cor(as.vector(tt1), as.vector(tt2))
        # [1] 0.2177539
        # > mean(corrT1T2)
        # [1] 0.2186725
        ttmultiplier = 1 / sqrt(2 + 2 * corrT1T2);
        rm(rhoMZ, rhoDZ, corrT1T2)
        rm(bestace);
    } # ttmultiplier
    
    # Preprocess set 1 (gene1, cvrt1)
    {
        q = qr(t(cvrt1));
        if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * ncol(cvrt1) ){
            stop("Colinear or zero covariates detected");
        }
        cvrt1 = t( qr.Q(q) );
        rm(q);
        
        
        rowsq1 = rowSums(gene1^2);
        gene1 = gene1 - tcrossprod(gene1, cvrt1) %*% cvrt1;
        rowsq2 = rowSums(gene1^2);
        
        # kill rows colinear with the covariates
        delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
        if(any(delete.rows)){
            stop('Genes colinear with covariates: ', paste0(rownames(gene1)[delete.rows], collapse = ', '))
        }
        gene1 = gene1 / sqrt(rowsq2);
        rm(rowsq1, rowsq2, delete.rows);
    }
    
    # Preprocess set 2 (gene2, cvrt2)
    {
        q = qr(t(cvrt2));
        if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * ncol(cvrt2) ){
            stop("Colinear or zero covariates detected");
        }
        cvrt2 = t( qr.Q(q) );
        rm(q);
        
        
        rowsq1 = rowSums(gene2^2);
        gene2 = gene2 - tcrossprod(gene2, cvrt2) %*% cvrt2;
        rowsq2 = rowSums(gene2^2);
        
        # kill rows colinear with the covariates
        delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
        if(any(delete.rows)){
            stop('Genes colinear with covariates: ', paste0(rownames(gene2)[delete.rows], collapse = ', '))
        }
        gene2 = gene2 / sqrt(rowsq2);
        rm(rowsq1, rowsq2, delete.rows);
    }
    
    # Create .listBuilder objects to collect results
    {
        collect.geneid = new('.listBuilder');
        collect.snpsid = new('.listBuilder');
        collect.zscore = new('.listBuilder');
        collect.pvalue = new('.listBuilder');
    } # collect.geneid, collect.snpsid, collect.zscore, collect.pvalue

    # Loop over SNPs, in blocks
    {
        # Functions consistent through the loop
        {
            # Correlation to t-statistic for slice1
            dfFull1 = ncol(	gene1 ) - 1 - nrow(cvrt1);
            testfn1 = function(x){ return( x * sqrt( dfFull1 / (1 - pmin(x^2,1)))); }
            
            # Correlation to t-statistic for slice2
            dfFull2 = ncol(	gene2 ) - 1 - nrow(cvrt2);
            testfn2 = function(x){ return( x * sqrt( dfFull2 / (1 - pmin(x^2,1)))); }
            
            # Threshold for abs(z)
            absthr =  qnorm(pvthreshold/2, lower.tail = FALSE);
            absthr[pvthreshold >= 1] = 0;
            absthr[pvthreshold <= 0] = Inf;

            # Get p-value from abs(z)
            pvfun = function(absz){
                return( .pv.nz(pnorm(absz,lower.tail = FALSE)*2));
            }
            
        } # dfFull1, testfn1, dfFull2, testfn2, absthr, pvfun
        
        blocksize = 1024L %/% 8L;
        Nsnps = nrow(snps1);
        nsteps = ceiling(Nsnps/blocksize);
        for( part in seq_len(nsteps) ){ # part = 1L
            fr = (part-1L)*blocksize + 1L;
            to = min(part*blocksize, Nsnps);
            message('Testing SNPs ', fr, ' - ', to);
            
            
            # Extract SNP slices
            {
                slice1 = snps1[fr:to, , drop = FALSE];
                slice2 = snps2[fr:to, , drop = FALSE];
            } # slice1, slice2
            
            # Residualize and standardize slice1
            {
                rowsq1 = rowSums(slice1^2);
                slice1 = slice1 - tcrossprod(slice1, cvrt1) %*% cvrt1;
                rowsq2 = rowSums(slice1^2);
                
                # kill rows colinear with the covariates
                delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
                if(any(delete.rows)){
                    slice1[delete.rows,] = 0;
                    rowsq2[delete.rows] = 1;
                }
                slice1 = slice1 / sqrt(rowsq2);
                
                rm(rowsq1, rowsq2, delete.rows);
            } # slice1
            
            # Residualize and standardize slice2
            {
                rowsq1 = rowSums(slice2^2);
                slice2 = slice2 - tcrossprod(slice2, cvrt1) %*% cvrt1;
                rowsq2 = rowSums(slice2^2);
                
                # kill rows colinear with the covariates
                delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
                if(any(delete.rows)){
                    slice2[delete.rows,] = 0;
                    rowsq2[delete.rows] = 1;
                }
                slice2 = slice2 / sqrt(rowsq2);
                
                rm(rowsq1, rowsq2, delete.rows);
            } # slice2
            
            # Getting t-statistics for slice1, slice2, combining them
            {
                tt1 = testfn1(tcrossprod(gene1,slice1));
                tt2 = testfn2(tcrossprod(gene2,slice2));
                
                zstat = (tt1 + tt2) * ttmultiplier;
                abszz = abs(zstat);
                rm(tt1,tt2);
                # qqnorm(as.vector(zstat))
            } # zstat
            
            # Select z-stats exceeding threshold, record them in collect*
            {
                selected = which(abszz > absthr, arr.ind = TRUE, useNames = FALSE);
                if( length(selected) > 0 ){
                    collect.geneid$add(selected[,1]);
                    collect.snpsid$add(selected[,2] + (fr - 1L));
                    collect.zscore$add(zstat[selected]);
                    collect.pvalue$add(pvfun(abszz[selected]));
                }
                rm(selected);
            }
            rm(slice1, slice2, zstat, abszz)
        }
        rm(part, blocksize, Nsnps, nsteps, fr, to);
        rm(dfFull1, testfn1, dfFull2, testfn2, absthr, pvfun)
    } # collect*
    
    # Form the resulting data frame
    {
        # gene names factor
        gene.factor = collect.geneid$unlist();
        levels(gene.factor) = rownames(gene);
        class(gene.factor) = 'factor';
        
        # SNPs names factor
        snps.factor = collect.snpsid$unlist();
        levels(snps.factor) = rownames(snps);
        class(snps.factor) = 'factor';
        
        results = data.frame( 
            gene = gene.factor,
            snps = snps.factor,
            zscore = collect.zscore$unlist(),
            pvalue = collect.pvalue$unlist());
    } # results
    
    return(results);
}

TwinMeta_testAll_old = function(
	snps1, gene1, cvrt1,
	snps2, gene2, cvrt2,
	pvThreshold){

	gene.offsets = .SlicedDataOffsets(gene1);
	snps.offsets = .SlicedDataOffsets(snps1);
	
	gene.names = rownames(gene1);
	snps.names = rownames(snps1);
	
	### Test all the conditions
	{
		stopifnot(all( gene.names == rownames(gene2) ))
		stopifnot(all( snps.names == rownames(snps2) ))
		
		stopifnot(all( gene.offsets == .SlicedDataOffsets(gene2) ))
		stopifnot(all( snps.offsets == .SlicedDataOffsets(snps2) ))
		
		stopifnot( ncol(gene1) == ncol(snps1) );
		stopifnot( ncol(gene2) == ncol(snps2) );
		
		if(nrow(cvrt1) > 0)
			stopifnot( ncol(cvrt1) == ncol(snps1) );
			
		if(nrow(cvrt2) > 0)
			stopifnot( ncol(cvrt2) == ncol(snps2) );
	}
	
	ngene = tail(gene.offsets,1);
	nsnps = tail(snps.offsets,1);
	nsamples = ncol(gene1) + ncol(gene2);
	dfFull = nsamples - 2 - nrow(cvrt1) - nrow(cvrt2);
# 	crthreshfun = function(pv){
# 		thr = qt(pv/2, dfFull, lower.tail = FALSE);
# 		thr = thr^2;
# 		thr = sqrt(  thr / (dfFull + thr) );
# 		thr[pv >= 1] = 0;
# 		thr[pv <= 0] = 1;
# 		return( thr );
# 	}
# 	crthreshold = crthreshfun(1e-5);
	ttthreshfun = function(pv){
		thr = qt(pv/2, dfFull, lower.tail = FALSE);
		return( thr );
	}
	ttthreshold = ttthreshfun(pvThreshold);
	pvfun = function(x){ return( (pt(-abs(x),dfFull)*2)); }

	cat('Preprocessing data set 1','\n');
	worker1 = .SingleMatrixEQTL$new(snps1, gene1, cvrt1);
	cat('Preprocessing data set 2','\n');
	worker2 = .SingleMatrixEQTL$new(snps2, gene2, cvrt2);

	cat('Running analysis','\n');
	{
		collect.geneid = new('.listBuilder1');
		collect.snpsid = new('.listBuilder1');
		collect.tvalue = new('.listBuilder1');
	
		tic = proc.time();
	
		for( gsl in 1:gene1$nSlices() ){ # gsl = 1
			# cat( ssl, 'of', snps1$nSlices(), '\n');
		 	
			tt1 = worker1$get_tstats(gsl);
			tt2 = worker2$get_tstats(gsl);
			
			sum1 = sum2 = sum11 = sum22 = sum12 = 0;
			for( ssl in seq_along(tt1) ){ # ssl = 1
				sum1  = sum1  + rowSums(tt1[[ssl]])
				sum2  = sum2  + rowSums(tt2[[ssl]])
				sum11 = sum11 + rowSums(tt1[[ssl]]^2)
				sum22 = sum22 + rowSums(tt2[[ssl]]^2)
				sum12 = sum12 + rowSums(tt1[[ssl]]*tt2[[ssl]])
			}
			crosscor = (sum12 - sum1*sum2/nsnps) / sqrt( ( sum11 - sum1^2/nsnps)*( sum22 - sum2^2/nsnps) );
			
			### Rows - genes (slice gsl)
			### Columns - snps (slice ssl)

			for( ssl in seq_along(tt1)){ # sl=1
				tt = (tt1[[ssl]]+tt2[[ssl]])/sqrt(2+2*crosscor);
				selected = which(abs(tt) > ttthreshold, arr.ind = TRUE, useNames=FALSE);
				if(NROW(selected)>0){
					collect.geneid$add(selected[,1] + gene.offsets[ssl]);
					collect.snpsid$add(selected[,2] + snps.offsets[ssl]);
					collect.tvalue$add(tt[selected]);
				}
			}
			
			cat(floor(100*gene.offsets[gsl+1]/tail(gene.offsets,1)),'% done, ',.s(collect.tvalue$count),' eQTLs found','\n',sep = '');
			flush.console()
			rm(tt1,tt2,tt,selected);
			# gc();
		}
	
		toc = proc.time();
		
		cat('Analysis finished in',(toc-tic)[3],'seconds','\n');
	}
	
	gene.factor = collect.geneid$unlist();
	levels(gene.factor) = gene.names;
	class(gene.factor) = 'factor'
	
	snps.factor = collect.snpsid$unlist();
	levels(snps.factor) = snps.names;
	class(snps.factor) = 'factor'
	
	result = data.frame( gene = gene.factor, snps = snps.factor, tstat = collect.tvalue$unlist(), 
								check.rows=FALSE, stringsAsFactors=FALSE);
	return(result);
}	
