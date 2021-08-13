library(MatrixEQTL)

### Format a number with comma delimited thousands
.s = function(x)formatC(x=x,digits=ceiling(log10(max(x)+1)),big.mark=',',big.interval=3);

### Single Simplified Matrix eQTL engine
.SingleMatrixEQTL <- setRefClass( ".SingleMatrixEQTL",
	fields = list( 
		zgene = "SlicedData",
		zsnps = "SlicedData",
		zcvrt = "matrix"
	),
	methods = list(
		initialize = function( snps, gene, cvrt ){
			### Tests for consistency
			# Done before being called by the TwinMeta_testAll

			### Initialize
			{
				# Standard deviations before standardization
				gene.std = vector('list',gene$nSlices());
				snps.std = vector('list',snps$nSlices());
			}

			### Covariate Processing
			{
				# Combine and add constant element
				if( cvrt$nRows()>0 ){
					cvrt$SetNanRowMean();
					cvrt$CombineInOneSlice();
					cvrt = rbind(matrix(1,1,snps$nCols()),cvrt$getSlice(1));
				} else {
					cvrt = matrix(1,1,snps$nCols());
				}

				# Standardize and orthogonolize via QR decomposition
				q = qr(t(cvrt));
				if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * snps$nCols() ){
					stop("Colinear or zero covariates detected");
				}
				.self$zcvrt = t( qr.Q(q) );
			}
			
			### Expression
			{
				gene$SetNanRowMean();
				# Orthogonolize expression w.r.t. covariates
				for( sl in 1:gene$nSlices() ){ # sl = 1L
					slice = gene$getSlice(sl);
					rowsq1 = rowSums(slice^2);
					slice = slice - tcrossprod(slice,zcvrt) %*% zcvrt;
					rowsq2 = rowSums(slice^2);
					# kill rows colinear with the covariates
					delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
					if(any(delete.rows)){
						slice[delete.rows,] = 0;
						rowsq2[delete.rows] = 1;
					}
					div = sqrt(rowsq2); #sqrt( rowSums(slice^2) );
		# 			div[ div == 0 ] = 1;
					gene.std[[sl]] = div;
					gene$setSliceRaw(sl, slice / div);
				}
				rm(rowsq1, rowsq2, delete.rows, div);
				rm( sl, slice );
				#gene$RowRemoveZeroEps();
			}
			
			### Genotype 
			{
				snps$SetNanRowMean();
				# Orthogonolize expression w.r.t. covariates
				# status("Orthogonolizing expression w.r.t. covariates");
				for( sl in 1:snps$nSlices() ){ # sl = 1L
					slice = snps$getSlice(sl);
					rowsq1 = rowSums(slice^2);
					slice = slice - tcrossprod(slice,zcvrt) %*% zcvrt;
					rowsq2 = rowSums(slice^2);
					# kill rows colinear with the covariates
					delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
					if(any(delete.rows)){
						slice[delete.rows,] = 0;
						rowsq2[delete.rows] = 1;
					}
					div = sqrt(rowsq2); #sqrt( rowSums(slice^2) );
		# 			div[ div == 0 ] = 1;
					snps.std[[sl]] = div;
					snps$setSliceRaw(sl, slice / div);
				}
				rm(rowsq1, rowsq2, delete.rows, div);
				rm( sl, slice );
				#snps$RowRemoveZeroEps();
			}
			
			### Save preprocessed data sets
			.self$zsnps = snps;
			.self$zgene = gene;
		},
		get_tstats = function(geneslice){

			dfFull = ncol(	zsnps ) - 1 - nrow(zcvrt);
			testfun = function(x){ return( x * sqrt( dfFull / (1 - pmin(x^2,1))));	}

			gslice = zgene[[geneslice]];
			
			rez = vector('list',zsnps$nSlices());
			for( sl in 1:zsnps$nSlices() ){ # sl = 1L
				rez[[sl]] = testfun(tcrossprod(gslice,zsnps[[sl]]));
			}
			return(rez);
		}
	)
)

.listBuilder1 = setRefClass(".listBuilder1",
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
			cat(".listBuilder1 object.\nIternal object in TwinMeta package.\n");
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

# sim = TwinMeta_simulate( Nm = 1000, Nd = 2000, Ns = 3000, Ngene = 10000, Nsnps = 100000, Ncvrt = 10);
TwinMeta_simulate = function(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt){

    # Default parameters
    if(FALSE){
        Nm = 1000
        Nd = 2000
        Ns = 3000
        
        Ngene = 10000
        Nsnps = 100000
        Ncvrt = 10
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
    
    twininfo = data.frame(
        first  = colnames(gene)[c(idsMZ1,idsDZ1)],
        second = colnames(gene)[c(idsMZ2,idsDZ2)],
        type = rep(c("MZ","DZ"), c(Nm, Nd)));
    
    return(list(gene = gene, snps = snps, cvrt = cvrt, twininfo = twininfo));
}

# gene = sim$gene; snps = sim$snps; cvrt = sim$cvrt;

TwinMeta_testAll = function(gene, snps, cvrt, twininfo, pvthreshold){
   
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
