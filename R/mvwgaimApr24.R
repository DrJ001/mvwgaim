
mvwgaim <- function (baseDiag, ...)
UseMethod("mvwgaim")

mvwgaim.default <- function(baseDiag, ...)
  stop("Currently the only supported method is \"asreml\"")

mvwgaim.asreml <- function (baseDiag, baseModel, genObj, merge.by = NULL, fix.lines = TRUE, Trait = NULL, gen.type = "interval", method = "fixed", selection = "interval", n.fa = 0, exclusion.window = 20, breakout = -1, TypeI = 0.05, trace = TRUE, verboseLev = 0, ...)
{
    asreml.options(Cfixed = TRUE)
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
        if(!baseModel$converge)
            stop("Base model not converged: Check base model before proceeding with QTL analysis.")
    }
    asremlEnv <- lapply(baseModel$formulae, function(el) attr(el, ".Environment"))
    phenoData <- eval(baseModel$call$data)
    if (missing(phenoData))
        stop("phenoData is a required argument.")
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "interval"))
        stop("genObj is not of class \"interval\"")
    if (is.null(merge.by))
        stop("Need name of matching column to merge datasets.")
    if (is.null(glines <- genObj$pheno[, merge.by]))
        stop("Genotypic data does not contain column \"", merge.by,
             "\".")
    if (is.null(plines <- phenoData[, merge.by]))
        stop("Phenotypic data does not contain column \"", merge.by,
             "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in genotypic \"", merge.by, "\" column do not match any names in phenotypic \"",
             merge.by, "\" column.")
    if (!(method %in% c("fixed","random")))
        stop("Method has to be either \"fixed\" or \"random\" (see ?wgaim.asreml).")
    if (!(selection %in% c("interval","chromosome")))
        stop("Selection method has to be either \"interval\" or \"chromosome\" (see ?wgaim.asreml).")
    if(!is.numeric(breakout) | breakout < -1 | breakout == 0)
        stop("breakout argument must be -1 or a positive integer.")
    if(is.null(Trait))
        stop("Trait needs to be non-NULL.")
    if((n.trait <- length(levels(phenoData[, Trait]))) > 2) {
        n.par.fa <- (n.fa+1)*n.trait - n.fa*(n.fa-1)/2
        n.par.us <- n.trait*(n.trait+1)/2
        if(n.par.fa > n.par.us)
            stop('n.fa set too high: reset and try again\n')
    }
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }
    if(gen.type %in% "interval")
        gdat <- lapply(genObj$geno, function(el) el$interval.data)
    else gdat <- lapply(genObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    nint <- lapply(gdat, function(el) 1:ncol(el))
    lint <- unlist(lapply(nint, length))
    mnams <- paste("Chr", rep(names(genObj$geno), times = lint), unlist(nint), sep = ".")
    dimnames(genoData) <- list(as.character(glines), mnams)
    genoData <- genoData[rownames(genoData) %in% as.character(plines),]
    rterms <- unlist(strsplit(deparse(baseModel$call$random[[2]], width.cutoff = 500), " \\+ "))
    plabs <- c(merge.by, Trait)
    sgrep <- paste(c("(",paste(plabs[1:2], collapse = ".*"),"|",paste(plabs[2:1], collapse = ".*"),")"), collapse = "")
    sterms <- baseModel$factor.names
    pterm <- rterms[grep(sgrep, rterms)]
    rterms <- rterms[!(rterms %in% pterm)]
    dterms <- unlist(strsplit(deparse(baseDiag$call$random[[2]], width.cutoff = 500), " \\+ "))
    dterm <- dterms[grep(sgrep, dterms)]
    whg <- levels(phenoData[[merge.by]]) %in% rownames(genoData)
    genetic.term <- merge.by
    vm <- FALSE
    if(!all(whg) & fix.lines){
#        !all(whg <- levels(phenoData[[merge.by]]) %in% rownames(genoData)))
        phenoData$Gomit <- phenoData$Gsave <- plines
        levels(phenoData$Gsave)[!whg] <- NA
        levels(phenoData$Gomit)[whg] <- "GEN"
        fix.form <- as.formula(paste(". ~ ", paste(Trait, "Gomit", sep = ":")," + .", sep = ""))
        pterm <- gsub(merge.by, "Gsave", pterm)
        ran.base <- formula(paste("~ ", paste(c(pterm, rterms), collapse = " + "), sep = ""))
        baseModel$call$data <- quote(phenoData)
        cat("\nFixing lines and updating initial base model:\n")
        cat("============================================\n")
        baseModel <- update(baseModel, fixed. = fix.form, random. = ran.base, ...)
        dterm <- gsub(merge.by, "Gsave", dterm)
        ran.based <- formula(paste("~ ", paste(c(dterm, rterms), collapse = " + "), sep = ""))
        baseDiag$call$data <- quote(phenoData)
        baseDiag <- update(baseDiag, fixed. = fix.form, random. = ran.based, ...)
        merge.by <- "Gsave"
    }
    message("\nQTL x ",Trait," Diagonal Random effects model.")
    cat("========================================\n")
    qtlModel <- baseModel
    if(ncol(genoData) > nrow(genoData)){
        cov.env <- wgaim:::constructCM(genoData)
        covObj <- cov.env$relm
        qterm <- paste("vm","(",merge.by,", covObj):","diag(", Trait, ")", sep = "")
        ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
        attr(genObj, "env") <- cov.env
        vm <- TRUE
    } else {
        covObj <- cbind.data.frame(rownames(genoData), genoData)
        names(covObj)[1] <- merge.by
        qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
        qtlModel$call$mbf$ints$cov <- "covObj"
        qterm <- paste("mbf('ints'):", "diag(", Trait, ")", sep = "")
        ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
    }
    assign("covObj", covObj, envir = parent.frame())
    qtlModel$call$data <- quote(phenoData)
    qtlModel <- update(qtlModel, random. = ran.form, ...)
    ## added diag model test

    ran.formd <- as.formula(paste(c("~", qterm, dterm, rterms), collapse = " + "))
    qtlDiag <- update(baseDiag, random. = ran.formd, ...)
    LRT <- 2*(qtlModel$loglik - baseDiag$loglik)
    pvalue <- 1 - pchisq.mixture(LRT, ntrait = n.trait)
    if(pvalue > TypeI) {
        cat(' Residual likelihood ratio test: p-value = ', pvalue, '\n')
        cat("No QTL to be selected - no significant Trait QTL variation\n")
        return()
    }
    if(n.fa > 0){
        qsp <- strsplit(qterm, ":")
        rhs <- sapply(qsp, "[", 2)
        if(n.trait == 2){
            rhs <- gsub("diag", "corh", rhs)
            qterm <- paste(sapply(qsp, "[", 1), rhs, sep = ":")
            ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
            message("\nQTL x ",Trait," Bivariate Random Effects model.")
            cat("===============================================\n")
            qtlModel <- update(qtlModel, random. = ran.form, ...)
        } else {
            for(i in 1:n.fa){
                if(length(grep("diag", rhs))){
                    rhs <- gsub("diag", "fa", rhs)
                    rhs <- gsub(")", ",1)", rhs)
                } else {
                    old <- paste(i - 1, ")", sep = "")
                    new <- paste(i, ")", sep = "")
                    rhs <- gsub(old, new, rhs)
                }
                qterm <- paste(sapply(qsp, "[", 1), rhs, sep = ":")
                ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
                message("\nQTL x ",Trait," Factor Analytic (",i,") Random Effects model.")
                cat("===================================================\n")
                qtlModel <- update(qtlModel, random. = ran.form, ...)
            }
        }
    }
    ldiag <- iclist <- coef.list <- vcoef.list <- list()
    qtl <- c(); iter <- 1
    state <- rep(1, ncol(genoData))
    names(state) <- mnams
    repeat {
        cat("\nQTL Selection Iteration (", iter,"):\n")
        cat("============================\n")
        assign("qtlModel", qtlModel, .GlobalEnv)
        assign("phenoData", phenoData, .GlobalEnv)
        selq <- qtlMSelect(qtlModel, phenoData, genObj, gen.type, selection, n.fa, Trait, exclusion.window, state, verboseLev)
    #    baseModel <- qtlModel
        state <- selq$state
        ldiag$oint[[iter]] <- selq$oint
        ldiag$ochr[[iter]] <- selq$ochr
        ldiag$blups[[iter]] <- selq$blups
        qtl[iter] <- selq$qtl
        cqtl <- strsplit(qtl[iter], "\\.")
        wchr <- sapply(cqtl, "[", 2)
        wint <- sapply(cqtl, "[", 3)
        message("Selecting QTL on chromosome ", wchr, " ", gen.type, " ", wint)
        tmp <- cbind.data.frame(rownames(genoData), genoData[, qtl[iter]])
        qtl.x <- gsub("Chr\\.", "X.", qtl[iter])
        names(tmp) <- c(merge.by, qtl.x)
        phenoData <- cbind.data.frame(ord = 1:nrow(phenoData), phenoData)
        phenoData <- merge(phenoData, tmp, by = merge.by, all.x = TRUE, all.y = FALSE)
        phenoData <- phenoData[order(phenoData$ord), ]
        phenoData <- phenoData[, -2]
        mout <- (1:ncol(genoData))[!as.logical(state)]
        genoSub <- genoData[, -mout]
        if(ncol(genoSub) > nrow(genoSub)){
            cov.env <- wgaim:::constructCM(genoSub)
            covObj <- cov.env$relm
            attr(genObj, "env") <- cov.env
        }
        else {
            covObj <- cbind.data.frame(rownames(genoSub), genoSub)
            names(covObj)[1] <- merge.by
            if(is.null(qtlModel$call$mbf$ints) & vm){
                attr(genObj, "env") <- NULL
                rterms <- unlist(strsplit(deparse(qtlModel$call$random[[2]]), " \\+ "))
                rterms <- rterms[!(rterms %in% qterm)]
                qsp <- strsplit(qterm, ":")
                rhs <- sapply(qsp, "[", 2)
                qterm <- paste(c("vm","(",merge.by,", covObj):", rhs), sep = "")
                qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
                qtlModel$call$mbf$ints$cov <- "covObj"
                ran.form <- as.formula(paste(c("~ ", qterm, rterms), collapse = " + "))
                qtlModel$call$random <- ran.form
            }
        }
        assign("covObj", covObj, envir = parent.frame())
        qtlModel$call$data <- baseModel$call$data <- quote(phenoData)
        fix.terms <- paste(c(qtl.x, paste(Trait, qtl.x, sep = ":")), collapse = " + ")
#        fix.terms <- paste(Trait, qtl.x, sep = ":")
        fix.form <- as.formula(paste(". ~ . +", fix.terms, sep = ""))
        cat("\nQTL x ",Trait," Fixed Effects Model Iteration (", iter, "):\n")
        cat("=================================================\n")
        qtlModel <- update(qtlModel, fixed. = fix.form, ...)
        list.coefs <- qtlModel$coefficients$fixed
        zind <- grep("X\\.", rownames(list.coefs))
        sub.list <- rev(list.coefs[zind, 1])
        names(sub.list) <- rev(rownames(list.coefs)[zind])
        coef.list[[iter]] <- sub.list
        vcoef.list[[iter]] <- rev(qtlModel$vcoeff$fixed[zind])

        cat("\nDiagonal QTL and Polygenic Model Iteration (", iter, "):\n")
        cat("================================================\n")
        qtlDiag$call$data <- baseDiag$call$data <- quote(phenoData)
        baseDiag <- update(baseDiag, fixed. = fix.form, ...)
        qtlDiag <- update(qtlDiag, fixed. = fix.form, ...)
        LRT <- 2*(qtlDiag$loglik - baseDiag$loglik)
        pvalue <- 1 - pchisq.mixture(LRT, ntrait=n.trait)
        cat("\nLikelihood Ratio Test Statistic: ", LRT, ", P-value: ", pvalue,"\n")
        iter <- iter + 1
        if(pvalue > TypeI)
           break
    }
    qtl.list <- list()
    qtl.list$selection <- selection
    qtl.list$method <- method
    qtl.list$type <- gen.type
    qtl.list$trait <- Trait
    qtl.list$trait.levels <- levels(phenoData[[Trait]])
    qtl.list$diag <- ldiag
    if (length(qtl)) {
        list.coefs <- qtlModel$coefficients$fixed
        trms <- attr(terms(qtlModel$call$fixed), "term.labels")
        marks <- trms[grep("X\\.", trms)]
        imarks <- marks[grep(":", marks)]
        mmarks <- marks[-grep(":", marks)]
        int.test <- list()
        for(i in 1:length(imarks)){
            forw <- paste(Trait,".*", mmarks[i], sep = "")
            reve <- paste(mmarks[i],".*", Trait, sep = "")
            zind <- grep(paste(forw, reve, sep = "|"), rownames(list.coefs))
            ci <- list.coefs[zind, 1]
            int.test[[i]] <- list(coef = zind[ci != 0], type = "zero")
        }
        wt <- waldTest(qtlModel, cc = int.test)
        final.terms <- ifelse(wt$Zero[,2] > 0.05, mmarks, imarks)
        other.terms <- trms[-grep("X\\.", trms)]
        fix.form <- as.formula(paste(". ~ ", paste(c(other.terms, final.terms), collapse = " + "), sep = ""))
        qtlModel <- update(qtlModel, fixed. = fix.form, ...)
#        qtl.list$diag$ic.list <- icbind
        qtl.list$diag$coef.list <- coef.list
        qtl.list$diag$vcoef.list <- vcoef.list
        qtl.list$diag$state <- state
        qtl.list$diag$genetic.term <- genetic.term
        qtl.list$diag$rel.scale <- 1
        if(exists("cov.env")) qtl.list$diag$rel.scale <- cov.env$scale
        qtl.list$iterations <- iter - 1
        qtl.list$breakout <- ifelse(breakout != -1, TRUE, FALSE)
        qtl.list$qtl <- qtl
        qtl.list$effects <- coef.list[[iter - 1]]
        qtl.list$veffects <- vcoef.list[[iter - 1]]
    }
    data.name <- paste(as.character(qtlModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = parent.frame())
    qtlModel <- wgaim:::envFix(qtlModel, asremlEnv)
    ## baseModel$call$data <- as.name(data.name)
    qtlModel$QTL <- qtl.list
    class(qtlModel) <- c("mvwgaim", "asreml")
    qtlModel
}

summary.mvwgaim <- function (object, genObj, LOD = TRUE, ...)
{
    if (missing(genObj))
        stop("genObj is a required argument")
    if (!inherits(genObj, "cross"))
        stop("genObj is not of class \"cross\"")
    if (is.null(object$QTL$qtl)) {
        cat("There are no significant putative QTL's above the threshold.\n")
        return(character(0))
    }
    trait <- object$QTL$Trait
    coefs <- object$coefficients$fixed
    inds <- grep("X\\.", rownames(coefs))
    nams <- rownames(coefs)[inds]
    coefs <- coefs[inds,]
    vcoef <- object$vcoeff$fixed[inds]
#    sigma2 <- ifelse(object$QTL$section, 1, object$sigma2)
    zrat <- coefs/sqrt(vcoef*object$sigma2)
    enams <- strsplit(nams, ":")
    object$QTL$effects <- sapply(enams, function(el) el[grep("X\\.", el)])
    traits <- sapply(enams, function(el){
        if(length(el) > 1) el[-grep("X\\.", el)]
        else "MAIN"
    })
    prefix <- paste(trait, "_", sep = "")
    traits <- gsub(prefix, "", traits)
    qtlm <- as.data.frame(getQTL(object, genObj))
    if(object$QTL$type == "interval")
        names(qtlm) <- c("Chromosome", "Interval", "Left Marker", "dist(cM)", "Right Marker", "dist(cM)")
    else names(qtlm) <- c("Chromosome", "Interval", "Marker", "dist(cM)")
    qtlm <- cbind.data.frame(Env = traits, qtlm)
    qtlm$Size <- round(coefs, 4)
    qtlm$Pvalue <- round(2 * (1 - pnorm(abs(zrat))), 4)
    if(LOD) qtlm$LOD <- round(0.5*log(exp(zrat^2), base = 10), 4)
    nints <- as.numeric(as.character(qtlm$Interval))
    qtlm <- qtlm[order(qtlm$Chromosome, nints, qtlm$Env),]
    qtlm
}

getQTL <- function (object, genObj)
{
    spe <- strsplit(object$QTL$effects, "\\.")
    wchr <- sapply(spe, "[", 2)
    wint <- as.numeric(sapply(spe, "[", 3))
    qtlm <- matrix(ncol = 6, nrow = length(wchr))
    for (i in 1:length(wchr)) {
        lhmark <- genObj$geno[[wchr[i]]]$map[wint[i]]
        qtlm[i,1:4] <- c(wchr[i], wint[i], names(lhmark), round(lhmark,
                                                                2))
        if (object$QTL$type == "interval") {
            if (length(genObj$geno[[wchr[i]]]$map) > 1)
                rhmark <- genObj$geno[[wchr[i]]]$map[wint[i] +
                  1]
            else rhmark <- genObj$geno[[wchr[i]]]$map[wint[i]]
            qtlm[i, 5:6] <- c(names(rhmark), round(rhmark, 2))
        }
        else qtlm <- qtlm[, -c(5:6),drop=FALSE]
    }
    qtlm
}

pchisq.mixture <- function(x, ntrait=2) {
    df <- 0:ntrait
    mixprobs <- dbinom(df, size=ntrait, prob=0.5)
    p <- c()
    for (i in 1:length(x)){
        p[i] <- sum(mixprobs*pchisq(x[i], df))
    }
    p
}

qchisq.mixture <- function(prob, ntrait=2, trace=F, maxiter=10) {

#   Iterative procedure to calculate percentiles for mixtures of chisquare variates
#
#   df the degrees of freedom for the constituent chi-squared variates
#   mixprobs  the mixing probabilities
#
    df <- 0:ntrait
    mixprobs <- dbinom(df, size=ntrait, prob=0.5)

#  Start value

    cv <- qchisq(prob, df=ntrait)
    if(trace) cat(" Starting value: ", cv, "\n")
    obj.fn <- function(cv=cv, df=cv, mixprobs=mixprobs) {
        obj <- ifelse(df[1] == 0, mixprobs[1] + sum(pchisq(cv, df=df[-1])*mixprobs[-1]), sum(pchisq(cv, df=df)*mixprobs)) - prob
        d.obj <- ifelse(df[1] == 0, sum(dchisq(cv, df=df[-1])*mixprobs[-1]), sum(dchisq(cv, df=df)*mixprobs))
        list(f=obj, df=d.obj)
    }
    convergence <- F
    iteration <- 0
    while(!convergence & iteration <= maxiter) {
        obj <- obj.fn(cv=cv, df=df, mixprobs=mixprobs)
        corr <- (obj$f/obj$df)
        cv <- cv - corr
        iteration <- iteration + 1
        if(trace) cat(" Iteration ", iteration, ": cv = ", cv, "\n")
        convergence <- abs(corr) < 0.00001
    }
    if(!convergence) warning(' non-convergence: result incorrect\n')
    cv
}

qtlMSelect <- function(asm, phenoData, genObj, gen.type, selection, n.fa, Trait, exclusion.window, state, verboseLev) {
    n.trait <- length(levels(phenoData[[Trait]]))
    sigma2 <- asm$sigma2
    if(asm$vparameters.con[length(asm$vparameters.con)] == 4)
        sigma2 <- 1
    cat("Predict step for outlier statistics \n")
    cat("====================================\n")
    sterms <- "(vm.*covObj)|(mbf.*ints)"
    rterms <- attr(terms.formula(asm$call$random), "term.labels")
    mterm <- rterms[grep(sterms, rterms)]
    cterms <- all.vars(as.formula(paste( "~", mterm)))[c(1,3)]
    oterm <- cterms[-grep(Trait, cterms)]
    spc <- c("diag","us","corgh","corh")
    spc <- spc[as.logical(sapply(spc, function(el, mterm) length(grep(el, mterm)), mterm))]
    if(length(spc)){
        rept <- paste(spc, "\\(", Trait, "\\)", sep = "")
        mterm <- gsub(rept, Trait, mterm)
    }
    pv <- predict(asm, classify = paste(cterms, collapse = ":"), only = mterm, vcov=TRUE, data = phenoData, maxit = 1)
    ord <- order(pv$pvals[[Trait]], pv$pvals[[oterm]])
    if(n.fa == 0)
        Ga <- asm$sigma2*diag(asm$vparameters[grep(sterms, names(asm$vparameters))])
    else if(n.fa > 0){
        if((n.trait != 2)) {
            vpars <- asm$vparameters[grep(sterms, names(asm$vparameters))]
            psi <- vpars[1:n.trait]
            Lam <- matrix(vpars[(n.trait + 1):((n.fa + 1) * n.trait)], ncol = n.fa)
            Ga <- asm$sigma2*(Lam %*% t(Lam) + diag(psi))
        }
        else {
            Ga <- asm$sigma2*summary(asm, vparameters = TRUE)$vparameters[[mterm]]
            Ga[1,2] <- Ga[2,1] <- Ga[1,2]*sqrt(Ga[1,1]*Ga[2,2])
        }
    }
    atilde <- pv$pvals[ord, 'predicted.value']
    atilde[is.na(atilde)] <- 0
    atilde <- matrix(atilde, ncol=n.trait, byrow=FALSE)
    pev <- pv$vcov[ord, ord]
    pev[is.na(pev)] <- 0
    Ginv <- MASS::ginv(as.matrix(Ga))
    if(!is.null(cov.env <- attr(genObj, "env"))) {
        trans <- cov.env$trans
        qtilde <- trans %*% atilde
        vatilde <- kronecker(Ga, cov.env$relm) - pev
    } else {
        trans <- diag(dim(atilde)[1])
        qtlide <- atilde
        vatilde <- kronecker(Ga, diag(nrow(atilde))) - pev
    }
    qtilde <- apply(qtilde, 1, function(el, Ginv) sum(c(t(el * Ginv)*el)), Ginv)
    vqtilde <- apply(trans, 1, function(el, Ginv, vatilde){
        tmp1 <- kronecker(diag(n.trait), el)
        tmp2 <- t(tmp1) %*% vatilde %*% tmp1
        sum(diag(Ginv %*% tmp2))
    }, Ginv, vatilde)
    gnams <- names(state)[as.logical(state)]
    names(qtilde) <- names(vqtilde) <- gnams
    oint <- ifelse(!is.na(qtilde/vqtilde), qtilde/vqtilde, 0)
    names(oint) <- gnams
    ochr <- NULL
    if(selection == "chromosome"){
        chr.names <- names(genObj$geno)
        nochr <- length(chr.names)
        allc <- sapply(strsplit(gnams, '\\.'), "[", 2)
        ochr <- c()
        for(c in 1:nochr){
            whc <- allc %in% chr.names[c]
            cqtilde <- qtilde[whc]
            nums <- cqtilde * cqtilde
            dens <- vqtilde[whc]
            ochr[c] <- ifelse(!is.na(sum(nums)/sum(dens)),sum(nums)/sum(dens),0)
        }
        names(ochr) <- chr.names
        mchr <- chr.names[ochr == max(ochr)]
        cint <- allc %in% mchr
        chri <- oint[cint]
        mint <- (1:length(chri))[chri == max(chri)]
        qtl <- names(chri)[mint]
        if(verboseLev > 0) {
            cat("\n Selection of chromosome using the AOM statistic\n")
            cat("=============================================== \n")
            for(i in 1:nochr)
                cat(" Chromosome ", chr.names[i], "Outlier Statistic ", ochr[i], "\n")
            cat("============================================= \n\n")
            cgen <- "Interval"
            if(gen.type == "marker") cgen <- "Marker"
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for(i in 1:length(chri))
                cat(cgen, names(chri)[i], "Outlier Statistic ", chri[i],"\n")
            cat("=============================================== \n\n")
        }
    } else {
        qtl <- names(oint)[oint == max(oint)]
        qsp <- unlist(strsplit(qtl, split="\\."))
        mint <- as.numeric(qsp[3]); mchr <- qsp[2]
        if(verboseLev > 0) {
            cgen <- "Interval"
            if(gen.type == "marker") cgen <- "Marker"
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for(i in 1:length(oint))
                cat(cgen, names(oint)[i], "Outlier Statistic ", oint[i],"\n")
            cat("=============================================== \n\n")
        }
    }
    ## fill out interval stats and update state
    qtl <- qtl[1]
    blups <- tint <- state
    tint[as.logical(state)] <- oint
    blups[as.logical(state)] <- qtilde/sqrt(abs(vqtilde))
    oint <- tint
    ## exclusion window
    schr <- sapply(strsplit(names(state), "\\."), "[", 2)
    wnams <- names(state)[schr %in% mchr]
    inums <- as.numeric(sapply(strsplit(wnams, "\\."),"[", 3))
    dists <- genObj$geno[[mchr]]$map
    if((gen.type == "interval") & (length(dists) > 1))
        dists <- dists[2:length(dists)] - diff(dists)/2
    dists <- dists[inums]
    exc <- wnams[abs(dists - dists[mint]) <= exclusion.window]
    state[exc] <- 0
    list(state = state, qtl = qtl, ochr = ochr, oint = oint, blups = blups)
}

fa.modify <- function(model, Trait, by) {
  which.term <- igrep(list(paste('fa\\(', Trait, sep=''), 'qtls'), names(model$G.param))
  if(length(which.term) > 0) {
    model$G.param <- favar.init(Trait, 'qtls', model$G.param)
  }
  which.term <- igrep(list(paste('fa\\(', Trait, sep=''), 'ints'), names(model$G.param))
  if(length(which.term) > 0) {
    model$G.param <- favar.init(Trait, 'ints', model$G.param)
  }
  which.term <- igrep(list(paste('fa\\(', Trait, sep=''), by), names(model$G.param))
  if(length(which.term) > 0) {
    model$G.param <- favar.init(Trait, by, model$G.param)
  }
  model
}

favar.init <- function(char1, char2, G.param) {
    which.term <- igrep(list(char1, char2), names(G.param))
    which.subterm <- grep(char1, names(G.param[[which.term]]))
    var.terms <- grep('var', names(G.param[[which.term]][[which.subterm]]$initial))
    con.terms <- G.param[[which.term]][[which.subterm]]$con[var.terms] == 'F'
    if(sum(con.terms) > 0) {
        G.param[[which.term]][[which.subterm]]$con[con.terms] <- 'P'
        G.param[[which.term]][[which.subterm]]$initial[con.terms] <- 0.001
    }
    G.param
}

igrep <- function(patterns, x, ...){
    if(length(patterns) == 1)
      grep(patterns[[1]], x, ...)
    else {
      xt <- x
      for(i in 1:length(patterns)){
        ind <- grep(patterns[[i]], x, ...)
        x <- x[ind]
        if(!length(x))
            return(integer(0))
      }
    pmatch(x, xt)
  }
 }

waldTest.asreml <- function(object, cc, keep.fac = TRUE, ...)
{
    if(oldClass(object) != "asreml")
        stop("Requires an object of class asreml\n")
    if(is.null(object$Cfixed)) {
        warning("Requires C matrix from model object. Refitting test model with argument \"Cfixed = T\"\n")
        asreml.options(Cfixed = TRUE)
        object <- update(object)
    }
    vrb <- as.matrix(object$Cfixed)
    tau <- c(object$coefficients$fixed)
    names(tau) <- rownames(object$coefficients$fixed)
    nc <- length(tau)
    sigma2 <- object$sigma2
    vrb <- vrb/sigma2
    ccnams <- names(tau)
    zdf <- cdf <- NULL
    cc <- lapply(cc, function(el, ccnams){
        if(!all(names(el) %in% c("coef","type","comp","group")))
            stop("Inappropriately named argument for comparison object.")
        if(is.numeric(el$coef)) {
            if(max(el$coef) > length(ccnams))
                stop("coefficient subscript out of bounds")
            names(el$coef) <- ccnams[el$coef]
        }
        else {
            if(any(is.na(pmatch(el$coef, ccnams))))
                  stop("Names in contrast do not match the names of coefficients of object")
            temp <- pmatch(el$coef, ccnams)
            names(temp) <- el$coef
            el$coef <- temp
        }
        el
    }, ccnams)
   ## split contrasts and other available tests
    ctype <- unlist(lapply(cc, function(el) el$type))
    if(!all(ctype %in% c("con","zero")))
        stop("Contrast types must be either \"con\" for treatment comparisons or \"zero\" for testing zero equality")
    cons <- cc[ctype %in% "con"]
    zero <- cc[ctype %in% "zero"]
    cse <- ctau <- zwtest <- cwtest <- zpval <- c()
    if(length(cons)) {
       CRows <- lapply(cons, function(el, nc){
           if(length(el) < 3){
               con <- contr.helmert(length(el$coef))[, (length(el$coef) - 1)]
               names(con) <- cnam <- names(el$coef)
               cat("Warning: default contrast being taken for", cnam, "is", con, "\n")
               row <- rep(0, nc)
               row[el$coef] <- con
               row
           }
           else {
               if(is.matrix(el$comp)) {
                   if(length(el$coef) != ncol(el$comp))
                       stop("Length of contrast does not match the number of specified coefficients")
                   cons <- split(el$comp, 1:nrow(el$comp))
                   rows <- lapply(cons, function(ell, first = el$coef, nc){
                       row <- rep(0, nc)
                       row[first] <- ell
                       row
                   }, first = el$coef, nc)
                   rows <- unlist(rows, use.names = F)
                   matrix(rows, nrow = nrow(el$comp), byrow = T)
               }
               else {
                   if(length(el$coef) != length(el$comp))
                       stop("Length of contrast does not match the number of specified coefficients")
                   row <- rep(0, nc)
                   row[el$coef] <- el$comp
                   row
               }
           }
       }, nc)
       Cmat <- do.call("rbind", CRows)
       if(!keep.fac)
           ccnams <- substring(ccnams, regexpr("\\_", ccnams) + 1, nchar(ccnams))
       cnam <- lapply(split(Cmat, 1:nrow(Cmat)), function(el, ccnams){
           namr <- ccnams[ifelse(el < 0, T, F)]
           naml <- ccnams[ifelse(el > 0, T, F)]
           c(paste(naml, collapse = ":"), paste(namr, collapse = ":"))
       }, ccnams)
       Cnam <- do.call("rbind", cnam)
       gnams <- lapply(cons, function(el){
           if(!is.null(el$group)){
               if(!any(names(el$group) %in% c("left","right")))
                   stop("group names must be \"left\" and \"right\".")
               if(is.null(el$group$left)){
                   if(is.matrix(el$comp))
                       el$group$left <- rep(NA, nrow(el$comp))
                   else el$group$left <- NA
               } else {
                   if(is.matrix(el$comp)){
                       if(length(el$group$left) == 1)
                           el$group$left <- rep(el$group$left, nrow(el$comp))
                       if(length(el$group$left) != nrow(el$comp))
                          stop("No. of group names do not match the number of comparisons in object")
                   }
               }
                if(is.null(el$group$right)){
                   if(is.matrix(el$comp))
                       el$group$right <- rep(NA, nrow(el$comp))
                   else el$group$right <- NA
               } else {
                   if(is.matrix(el$comp)) {
                       if(length(el$group$right) == 1)
                           el$group$right <- rep(el$group$right, nrow(el$comp))
                       if(length(el$group$right) != nrow(el$comp))
                          stop("No. of group names do not match the number of comparisons in object")
                   }
               }
           } else {
               if(is.matrix(el$comp))
                   el$group$left <- el$group$right <- rep(NA, nrow(el$comp))
               else el$group$left <- el$group$right <- NA
           }
           cbind(el$group$left, el$group$right)
       })
       Gnam <- do.call("rbind", gnams)
       Cnam[!is.na(Gnam[,1]), 1] <- Gnam[!is.na(Gnam[,1]),1]
       Cnam[!is.na(Gnam[,2]), 2] <- Gnam[!is.na(Gnam[,2]),2]
       for(i in 1:nrow(Cmat)) {
           varmat <- sum(Cmat[i,  ]*crossprod(vrb, t(Cmat)[, i]))
           cse[i] <- sqrt(varmat * sigma2)
           ctau[i] <- sum(Cmat[i,  ]*tau)
           cwtest[i] <- (ctau[i]/cse[i])^2
       }
       cdf <- data.frame(wald = round(cwtest, 6), pval = round(1 - pchisq(cwtest, 1), 6),
                         coef = round(ctau, 6), se = round(cse, 6))
       attr(cdf, "names") <- c("Wald Statistic", "P-Value", "Cont. Coef.", "Std. Error")
       attr(cdf, "row.names") <- paste(Cnam[,1], Cnam[,2],  sep = " vs ")
       oldClass(cdf) <- "data.frame"
   }
      if(length(zero)) {
        ZRows <- lapply(zero, function(el, nc){
            rows <- rep(rep(0, nc), length(el$coef))
            dum <- seq(0, (length(el$coef) - 1) * nc, by = nc)
            rows[el$coef + dum] <- 1
            matrix(rows, nrow = length(el$coef), byrow = T)
        }, nc)
        znam <- unlist(lapply(zero, function(el, ccnams) {
            if(is.null(el$group))
                paste(ccnams[el$coef], collapse = ":")
            else el$group
            }, ccnams))
        if(any(table(znam) > 1))
            stop("Duplicate names in group structures for zero equality tests.")
        for(i in 1:length(ZRows)) {
            varmat <- ZRows[[i]] %*% crossprod(vrb, t(ZRows[[i]]))
            Ctau <- ZRows[[i]] %*% tau
            zwtest[i] <- sum(Ctau*crossprod(solve(varmat), Ctau))/sigma2
            zpval[i] <- 1 - pchisq(zwtest[i], nrow(ZRows[[i]]))
        }
        zdf <- data.frame(wald = round(zwtest, 6), pval = round(zpval, 6))
        attr(zdf, "names") <- c("Wald Statistic", "P-Value")
        attr(zdf, "row.names") <- znam
        oldClass(zdf) <- "data.frame"
      }
    res <- list(Contrasts = cdf, Zero = zdf)
    invisible(res)
}

waldTest <- function(object, ...)
    UseMethod("waldTest")


