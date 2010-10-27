# �Q�ߓI�ɕ��z�Ɏ��R�Ȑ����
# ADF: Asymptotically Distribution Free estimator

library(MASS)

# nlm()�ւ̑�����
ADF.fmin <- function(theta, s, w){
    # ���f���̖��m���i�O���ϐ��ɑΉ�)
    a   <- theta[1]
    b   <- theta[2]
    Mux <- theta[3]
    Sx2 <- theta[4]
    Se2 <- theta[5]
    Sx3 <- theta[6]
    Se3 <- theta[7] 

    # ���f���̍Đ��l(�����ϐ��ɑΉ�)
    Muy  <- a   * Mux + b
    Sxy  <- a   * Sx2
    Sy2  <- a^2 * Sx2 + Se2

    Sx2y <- a   * Sx3
    Sxy2 <- a^2 * Sx3
    Sy3  <- a^3 * Sx3 + Se3

    S <- matrix(c(Mux, Muy, Sx2, Sxy, Sy2, Sx3, Sx2y, Sxy2, Sy3))
    d <- (s - S)

    t(d) %*% w %*% (d) # text p.7, expr 1.39
}

ADF.reg <- function(x, y){
    n     <- length(x)

    mux   <- mean(x)
    muy   <- mean(y)
    dx    <- (x - mux)
    dy    <- (y - muy)

    sx2y0 <- mean(dx^2)
    sx0y2 <- mean(dy^2)
    sx1y1 <- mean(dx   * dy )

    sx3y0 <- mean(dx^3 * dy^0)
    sx2y1 <- mean(dx^2 * dy^1)
    sx1y2 <- mean(dx^1 * dy^2)
    sx0y3 <- mean(dx^0 * dy^3)

    sx4y0 <- mean(dx^4 * dy^0)
    sx3y1 <- mean(dx^3 * dy^1)
    sx2y2 <- mean(dx^2 * dy^2)
    sx1y3 <- mean(dx^1 * dy^3)
    sx0y4 <- mean(dx^0 * dy^4)

    sx5y0 <- mean(dx^5 * dy^0)
    sx4y1 <- mean(dx^4 * dy^1)
    sx3y2 <- mean(dx^3 * dy^2)
    sx2y3 <- mean(dx^2 * dy^3)
    sx1y4 <- mean(dx^1 * dy^4)
    sx0y5 <- mean(dx^0 * dy^5)

    sx6y0 <- mean(dx^6 * dy^0)
    sx5y1 <- mean(dx^5 * dy^1)
    sx4y2 <- mean(dx^4 * dy^2)
    sx3y3 <- mean(dx^3 * dy^3)
    sx2y4 <- mean(dx^2 * dy^4)
    sx1y5 <- mean(dx^1 * dy^5)
    sx0y6 <- mean(dx^0 * dy^6)

    # �d�ݕt��
    v <- matrix(c(
        sx2y0,    # w[x,x]
        sx1y1,    # w[y,x]
        sx3y0,    # w[sx2,x]
        sx2y1,    # w[sxy,x]
        sx1y2,    # w[sy2,x]
        sx4y0 - 3 * sx2y0^2,            # w[sx3,x]
        sx3y1 - 3 * sx2y0 * sx1y1,        # w[sx2y1, x]
        sx2y2 - sx2y0 * sx0y2 - 2 * sx1y1^2,    # w[sxy2,x]
        sx1y3 - 3 * sx0y2 * sx1y1,        # w[sy3,x]

        sx0y2,    # w[y,y]
        sx2y1,    # w[sx2,y]
        sx1y2,    # w[sxy,y]
        sx0y3,    # w[sy2,y]
        sx3y1 - 3 * sx2y0 * sx1y1,        # w[sx3,y]
        sx2y2 - sx2y0 * sx0y2 - 2 * sx1y1^2,    # w[sx2y,y]
        sx1y3 - 3 * sx0y2 * sx1y1,        # w[sxy2,y]
        sx0y4 - 3 * sx0y2^2,            # w[sy3,y]

        sx4y0 - sx2y0^2,        # w[sx2,sx2]
        sx3y1 - sx2y0 * sx1y1,        # w[sxy,sx2]
        sx2y2 - sx2y0 * sx0y2,        # w[sy2,sx2]
        sx5y0 - 4 * sx3y0 * sx2y0,    # w[sx3,sx2]
        sx4y1 - 2 * sx2y1 * sx2y0 - 2 * sx3y0 * sx1y1,            # w[sx2y,sx2]
        sx3y2 - sx1y2 * sx2y0 - 2 * sx2y1 * sx1y1 - sx3y0 * sx0y2,    # w[sxy2,sx2]
        sx2y3 - sx2y0 * sx0y3 - 3 * sx2y1 * sx0y2,            # w[sy3,sx2]

        sx2y2 - sx1y1^2,                # w[sxy,sxy]
        sx1y3 - sx1y1 * sx0y2,                # w[sy2,sxy]
        sx4y1 - sx3y0 * sx1y1 - 3 * sx2y0 * sx2y1,    # w[sx3,sxy]
        sx3y2 - sx2y0 * sx1y2 - 3 * sx2y1 * sx1y1,    # w[sx2y,sxy]
        sx2y3 - sx0y2 * sx2y1 - 3 * sx1y2 * sx1y1,    # w[sxy2,sxy]
        sx1y4 - sx0y3 * sx1y1 - 3 * sx0y2 * sx1y2,    # w[sy3,sxy]

        sx0y4 - sx0y2^2,                # w[sy2,sy2]
        sx3y2 - sx0y2 * sx3y0 - 3 * sx1y2 * sx2y0,    # w[sx3,sy2]
        sx2y3 - sx2y1 * sx0y2 - sx0y3 * sx2y0 - 2 * sx1y2 * sx1y1,    # w[sx2y,sy2]
        sx1y4 - 2 * sx1y2 * sx0y2 - 2 * sx0y3 * sx1y1,        # w[sxy2,sy2]
        sx0y5 - 4 * sx0y3 * sx0y2,                # w[sy3,sy2]

        sx6y0 - 6 * sx4y0 * sx2y0 - sx3y0^2 + 9 * sx2y0^3,    # w[sx3,sx3]
        sx5y1 - 4 * sx2y0 * sx3y1 - sx3y0 * sx2y1
            - 2 * sx4y0 * sx1y1 + 9 * sx2y0^2 * sx1y1,    # w[sx2y,sx3]
        sx4y2 - 3 * sx2y0 * sx2y2 - sx3y0 * sx1y2 - 2 * sx3y1 * sx1y1
            - sx4y0 * sx0y2 + 6 * sx2y0 * sx1y1^2 + 3 * sx2y0^2 * sx0y2,    # w[sxy2,sx3]
        sx3y3 - 3 * sx3y1 * sx0y2 - sx3y0 * sx0y3 - 3 * sx2y0 * sx1y3
            + 9 * sx2y0 * sx0y2 * sx1y1,            # w[sy3,sx3]

        sx4y2 - 4 * sx3y1 * sx1y1 - 2 * sx2y2 * sx2y0 - sx2y1^2
            + 8 * sx2y0 * sx1y1^2 + sx2y0^2 * sx0y2,     # w[sx2y,sx2y]
        sx3y3 - sx2y0 * sx1y3 - 4 * sx2y2 * sx1y1 - sx0y2 * sx3y1
            - sx2y1 * sx1y2 + 5 * sx2y0 * sx1y1 * sx0y2 + 4 * sx1y1^3,    # w[sxy2,sx2y]
        sx2y4 - 3 * sx0y2 * sx2y2 - sx0y3 * sx2y1 - 2 * sx1y3 * sx1y1
            - sx2y0 * sx0y4 + 6 * sx0y2 * sx1y1^2 + 3 * sx0y2^2 * sx2y0,    # w[sy3,sx2y]

        sx2y4 - 4 * sx1y3 * sx1y1 - 2 * sx2y2 * sx0y2 - sx1y2^2
            + 8 * sx0y2 * sx1y1^2 + sx0y2^2 * sx2y0,    # w[sxy2,sxy2]
        sx1y5 - 4 * sx0y2 * sx1y3 - sx0y3 * sx1y2
            - 2 * sx0y4 * sx1y1 + 9 * sx0y2^2 * sx1y1,    # w[sy3,sxy2]

        sx0y6 - 6 * sx0y4 * sx0y2 - sx0y3^2 + 9 * sx0y2^3    # w[sy3,sy3]
    ))

    # �Ίp���̔��Α������߂�
    w  <- matrix(0, 9, 9)
    ic <- 0
    for(j in 1:9){
        for(i in j:9){
            ic <- ic + 1
            w[i, j] = v[ic]
        }
    }

    # �Ίp���̒l��␳
    w <- w + t(w) - diag(diag(w))



    # �W�{�̐ϗ�(1��,2��,3��)��Vector 
    # �� �cVector�ɂ��邽�߂�matrix()�����Ă���
    s <- matrix(c(mux, muy, sx2y0, sx1y1, sx0y2, sx3y0, sx2y1, sx1y2, sx0y3))

    # ���m���̏����l (a, b, Mux, Sx2, Se2, Sx3, Se3)
    regr  <- lm(y~x)
    a_0   <- regr$coeff[2]
    b_0   <- regr$coeff[1]

    se2_0 <- mean((regr$resid)^2)
    se3_0 <- mean((regr$resid)^3)

    theta_0 <- matrix(c(a_0, b_0, mux, sx2y0, se2_0, sx3y0, se3_0))
    rnlm    <- nlm(ADF.fmin, theta_0, s, ginv(w), iterlim=1000, hessian=TRUE)
    retval <- list();

    # ��^2 
    #chi2 <- rnlm$minimum * n

    # �C����^2 (Yuan-Bentler proposal) (Yuan & Bentler, 1997) 
    # �� X^2���z�ɂ��߂��ߎ��l�ɂ��邽�߁C��^2�l�����Ȃ߂ɏC������i���p�������Ȃ��Ȃ�j
    chi2 <- rnlm$minimum * n / (1 + rnlm$minimum)

    # H0: �uy <- x ���������v
    # H1: �uH0�ł͂Ȃ��v
    # ���R�x df=2 (�ϗ�9�� - ���m��2��)
    df <- 2
    # p�l
    p <- pchisq(df=df, chi2, lower.tail=F)

    # �W���덷�̐��� (p.13)
    #se <- sqrt( diag(solve(rnlm$hessian)) / n )

    # ���f���̓K�����̐�Ε]��
    # RMSEA
#    sqrt( max( (chi2/n) / df - (1/ (n-1)), 0) )

    # ���m���̐��茋��
#    retval$a   <- rnlm$estimate[1]
#    retval$b   <- rnlm$estimate[2]
#    retval$Mux <- rnlm$estimate[3]
#    retval$Sx2 <- rnlm$estimate[4]
#    retval$Se2 <- rnlm$estimate[5]
#    retval$Sx3 <- rnlm$estimate[6]
#    retval$Se3 <- rnlm$estimate[7]

    retval$rnlm      <- rnlm

    retval$chi2      <- chi2
    retval$df        <- df
    retval$pvalue    <- p
    #retval$se        <- se

    # ���f���̓K�����̑��Ε]��
    # �l���������قǓK���x������

    # AIC: Akaike's Information Criterion �Ԓr���ʊ
#    retval$AIC <- chi2 - 2 * df

    # CAIC: Consistent AIC
#    retval$CAIC <- chi2 - (log(n)+1) * df

    # SBC: Schwarz Bayesian Criterion
#    retval$SBC <- chi2 - log(n) * df

    return(retval)
}

skewness <- function(v){
    mean( (v - mean(v))^3 ) / sd(v)^3
}

ADF.plot <- function(x, y, filename = "", title = "", width = 800, height = 600){

    if(filename != ""){ # �t�@�C�������w�肳��Ă���΃t�@�C���ɏo��
        png(filename = filename, width = width, height = height,
            units = "px", pointsize = 12, bg = "white", res = NA)
    }

    plot(x, y, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)));

    retval = ADF.reg(x, y);

    abline(b=retval$rnlm$estimate[1], a=retval$rnlm$estimate[2], col="blue")

    b <- retval$rnlm$estimate[1] # �P��A�W������l

# �ʏ�̒P��A����
#    rlm = lm(y~x)
#    abline(a=rlm$coefficients[1], b=rlm$coefficients[2], col="red", lty=4)

    sign <- ""
    if(retval$pvalue < 0.01){
        sign <- "**"
    }
    else if(retval$pvalue < 0.05){
        sign <- "*"
    }

    notes <- sprintf("n=%d, b=%.02f, skewness of x=%.02f and y=%.02f, chi2(%d)=%5.02f, p=%.03f%s",
                length(x), b, skewness(x), skewness(y), retval$df, retval$chi2, retval$pvalue, sign)

    title(main = paste(title, notes, sep="\n"));

    if(filename != ""){
        dev.off();
    }

    retval$notes <- notes

    return(retval)
}