# シミュレーションに使用した歪度を出力する
# (text p.43)

source("ADF.R")

B <- c(0.8)                           # 回帰係数(歪度には無関係なため一つ)
K <- c(2.5, 7.5)                      # 分布の歪み
N <- c(50, 100, 200, 300, 500, 1000)  # サンプルの数

simu_n  <- 10000

for(b in B){
    for(k in K){
        for(n in N){
            cat(sprintf("b = %g k = %g n = %4d  ", b, k, n));
            flush.console(); # print()の結果を強制的に表示させる

            i       <- 0
            dx      <- c()
            dy      <- c()

            repeat{
                theta <- 1/sqrt(k) # xの分散を1にするためにthetaを設定

                # ガンマ分布の乱数を作る
                x <- rgamma(n, k, theta)
                e <- sqrt(1-b^2) * rgamma(n, k, theta)
                y <- b * x + e # 定数は結果に影響しないため不要

                dx <- c(dx, skewness(x))
                dy <- c(dy, skewness(y))

                i <- i + 1

                if(i >= simu_n){
                    break # repeatループから抜ける
                }
            }

            cat(sprintf("\n\tskewness of x = %.02f(%.02f), skewness of y = %.02f(%.02f)\n",
                 mean(dx), sd(dx), mean(dy), sd(dy) ))

        }
    }
}
