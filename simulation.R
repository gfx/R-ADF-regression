# シミュレーション
# (text p.43)

source("ADF.R")

B <- c(0.8, 0.5, 0.2)                 # 回帰係数
K <- c(2.5, 7.5)                      # 分布の歪み
N <- c(50, 100, 200, 300, 500, 1000)  # サンプルの数

simu_n  <- 100 # シミュレーションをそれぞれ何回行うか

#estimates_title <- c("a", "b", "Mu x", "Sx2", "Se2", "Sx3", "Se3");

for(b in B){
    for(k in K){
        for(n in N){
            cat(sprintf("b = %g k = %g n = %4d  ", b, k, n));
            flush.console(); # print()の結果を強制的に表示させる

            success <- 0
            i       <- 0

            xy_rejection <- 0
            yx_rejection <- 0

            e_xy  <- NULL # estimates for x -> y (correct)
            e_yx  <- NULL # estimates for y -> x (wrong)
            se_xy <- NULL # SE for x -> y (c)
            se_yx <- NULL # SE for y -> x (w)

            chi2_xy <- NULL
            chi2_yx <- NULL

            repeat{
                theta <- 1/sqrt(k) # xの分散を1にするためにthetaを設定

                # ガンマ分布の乱数を作る
                x <- rgamma(n, k, theta)
                e <- sqrt(1-b^2) * rgamma(n, k, theta)
                y <- b * x + e # 定数は結果に影響しないため不要

                # 推定
                rnlm_xy <- ADF.reg(x, y)
                rnlm_yx <- ADF.reg(y, x)

#                print(c(
#                    i,                                  # 実験番号
#                    (retval_xy$chi2 < retval_yx$chi2),  # モデル適合 成功(1)/失敗(0)
#                    retval_xy$rnlm$code,                # x->yのコード値 (1か2なら可)
#                    retval_yx$rnlm$code                 # y->xのコード値
#                ))

                # xy/yx双方のコード値が2以下なら次に進む
                if(rnlm_xy$code <= 2 && rnlm_yx$code <= 2){
                    retval_xy <- ADF.res(rnlm_xy);
                    retval_yx <- ADF.res(rnlm_yx);

                    # save
                    e_xy  <- rbind( e_xy,  rnlm_xy$estimate )
                    e_yx  <- rbind( e_yx,  rnlm_yx$estimate )
                    se_xy <- rbind( se_xy, retval_xy$se )
                    se_yx <- rbind( se_yx, retval_yx$se )

                    chi2_xy <- rbind( chi2_xy, retval_xy$chi2 )
                    chi2_yx <- rbind( chi2_yx, retval_yx$chi2 )

                    if(retval_xy$chi2 < retval_yx$chi2){
                        success <- success + 1
                    }

                    # 棄却数をカウントする
                    if(retval_xy$pvalue < 0.05){
                        xy_rejection <- xy_rejection + 1
                    }
                    if(retval_yx$pvalue < 0.05){
                        yx_rejection <- yx_rejection + 1
                    }

                    i <- i + 1

                    if(i >= simu_n){
                        break # repeatループから抜ける
                    }
                }
            }

            # モデルの棄却数（絶対評価）
            cat(sprintf("x->y reject = %5d x<-y reject = %5d ",
                 xy_rejection, yx_rejection))

            # モデルの妥当性の比較（相対評価）
            cat(sprintf("success = %5d / %5d\n", success, i))
            flush.console(); # print()の結果を強制的に表示させる

            # save to files

            id <- sprintf("b%.02f-k%.02f-n%d", b, k, n)
            write.csv(e_xy,  file = sprintf("data/e-xy-%s.csv", id))
            write.csv(e_yx,  file = sprintf("data/e-yx-%s.csv", id))
            write.csv(se_xy, file = sprintf("data/se-xy-%s.csv", id))
            write.csv(se_yx, file = sprintf("data/se-yx-%s.csv", id))

            write.csv(chi2_xy, file = sprintf("data/chi2-xy-%s.csv", id))
            write.csv(chi2_yx, file = sprintf("data/chi2-yx-%s.csv", id))
        }
    }
}
