# 反社会的行動の縦断的データ
source("ADF.R")

#anti[1:4] 反社会的行動値
#read[1:4] 読解力テスト
#gen       性別
#homecog   家庭からの働きかけ
#subjid    通し番号

antiread <- read.table("antiread.dat", header = TRUE)


file_tmpl = "plots/antiread-i%d-j%d.png"
for(i in 1:4){
	for(j in 1:4){
		if(i != j){
			correctness = if(i < j){ "T" }else{ "F" } # 因果の正しさ(T=TRUE, F=FALSE)

			file = sprintf(file_tmpl, i, j)
			title = sprintf("antiread(%d, %d) %s", i, j, correctness)

			retval <- ADF.plot(antiread[,i], antiread[,j], filename=file, title=title)
			cat(correctness, retval$notes, "\n");

			flush.console()
		}
	}
}
