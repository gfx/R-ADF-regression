# galton�̐e�q�̐g���f�[�^��ADF��K�p����
library(UsingR)
data(galton)

source("ADF.R")

# �C���`����Z���`���[�g���ɕϊ�
parent <- galton$parent * 2.54
child  <- galton$child  * 2.54

file <- paste("plots\\", "galton-p-c", ".png", sep = "");
retval <- ADF.plot(parent, child, file=file, title="Galton (parent->child)")
cat(retval$notes, "\n")

file <- paste("plots\\", "galton-c-p", ".png", sep = "");
retval <- ADF.plot(child, parent, file=file, title="Galton (child->parent)")
cat(retval$notes, "\n")