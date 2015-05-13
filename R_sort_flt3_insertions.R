# =======
#   License
# =======
#   This code is released under the GNU General Public License 3.0. A copy
# of this license is in the LICENSE.txt file.
# copyright Irina Krier 2015
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


rm(list=ls())

data=read.table("parsed_flt3_SIs_all_libs.csv",sep="\t",header=1,fill = T)

heatmap(as.matrix(data[,grep("upall",colnames(data))]+data[,grep("downall",colnames(data))]),Rowv = NA, Colv = NA)

barplot(colSums(as.matrix(data[,grep("upall",colnames(data))]+data[,grep("downall",colnames(data))])),las=2)

clinically_pos=paste("upall_Lib_",c(13,17,28,31,40,49,64,72,78,80,83,84),sep="")
clinically_neg=paste("upall_Lib_",c(1,2,4,5,7,8,9,10,14,18,20,21,27,33,35,36,37,42,46,50,53,54,56,57,58,60,62,71,74,76,77,81,82),sep="")

pos=list()
range=seq(1,5001,200)
for(i in 1:length(range))
{
  pos[[i]]=which(colSums(as.matrix(data[,grep("upall",colnames(data))]+data[,grep("downall",colnames(data))]))>range[i])
}

plot(range,lapply(pos,length),type="l",xlab="Threshold applied",ylab="Number of samples with insertions")


recall=lapply(pos,function(x){sum(clinically_pos%in%names(x))/length(clinically_pos)})
precision=lapply(pos,function(x){sum(names(x)%in%clinically_pos)/length(x)})
specificity=lapply(pos,function(x){sum(!clinically_neg%in%names(x))/length(clinically_neg)})
npv=lapply(pos,function(x){sum(!names(x)%in%clinically_pos)/length(x)})
FPtype1=lapply(pos,function(x){sum(clinically_neg%in%names(x))/length(clinically_neg)})
FNtype2=lapply(pos,function(x){sum(!clinically_pos%in%names(x))/length(clinically_pos)})

plot(recall,precision,type="l")
plot(FPtype1,recall,type="l")

pdf("PrecisionRecallSPec.pdf")
plot(range,recall,type="l",ylim=c(0,1),xlab="Cutoff sum",ylab="")
points(range,precision,type="l",col=2)
points(range,specificity,type="l",col=3)
#points(range,FNtype2,type="l",col=4)
#points(range,FPtype1,type="l",col=5)
#points(range,npv,type="l",col=6)
abline(v=1201)
legend("bottomleft",legend=c("Recall","Precision","Specificity"),col=c(1,2,3),lty = 1)
dev.off()

clinically_pos[!clinically_pos%in%names(pos[[2]])]
#We never can find 28, 31 and 80 with this filtering.
#Indeed in the filtered file they don't really show up.

names(pos[[7]])
recall[[7]]
precision[[7]]
specificity[[7]]

