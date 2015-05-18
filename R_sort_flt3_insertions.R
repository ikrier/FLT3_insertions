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

data=read.table("parsed_flt3_SIs_all_libs.csv",sep="\t",header=1,fill = T,stringsAsFactors = F)

heatmap(as.matrix(data[,grep("upall",colnames(data))]+data[,grep("downall",colnames(data))]),Rowv = NA, Colv = NA)

barplot(colSums(as.matrix(data[,grep("upall",colnames(data))]+data[,grep("downall",colnames(data))])),las=2)

clinically_pos=paste("upall_Lib_",c(13,17,28,31,40,49,64,72,78,80,83,84),sep="")
clinically_neg=paste("upall_Lib_",c(1,2,4,5,7,8,9,10,14,18,20,21,27,33,35,36,37,42,46,50,53,54,56,57,58,60,62,71,74,76,77,81,82),sep="")

pos=list()
range=seq(1,5001,200)
for(i in 1:length(range))
{
  pos[[i]]=which(colSums(as.matrix(data[,grep("upall",colnames(data))]+data[,grep("downall",colnames(data))]))>=range[i])
}

plot(range,lapply(pos,length),type="l",xlab="Threshold applied",ylab="Number of samples with insertions")


recall=lapply(pos,function(x){sum(clinically_pos%in%names(x))/length(clinically_pos)})
precision=lapply(pos,function(x){sum(names(x)%in%clinically_pos)/length(x)})
specificity=lapply(pos,function(x){sum(!clinically_neg%in%names(x))/length(clinically_neg)})
npv=lapply(pos,function(x){sum(!names(x)%in%clinically_pos)/length(x)})
FPtype1=lapply(pos,function(x){sum(clinically_neg%in%names(x))/length(clinically_neg)})
FNtype2=lapply(pos,function(x){sum(!clinically_pos%in%names(x))/length(clinically_pos)})

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

clinically_pos[!clinically_pos%in%names(pos[[7]])]
#We never can find 28, 31 and 80 with this filtering.
#Indeed in the filtered file they don't really show up.

positives=as.character(sapply(names(pos[[7]]),function(x){strsplit(x,"_")[[1]][3]}))
recall[[7]]
precision[[7]]
specificity[[7]]

dataout=as.data.frame(data[,c("chromosome","startrange","endrange","sequence","length")])

selectall=grep("all",colnames(data))
selectuniq=grep("uniq",colnames(data))
libnames=sapply(colnames(data)[grep("upall",colnames(data))],function(x){{strsplit(x,"_")[[1]][3]}})

for(lib in 1:length(libnames))
{
  dataout[,paste("Lib_",libnames[lib],"_all",sep="")]=rowSums(data[,selectall[lib*2+c(-1,0)]])
  dataout[,paste("Lib_",libnames[lib],"_unique",sep="")]=rowSums(data[,selectuniq[lib*2+c(-1,0)]])
  dataout[,paste("Lib_",libnames[lib],"_MutationState",sep="")]=as.numeric(rowSums(data[,selectall[lib*2+c(-1,0)]]))>=201
}

dataoutsmall=dataout[rowSums(dataout[grep("MutationState",colnames(dataout))])>0,]
dataoutsmall=rbind(dataoutsmall,c("TOTAL",rep("",4),colSums(dataout[,-(1:5)])))

for(lib in 1:length(libnames))
{
  dataoutsmall[nrow(dataoutsmall),paste("Lib_",libnames[lib],"_MutationState",sep="")]=c("Normal","Mutated")[1+(as.numeric(dataoutsmall[nrow(dataoutsmall),paste("Lib_",libnames[lib],"_all",sep="")])>=1201)]
  dataoutsmall[-nrow(dataoutsmall),paste("Lib_",libnames[lib],"_MutationState",sep="")]=c("Absent","Present")[1+as.numeric(as.logical(dataoutsmall[-nrow(dataoutsmall),paste("Lib_",libnames[lib],"_MutationState",sep="")]))]
}

write.table(dataoutsmall,file="Results_FLT3_summary.csv",sep="\t",quote=F,row.names=F)
 
#Something weird with lib84 where 447 supporting reads but 0 unique!
#42 is actually negative.