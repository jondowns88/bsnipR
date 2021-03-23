#Install the BiocManager
install.packages("BiocManager")
BiocManager::install(version = "3.12")

#Use biocmanager to install sva and bladderbatch packages
BiocManager::install(c("sva", "bladderbatch"))

#Load new libs
library(sva)
library(bladderbatch)
options(stringsAsFactors=FALSE)

#(Load bladder dataset)
data(bladderdata)
pheno = pData(bladderEset)
# add fake age variable for numeric
pheno$age = c(1:7, rep(1:10, 5))
write.table(data.frame(cel=rownames(pheno), pheno), row.names=F, quote=F, sep="\t", file="bladder-pheno.txt")

edata = exprs(bladderEset)
write.table(edata, row.names=T, quote=F, sep="\t", file="bladder-expr.txt")
# use dataframe instead of matrix
mod = model.matrix(~as.factor(cancer) + age, data=pheno)
t = Sys.time()

#(JD) The numCov argument seems to be no more
cdata = ComBat(dat=edata, batch=as.factor(pheno$batch), mod=mod #, numCov=match("age", colnames(mod)))
print(Sys.time() - t)
print(cdata[1:5, 1:5])

write.table(cdata, "r-batch.txt", sep="\t", quote=F)