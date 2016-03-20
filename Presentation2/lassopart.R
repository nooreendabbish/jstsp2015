
## Create list of unique region abbreviations from original list
 ## Adds integers (starting with 2) to non-unique values
 setwd('~/Documents/BejingZhang')
 abbrev <- read.table("abbreviations.txt")
 abbrev[] <- lapply(abbrev, as.character)
 v <- vector("list", 0)
 v <- abbrev[[1]][1]
 for (j in 2:dim(abbrev)[1]){
     x <- 1
     y <- abbrev[[1]][j]
     while (is.element(y, v)){
        x <- x+1
        y <- paste0(abbrev[[1]][j],as.character(x), sep="")
     }
     v <- c(v,y)
 }

 table(v) # they are now all unique
 
 ## create list of downloaded folders with data
 ## there were problems reading in 954,955,956,957
 array <- c(787:839,841:851,853,857:865,867:872,875:878,
            880:938,940:953,974:978)
 array2 <-list()

metadata <- read.csv("metadata-bejingZangOnly.csv", header = TRUE)
head(metadata)
metadata <- metadata[,c("upload_data.id","upload_data.subject_pool", "upload_data.group_size")]
colnames(metadata) <- c("id","age","gender")
metadata[metadata$id == 787,]$gender

vectorgender <- data.frame(matrix(nrow=length(array), ncol=((length(v)*(length(v)-1)/2+1+1+1) )))


##need to create a vector of labels
corrlabels <- c()
for(i in 2:(length(v))){
    for (j in 1:(i-1)){
        corrlabels <- c(corrlabels,paste(v[i],v[j],sep="-"))
    }
}

colnames(vectorgender) <- c("id","age","gender",corrlabels)
levels(vectorgender$gender) <- c("Male","Female")

## Read in connectivity matrices, add row/col labels, and add 1 to diagonal
 for (i in 1:length(array)){
 array2[[i]] <- read.table(paste0(paste0("connectivity_matrix",as.character(array[i])),".txt"))
 dimnames(array2[[i]])[[1]] <- as.vector(v)
 dimnames(array2[[i]])[[2]] <- as.vector(v)
 array2[[i]] <- array2[[i]] + diag(x=1,nrow = length(array2[[i]]))
 vectorgender[i,] <- c(metadata[metadata$id == array[i], ]$id,
                       metadata[metadata$id == array[i], ]$age,
                       as.factor(metadata[metadata$id == array[i], ]$gender),
                       unlist(array2[[i]][upper.tri(array2[[i]],diag=FALSE)]))
 }


setwd('~/Dropbox/STAT 9190 Project/jstsp2015/Presentation2')

## install.packages("glmnet", repos='http://cran.us.r-project.org')
library(glmnet)

fit = glmnet(x = as.matrix(vectorgender[,-(1:3)]),
             y = as.vector(vectorgender$gender))

plot(fit, label =TRUE)

cvfit = cv.glmnet(x = as.matrix(vectorgender[,-(1:3)]),
                  y = (vectorgender$gender))

plot(cvfit)

cvfit$lambda.min

cvfit$lambda.1se

b <- as.matrix(coef(cvfit, s= "lambda.1se"))
rownames(b)[b!=0]
