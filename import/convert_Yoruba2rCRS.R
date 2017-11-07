# Conversion from hg19 mtDNA (15571 bp) to rCRS (15569)
setwd('/Users/ereznik/gitrepos/MTVariantPipeline/import/')

# based on https://www.snpedia.com/index.php/MtDNA_Position_Conversions

res = data.frame(matrix(NA,16571,1))
res[,1] = 1:16571
rownames(res) = paste('Position',c(1:16571),sep = ':')

res[paste('Position',1:309,sep = ':'),1] = 1:309

# drop 310
res['Position:310',1] = NA

# from 311 to 316, drop one
startix = which(rownames(res) == 'Position:311')
endix = which(rownames(res) == 'Position:316')
res[startix:endix,1] = res[startix:endix,1] - 1

# drop 317
res['Position:317',1] = NA

# subtract 2 from 318 to 3107
startix = which(rownames(res) == 'Position:318')
endix = which(rownames(res) == 'Position:3107')
res[startix:endix,1] = res[startix:endix,1] - 2

# subtract 1 from 3109 to 16194
startix = which(rownames(res) == 'Position:3109')
endix = which(rownames(res) == 'Position:16194')
res[startix:endix,1] = res[startix:endix,1] - 1

# drop 16195
res['Position:16195',1] = NA

# subtract 2 from 16196 to 16571
startix = which(rownames(res) == 'Position:16196')
endix = which(rownames(res) == 'Position:16571')
res[startix:endix,1] = res[startix:endix,1] - 2

colnames(res) = 'Mapping'
write.csv(res,'Yoruba2rCRS.txt')
