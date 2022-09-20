library(data.table)
library(stringr)
library(magrittr)
library(abind)
library(pheatmap)

for(gl in c('kegg', 'hmark')){
for(do.sl.rdn in c(F, T)){
# genesets download from MSigDB
if(gl == 'kegg'){
  c2 = fread('gmt/c2.all.v7.5.1.symbols.gmt', head= F, sep = '/')[, c(1,7)]
  c2$V1 = gsub('\\thttp:', '', c2$V1)
  c2$V7 = str_split(c2$V7, '\\t') %>% lapply(function(i) i[-1])
  colnames(c2) = c('gsName', 'glist')
  c2 = c2[grep('KEGG_', gsName)] #
} else {
  c2 = fread(paste0('gmt/h.all.v7.5.1.symbols.gmt'), head= F, sep = '/')[, c(1,7)]
  c2$V1 = gsub('\\thttp:', '', c2$V1)
  c2$V7 = str_split(c2$V7, '\\t') %>% lapply(function(i) i[-1])
  colnames(c2) = c('gsName', 'glist')
}
# paper supplementary table
sl.gs0 = fread('data/SLlist.csv', head= F)
sl.gs0 = sl.gs0$V1 %>% sort() %>% str_split('-', simplify = T) %>% .[, c(1:2)] %>% data.table()
TOTN = 1
if(do.sl.rdn){
  TOTN = 200
}
res2 = lapply(1:TOTN, function(SMPN){
  sl.gs = sl.gs0
  if(do.sl.rdn){
    set.seed(SMPN)
    sl.gs$V2 = sample(sl.gs$V2, nrow(sl.gs), replace = F)
  }
  kg = c2$gsName
  n = kg %>% length()
  res4 = matrix(0, nrow = n, ncol = n)
  rownames(res4) = kg
  colnames(res4) = kg
  tmp = lapply(1:nrow(sl.gs), function(n){
    gp = sl.gs[n] %>% as.character()
    r = lapply(1:nrow(c2), function(j){
      c(gp[1] %in% c2$glist[[j]], gp[2] %in% c2$glist[[j]], j, c2$gsName[j])
    }) %>% do.call(rbind, .) %>% data.table()
    colnames(r) = c('i1', 'i2', 'j', 'gsn')
    x = r$j[r$i1 == 'TRUE'] %>% as.numeric()
    y = r$j[r$i2 == 'TRUE'] %>% as.numeric()
    if(length(x) * length(y) > 0){
      res4[x, y] <<- res4[x, y] + 1
      res4[y, x] <<- res4[y, x] + 1
    }
  })
  res4
})
if(do.sl.rdn){
  save(file= paste0('sl-', gl, '-ov-rdn-mat.rda'), res2)
} else {
  res4 = res2[[1]]
  save(file= paste0('sl-', gl, '-ov-mat.rda'), res4)
}
}
load(paste0('sl-', gl, '-ov-mat.rda'))
load(paste0('sl-', gl, '-ov-rdn-mat.rda'))
res4 = res4[colnames(res2[[1]]), colnames(res2[[1]])]
res3 = do.call(function(...) abind(..., along= 3), res2)
zs = lapply(1:nrow(res4), function(i){
  sapply(1:ncol(res4), function(j){
    (res4[i, j] - mean(res3[i, j, ])) / (1 + sd(res3[i, j, ]))
  })
}) %>% do.call(rbind, .)
colnames(zs) = colnames(res4)
rownames(zs) = rownames(res4)
pheatmap(zs)
save(file= paste0('sl-', gl, '-enrich.rda'), zs)
}

