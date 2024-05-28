
setwd('~/projects/Tcell-pancancer/cox')

out = read.csv('out.csv',header=T,row.names=1)


tran = function(temp)
{
  temp_p = c()
  for(i in temp)
  {
    if(i == 0) {temp_p = append(temp_p,0)}
    if(i > 1) {temp_p = append(temp_p,2)}
    if(i < 1&&i != 0) {temp_p = append(temp_p,1)}
  }
  return(temp_p)
}

out_temp = sapply(out,tran)
rownames(out_temp) = rownames(out)
out = as.data.frame(out_temp)

temp_c = c()
for(i in colnames(out))
{
  eval(parse(text=paste0('temp_c = c(temp_c,out$',i,')')))
}


# data1 <- data.frame(
#   x = rep(1:length(colnames(out)), each = length(out$ACC)),
#   y = rep(1:length(out$ACC), times = length(colnames(out))),
#   value = temp_c,
#   text = rownames(out),
#   fill_column = rep(c(rep("no_fill",(length(colnames(out))-1)), "fill"), each = length(out$ACC))
# )
# 
# color_scale <- scale_fill_manual(
#   name = NULL,
#   breaks = c("Not", "Protect","Risky"),
#   labels = c("Not   ", "Protect   ","Risky"),
#   values = c(
#     "Not" = "white",
#     "Protect" = "steelblue1",
#     "Risky" = "#BB362F"
#   ),
#   guide = guide_legend(override.aes = list(size = 0))
# )
# 
# pdf('cox queue.pdf',width=80,height=5)
# 
# ggplot(data1, aes(x = x, y = y, fill = value)) +
#   geom_tile(width = 0.9, height = 0.9,color = "black", size = 0.2) + 
#   geom_line(aes(x=1,y=1))+
#   theme_void()+geom_text(aes(label = text,x=x+0.01,y=y),data = subset(data1, fill_column == "fill"), color = "black", size = 3)+
#   #theme_void()+geom_text(aes(label = sample,x=x,y=y+0.01),data1 = subset(data1,fill_row), color = "black", size = 6)+
#   scale_x_continuous(expand = c(4, 0))+
#   color_scale
# 
# dev.off()

myColor <- colorRampPalette(c( "white","steelblue1" ,"#BB362F"))(3)

library(pheatmap)
pdf('test.pdf')
col_group = data.frame(Tumor=colnames(out))
rownames(col_group) = colnames(out)
ann_color_temp = c('#007bbb', '#4d4398', '#c38743', '#ad7e4e', '#192f60', '#884898', '#3eb370', '#80989b', '#c8c2c6', '#93ca76', '#e09e87', '#e5e4e6', '#a58f86', '#726250', '#a89dac', '#9790a4', '#824880', '#eae5e3', '#954e2a', '#bc763c', '#eb6ea5', '#ffec47', '#e17b34', '#c37854', '#b77b57', '#e4ab9b', '#eec362', '#9f563a')
ann_colors = c()
for(i in 1:length(colnames(out)))
{
  eval(parse(text=paste0('ann_colors = append(ann_colors,','c(',colnames(out)[i],'="',ann_color_temp[i],'"))')))
}
ann_colors = list(Tumor=ann_colors)
pheatmap(out,cluster_rows = FALSE,cluster_cols = FALSE,cellwidth = 15,cellheight = 8,annotation_col=col_group,annotation_colors=ann_colors,labels_col = F,border_color = "black",color=myColor)
dev.off()
