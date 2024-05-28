library(ggplot2)

setwd('~/projects/Tcell-pancancer')

colors = c('#007bbb', '#4d4398', '#c38743', '#ad7e4e', '#192f60', '#884898', '#3eb370', '#80989b', '#c8c2c6', '#93ca76', '#e09e87', '#e5e4e6', '#a58f86', '#726250', '#a89dac', '#9790a4', '#824880', '#eae5e3', '#954e2a', '#bc763c', '#eb6ea5', '#ffec47', '#e17b34', '#c37854', '#b77b57', '#e4ab9b', '#eec362', '#9f563a')
temp_final_data = final_data


custom_colors = c()
n = 1
for(i in unique(temp_final_data$CODE))
{
  eval(parse(text= paste('custom_colors <- append(custom_colors,c(',i,' = colors[',n,']))',sep='')))
  n = n +1
}



temp_c = c()
n = 1
for(i in unique(temp_final_data$CODE))
{
  temp_c = append(temp_c,rep(colors[n],length(rownames(temp_final_data[which(temp_final_data$CODE==i),]))))
  n = n +1
}

temp_final_data$color = temp_c

pdf('out.pdf',width=12,height=5)
ggplot(data = temp_final_data, aes(x = CODE, y = risk,fill=CODE))+   #指定数据集，设置坐标轴名称、类别颜色
  scale_fill_manual(values = custom_colors)+
  geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2, notch=T)+
  theme_classic()+
  theme(axis.text.x = element_text(size=20,color='black',angle=60,hjust=1),axis.text.y = element_text(size=20,color='black',hjust=1))
dev.off()
