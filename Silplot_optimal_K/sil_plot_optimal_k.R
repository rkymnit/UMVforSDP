library(cluster)    # clustering algorithms
library(factoextra) 
hp1_org<-read.csv("ethereum.csv")
hp2_org<-read.csv("hburg.csv")
hp3_org<-read.csv("mios.csv")
hp4_org<-read.csv("piet.csv")
set.seed(11)

#############################################################################

pdf("nbclust_hp1_sil.pdf",width =3.5, height = 2.5)
fviz_nbclust(hp1_org, kmeans, method='silhouette')
dev.off()

pdf("nbclust_hp2_sil.pdf",width =3.5, height = 2.5)
fviz_nbclust(hp2_org, kmeans, method='silhouette')
dev.off()

pdf("nbclust_hp3_sil.pdf",width =3.5, height = 2.5)
fviz_nbclust(hp3_org, kmeans, method='silhouette')
dev.off()

pdf("nbclust_hp4_sil.pdf",width =3.5, height = 2.5)
fviz_nbclust(hp4_org, kmeans, method='silhouette')
dev.off()


pdf("nbclust_hp3_sil.png", width = 4, height = 4, units = 'in', res = 500)
fviz_nbclust(hp3_org, kmeans, method='silhouette')
dev.off()

pdf("nbclust_hp4_sil.png", width = 4, height = 4, units = 'in', res = 500)
fviz_nbclust(hp4_org, kmeans, method='silhouette')
dev.off()


##############################################################################



png("nbclust_hp1_sil.png", width = 4, height = 4, units = 'in', res = 500)
fviz_nbclust(hp1_org, kmeans, method='silhouette')
dev.off()

png("nbclust_hp2_sil.png", width = 4, height = 4, units = 'in', res = 500)
fviz_nbclust(hp2_org, kmeans, method='silhouette')
dev.off()

png("nbclust_hp3_sil.png", width = 4, height = 4, units = 'in', res = 500)
fviz_nbclust(hp3_org, kmeans, method='silhouette')
dev.off()

png("nbclust_hp4_sil.png", width = 4, height = 4, units = 'in', res = 500)
fviz_nbclust(hp4_org, kmeans, method='silhouette')
dev.off()



pdf("nbclust_hp1_sil.pdf")
fviz_nbclust(horg, kmeans, method='silhouette')   #its working
dev.off()
pdf("nbclust_horg_gap.pdf")
fviz_nbclust(horg, kmeans, method='gap_stat')   #its working
dev.off()
pdf("nbclust_horg_wss.pdf")
fviz_nbclust(horg, kmeans, method='wss')   #its working
dev.off()