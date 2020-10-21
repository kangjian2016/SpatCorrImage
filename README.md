# Spatial Adaptive Varying Correlation Analysis for Multimodal Neuroimaging Data 

** An R package for performing the spatial adaptive varying correlation analysis for multimodal neuroimaging data

* ** Install and load the package **
```{r}
# Please make sure you have installed the following packages
# ggplot2, reshape, rlang, igraph
devtools::install_github("kangjian2016/BayesGPfit")
```

* ** Simulate images **
```{r}
set.seed(1000)
dim_img = c(10,10,10)
n = 50
grids <- lapply(1:3,function(i) seq(-round(dim_img[i]/2)+1,round(dim_img[i]/2),length=dim_img[i]))
img_1 <- array(rnorm(prod(dim_img)*n),dim=c(dim_img,n))
img_2 <- array(rnorm(prod(dim_img)*n),dim=c(dim_img,n))
cor_region <- create_sphere_mask(grids,radius=3)
cor_region_list <- array(cor_region,dim=dim(img_2))
img_2 <- ifelse(cor_region_list,img_1+rnorm(prod(dim_img)*n,sd=0.5),img_2)
```

* ** Compute correlation maps **
```{r}
mask <- create_sphere_mask(grids,radius=5)
cor_img <- comp_3D_images_corr(img_1,img_2,mask)
plot_3D_image_slices(cor_img,grids,c(-2,2))
```

* ** Adaptive correlation analysis**
```{r}
res <- Spat_Corr_3D_images(img_1,img_2,grids,mask,pos_prob=0.90,neg_prob=0.90,n_cor=0.5,size=10)
plot_3D_image_slices(list(cor=res$cor_img,smooth_cor=res$smooth_cor_est,pos_cluster=res$pos_cluster),grids,c(-3,3))
print(res)
```
