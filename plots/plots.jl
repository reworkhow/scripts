##############################################################################
#Comparison of prediction accuracies of two methods with different traing size
##############################################################################

#data
x1=["50","100","200","400","800","2000","4000","5000","6000","8000"]
d1=[0.341846,0.484607,0.661343,0.817135,0.90928,0.972444,0.977392,0.979423,0.982057,0.983476]
d2=[0.339902,0.47473,0.628285,0.782205,0.882108,0.959679,0.966109,0.971257,0.973503,0.975696];

#plot size
plot(size=(600,400))

#make plots with defined marker types and colors
plot!(x1,d1,line=(:dash,1),color=:red,marker=(3,0.8),label="MT-BayesC",legend = :topleft)#bottomright
plot!(x1,d2,line=(:dashdot,1),color=:blue,marker=(3,0.8),label="MT-BayesC-R")

#add xaxis,yaxis,title
xaxis!("training size",font(10, "Courier")) #grid=false
yaxis!("prediction accuracy",(0.3,1.1),0.3:0.1:1.1,font(10, "Courier")) #name,ylims,yticks,fonts
title!("Comparison of MT-BayesC and MT-BayesC-R",titlefont=font(10,"Courier"))

#add annotation
annotate!([(2.5,d1[3]+0.05,text("*",16,:red,:center))])
annotate!([(3.5,d1[4]+0.05,text("*",16,:red,:center))])
annotate!([(4.5,d1[5]+0.05,text("*",16,:red,:center))])
annotate!([(5.5,d1[6]+0.05,text("*",16,:red,:center))])
annotate!([(6.5,d1[7]+0.05,text("*",16,:red,:center))])
annotate!([(7.5,d1[8]+0.05,text("*",16,:red,:center))])
