# server.R
#install.packages("kernlab")
library(shiny)
library(kernlab)
#library(edgeR)

shinyServer(function(input, output){
# print title
    # output$caption <- renderPrint({
    #           input$caption
    #    })
# print TDS reference
output$tds.reference.view<-renderTable({
        file.tds<-input$TDS_reference
        if(is.null(file.tds))
        return(NULL)
        data<-read.csv(file.tds$datapath, header=T, row.names=1)
        data.gene<-data[data[,1]!=0,] # take out reference genes
        data.gene[1:5,]
    })
# print drug logFC
output$drug.data.view<-renderTable({
        file.drug<-input$Drug_data
        if(is.null(file.drug))
        return(NULL)
        drug.table<-read.csv(file.drug$datapath, header=T, row.names=1)
        drug.table[1:5,1:5]
    })
# calculating Z-score
output$plotZScore<-renderPlot({
    logFC<-as.numeric(input$LogFC.cutoff)   # read in logFC cut-off
    file.tds<-input$TDS_reference           # read in reference TDS file
    if(is.null(file.tds))
    return(NULL)
    data<-read.csv(file.tds$datapath, header=T,row.names=1)
    data.ref<-data[data[,1]==0,]      # For reference
    data.gene<-data[data[,1]!=0,]     # For z-score calculation
    file.drug<-input$Drug_data        # read in actual drug induced logFC (tgmel+durg/tgmel)
    if(is.null(file.drug))
    return(NULL)
    drug.table<-as.matrix(read.csv(file.drug$datapath, header=T, row.names=1))
    # sign any value between cut.off and -1*cut.off as 0
    for (i in 1:dim(drug.table)[1]){
        for (j in 1:dim(drug.table)[2]){
            if (drug.table[i,j]>0){ if (drug.table[i,j] < logFC){drug.table[i,j]<-0}}
            if (drug.table[i,j]<0){ if (drug.table[i,j] > (-1)*logFC){drug.table[i,j]<-0}}
        }
    }
    ref<-cbind(rownames(data.gene),data.gene)
    drug<-cbind(rownames(drug.table),drug.table)
    colnames(drug)[1]<-"rownames(drug.table)"
    all.tab<-as.matrix(merge(x=ref, y=drug, by.x="rownames(data.gene)", by.y="rownames(drug.table)", all=F))
    # in all.tab $logFC is the the expression change of tgmel/wt
    Tg.score<-matrix(rep(0, 4*dim(all.tab)[1]),ncol=4)   # Form disease signature logFC table and model a perfect drug. The weight value requires double check. It should be the ROC value.
    colnames(Tg.score)<-c("logFC(Tg.over.wt","weight","weight.Tg.logFC","perfect.drug.logFC")
    Tg.score[,1]<-as.matrix(all.tab[,5])
    Tg.score[,2]<-all.tab[,2]
    class(Tg.score)<-"numeric"
    Tg.score[,3]<-Tg.score[,1]*Tg.score[,2]
    Tg.score[,4]<-Tg.score[,3]*(-1)
    drug.score<-matrix(rep(0,6*(dim(all.tab)[2]-5)),ncol=6) # first 5 columns are TDS registration
    rownames(drug.score)<-colnames(all.tab)[6:dim(all.tab)[2]]
    colnames(drug.score)<-c("Zscore","adjusted.Zscore","scoreA","scoreB","scoreC","scoreD")
    for (i in 1:dim(drug.score)[1]){
        tab<-cbind(Tg.score, all.tab[,i+5], rep(0,dim(all.tab)[1]),rep(0,dim(all.tab)[1]))
        class(tab)<-"numeric"
        colnames(tab)[5]<-colnames(all.tab)[i+5]
        colnames(tab)[6]<-"gene direction"
        colnames(tab)[7]<-"gene score"
        tab.mod<-tab                                 # set logFC limit to +/-5
        for (j in 1:dim(tab.mod)[1]){
            if(tab.mod[j,4]>5){tab.mod[j,4]<-5}
            if(tab.mod[j,4]<(-5)){tab.mod[j,4]<-(-5)}
            if(tab.mod[j,5]>5){tab.mod[j,5]<-5}
            if(tab.mod[j,5]<(-5)){tab.mod[j,5]<-(-5)}
        }
        for (k in 1:dim(tab.mod)[1]){
            if (tab.mod[k,4]*tab.mod[k,5]<0){tab.mod[k,6]<-(1)
                tab.mod[k,7]<-abs(tab.mod[k,5]/tab.mod[k,4])}    # if different from expectation, do futher calculation;
            if (tab.mod[k,4]*tab.mod[k,5]==0){tab.mod[k,6]<-(0)
                tab.mod[k,7]<-0}                                 # if no change, got 0 score;
            if (tab.mod[k,4]*tab.mod[k,5]>0){tab.mod[k,6]<-(-1)   # if consistant, do further calculation
                tab.mod[k,7]<-tab.mod[k,5]/tab.mod[k,4]
            }
        }
        z.score<- sum((tab.mod[,6]*tab.mod[,7]))/(dim(tab.mod)[1]^0.5)    # calculate z-score
        drug.score[i,1]<-z.score
        z.score.2<- sum(tab.mod[,6])/(dim(tab.mod)[1]^0.5) # calculate absolute z-score
        drug.score[i,2]<-z.score.2
           }
    plot(drug.score[,1],
         main="TDS repression Z-score",
         pch=1,
         xlab="Drug Index",
         ylab="Repression Z-score")
    points(drug.score[,2],
           pch=16)
    legend("topright",
           legend=c("Repression Z-score","Adjusted Repression Z-score"),
           pch=c(1,16))
})


# calculating drug score
output$plotDrugScore<-renderPlot({
    logFC<-as.numeric(input$LogFC.cutoff)   # read in logFC cut-off
    file.tds<-input$TDS_reference           # read in reference TDS file
    if(is.null(file.tds))
    return(NULL)
    data<-read.csv(file.tds$datapath, header=T,row.names=1)
    data.ref<-data[data[,1]==0,]      # For reference
    data.gene<-data[data[,1]!=0,]     # For z-score calculation
    file.drug<-input$Drug_data        # read in actual drug induced logFC (tgmel+durg/tgmel)
    if(is.null(file.drug))
    return(NULL)
    drug.table<-as.matrix(read.csv(file.drug$datapath, header=T, row.names=1))
    # sign any value between cut.off and -1*cut.off as 0
    for (i in 1:dim(drug.table)[1]){
        for (j in 1:dim(drug.table)[2]){
            if (drug.table[i,j]>0){ if (drug.table[i,j] < logFC){drug.table[i,j]<-0}}
            if (drug.table[i,j]<0){ if (drug.table[i,j] > (-1)*logFC){drug.table[i,j]<-0}}
        }
    }
    ref<-cbind(rownames(data.gene),data.gene)
    drug<-cbind(rownames(drug.table),drug.table)
    colnames(drug)[1]<-"rownames(drug.table)"
    all.tab<-as.matrix(merge(x=ref, y=drug, by.x="rownames(data.gene)", by.y="rownames(drug.table)", all=F))
    # in all.tab $logFC is the the expression change of tgmel/wt
    Tg.score<-matrix(rep(0, 4*dim(all.tab)[1]),ncol=4)   # Form disease signature logFC table and model a perfect drug. The weight value requires double check. It should be the ROC value.
    colnames(Tg.score)<-c("logFC(Tg.over.wt","weight","weight.Tg.logFC","perfect.drug.logFC")
    Tg.score[,1]<-as.matrix(all.tab[,5])
    Tg.score[,2]<-all.tab[,2]
    class(Tg.score)<-"numeric"
    Tg.score[,3]<-Tg.score[,1]*Tg.score[,2]
    Tg.score[,4]<-Tg.score[,3]*(-1)
    drug.score<-matrix(rep(0,6*(dim(all.tab)[2]-5)),ncol=6) # first 5 columns are TDS registration
    rownames(drug.score)<-colnames(all.tab)[6:dim(all.tab)[2]]
    colnames(drug.score)<-c("Zscore","adjusted.Zscore","Ks.up","Ks.dn","S","score")
    for (i in 1:dim(drug.score)[1]){
        tab<-cbind(Tg.score, all.tab[,i+5], rep(0,dim(all.tab)[1]),rep(0,dim(all.tab)[1]))
        class(tab)<-"numeric"
        colnames(tab)[5]<-colnames(all.tab)[i+5]
        colnames(tab)[6]<-"gene direction"
        colnames(tab)[7]<-"gene score"
        tab.mod<-tab                                 # set logFC limit to +/-5
        for (j in 1:dim(tab.mod)[1]){
            if(tab.mod[j,4]>5){tab.mod[j,4]<-5}
            if(tab.mod[j,4]<(-5)){tab.mod[j,4]<-(-5)}
            if(tab.mod[j,5]>5){tab.mod[j,5]<-5}
            if(tab.mod[j,5]<(-5)){tab.mod[j,5]<-(-5)}
        }
        for (k in 1:dim(tab.mod)[1]){
            if (tab.mod[k,4]*tab.mod[k,5]<0){tab.mod[k,6]<-(1)
                tab.mod[k,7]<-abs(tab.mod[k,5]/tab.mod[k,4])}    # if different from expectation, do futher calculation;
            if (tab.mod[k,4]*tab.mod[k,5]==0){tab.mod[k,6]<-(0)
                tab.mod[k,7]<-0}                                 # if no change, got 0 score;
            if (tab.mod[k,4]*tab.mod[k,5]>0){tab.mod[k,6]<-(-1)   # if consistant, do further calculation
                tab.mod[k,7]<-tab.mod[k,5]/tab.mod[k,4]
            }
        }
        z.score<- sum((tab.mod[,6]*tab.mod[,7]))/(dim(tab.mod)[1]^0.5)    # calculate z-score
        drug.score[i,1]<-z.score
        z.score.2<- sum(tab.mod[,6])/(dim(tab.mod)[1]^0.5) # calculate absolute z-score
        drug.score[i,2]<-z.score.2
    }
    # calculating Kolmogorov-Smirnov Statistics ks.up and ks.dn
    # use "data.gene", "drug.table"
    rankMatrix<-drug.table
    for (i in 1:dim(rankMatrix)[2]){
        rankMatrix[,i]<-order(rankMatrix[,i]) # convert drug logFC to a rank matrix
    }
    tds.up<-as.matrix(subset(data.gene, logFC >0))
    tds.up<-cbind(rownames(tds.up), tds.up)
    tds.dn<-as.matrix(subset(data.gene, logFC <0))
    tds.dn<-cbind(rownames(tds.dn), tds.dn)
    tag.list.up<-tds.up
    tag.list.dn<-tds.dn
    
    scores<-matrix(rep(0,dim(rankMatrix)[2]*4), ncol=4)
    colnames(scores)<-c("Ks.up","Ks.dn","S","score")
    rownames(scores)<-colnames(rankMatrix)
    
    for (i in 1:dim(rankMatrix)[2]){ # loop each incidence
        n<-dim(rankMatrix)[1]
        t.up<-dim(tag.list.up)[1]
        t.dn<-dim(tag.list.dn)[1]
        # calculating ks.up
        a.list<-rep(0, dim(tag.list.up)[1])
        b.list<-rep(0, dim(tag.list.up)[1])
        instance <- matrix(rankMatrix[,i]) # take out instance i
        rownames(instance)<-rownames(rankMatrix)
        instance.i<-as.matrix(instance[order(instance[,1]),]) # sort instance i
        
        V<-as.matrix(sort(instance.i[rownames(tag.list.up),1]))
        for (j in 1:dim(tag.list.up)[1]){
            a.list[j]<-(j/t.up)-(V[j]/n)
            b.list[j]<-(V[j]/n)-((j-1)/t.up)
        }
        a<-max(a.list)
        b<-max(b.list)
        ks.up<-0
        if(a>b){ks.up=a}
        if(a<b){ks.up=-1*b}
        scores[i,1]<-ks.up
        
        #calculating ks.dn
        a.list<-rep(0, dim(tag.list.dn)[1])
        b.list<-rep(0, dim(tag.list.dn)[1])
        instance <- matrix(rankMatrix[,i]) # take out instance i
        rownames(instance)<-rownames(rankMatrix)
        instance.i<-as.matrix(instance[order(instance[,1]),]) # sort instance i
        
        V<-as.matrix(sort(instance.i[rownames(tag.list.dn),1]))
        for (j in 1:dim(tag.list.dn)[1]){
            a.list[j]<-(j/t.dn)-(V[j]/n)
            b.list[j]<-(V[j]/n)-((j-1)/t.dn)
        }
        a<-max(a.list)
        b<-max(b.list)
        ks.dn<-0
        if(a>b){ks.dn=a}
        if(a<b){ks.dn=-1*b}
        scores[i,2]<-ks.dn
    }
    
    for (x in 1:dim(scores)[1]){
        if (scores[x,1]*scores[x,2]>0){scores[x,3]<-0}
        if (scores[x,1]*scores[x,2]<0){scores[x,3]<-scores[x,1]-scores[x,2]}
    }
    p<-max(scores[,3])
    q<-min(scores[,3])
    
    for (x in 1:dim(scores)[1]){
        if (scores[x,3]>0){scores[x,4]<-scores[x,3]/p}
        if (scores[x,3]<0){scores[x,4]<-(-1)*(scores[x,3]/q)}
    }
    #scores<-cbind(scores, instances)
    #scores<-scores[order(scores[,4]),]
    
    drug.score[,3]<-scores[,1]
    drug.score[,4]<-scores[,2]
    drug.score[,5]<-scores[,3]
    drug.score[,6]<-scores[,4]

    ######
    plot(drug.score[,6],
    main="Drug Score",
    pch=1,
    xlab="Drug Index",
    ylab="TDS score")
    ylim=c(min(drug.score[,6])-1, max(drug.score[,6])+1)
    points(drug.score[,5],
    pch=16)
    legend("topright",
    legend=c("Drug score (relative)","Drug score (S; absolute)"),
    pch=c(1,16))
})

#### Download Data
output$downloadData <- downloadHandler(
filename = function() { paste(input$outputfile, '.csv', sep='') },
content = function(file) {
    
    
    logFC<-as.numeric(input$LogFC.cutoff)   # read in logFC cut-off
    file.tds<-input$TDS_reference           # read in reference TDS file
    if(is.null(file.tds))
    return(NULL)
    data<-read.csv(file.tds$datapath, header=T,row.names=1)
    data.ref<-data[data[,1]==0,]      # For reference
    data.gene<-data[data[,1]!=0,]     # For z-score calculation
    file.drug<-input$Drug_data        # read in actual drug induced logFC (tgmel+durg/tgmel)
    if(is.null(file.drug))
    return(NULL)
    drug.table<-as.matrix(read.csv(file.drug$datapath, header=T, row.names=1))
    # sign any value between cut.off and -1*cut.off as 0
    for (i in 1:dim(drug.table)[1]){
        for (j in 1:dim(drug.table)[2]){
            if (drug.table[i,j]>0){ if (drug.table[i,j] < logFC){drug.table[i,j]<-0}}
            if (drug.table[i,j]<0){ if (drug.table[i,j] > (-1)*logFC){drug.table[i,j]<-0}}
        }
    }
    ref<-cbind(rownames(data.gene),data.gene)
    drug<-cbind(rownames(drug.table),drug.table)
    colnames(drug)[1]<-"rownames(drug.table)"
    all.tab<-as.matrix(merge(x=ref, y=drug, by.x="rownames(data.gene)", by.y="rownames(drug.table)", all=F))
    # in all.tab $logFC is the the expression change of tgmel/wt
    Tg.score<-matrix(rep(0, 4*dim(all.tab)[1]),ncol=4)   # Form disease signature logFC table and model a perfect drug. The weight value requires double check. It should be the ROC value.
    colnames(Tg.score)<-c("logFC(Tg.over.wt","weight","weight.Tg.logFC","perfect.drug.logFC")
    Tg.score[,1]<-as.matrix(all.tab[,5])
    Tg.score[,2]<-all.tab[,2]
    class(Tg.score)<-"numeric"
    Tg.score[,3]<-Tg.score[,1]*Tg.score[,2]
    Tg.score[,4]<-Tg.score[,3]*(-1)
    drug.score<-matrix(rep(0,6*(dim(all.tab)[2]-5)),ncol=6) # first 5 columns are TDS registration
    rownames(drug.score)<-colnames(all.tab)[6:dim(all.tab)[2]]
    colnames(drug.score)<-c("Zscore","adjusted.Zscore","Ks.up","Ks.dn","S","score")
    for (i in 1:dim(drug.score)[1]){
        tab<-cbind(Tg.score, all.tab[,i+5], rep(0,dim(all.tab)[1]),rep(0,dim(all.tab)[1]))
        class(tab)<-"numeric"
        colnames(tab)[5]<-colnames(all.tab)[i+5]
        colnames(tab)[6]<-"gene direction"
        colnames(tab)[7]<-"gene score"
        tab.mod<-tab                                 # set logFC limit to +/-5
        for (j in 1:dim(tab.mod)[1]){
            if(tab.mod[j,4]>5){tab.mod[j,4]<-5}
            if(tab.mod[j,4]<(-5)){tab.mod[j,4]<-(-5)}
            if(tab.mod[j,5]>5){tab.mod[j,5]<-5}
            if(tab.mod[j,5]<(-5)){tab.mod[j,5]<-(-5)}
        }
        for (k in 1:dim(tab.mod)[1]){
            if (tab.mod[k,4]*tab.mod[k,5]<0){tab.mod[k,6]<-(1)
                tab.mod[k,7]<-abs(tab.mod[k,5]/tab.mod[k,4])}    # if different from expectation, do futher calculation;
            if (tab.mod[k,4]*tab.mod[k,5]==0){tab.mod[k,6]<-(0)
                tab.mod[k,7]<-0}                                 # if no change, got 0 score;
            if (tab.mod[k,4]*tab.mod[k,5]>0){tab.mod[k,6]<-(-1)   # if consistant, do further calculation
                tab.mod[k,7]<-tab.mod[k,5]/tab.mod[k,4]
            }
        }
        z.score<- sum((tab.mod[,6]*tab.mod[,7]))/(dim(tab.mod)[1]^0.5)    # calculate z-score
        drug.score[i,1]<-z.score
        z.score.2<- sum(tab.mod[,6])/(dim(tab.mod)[1]^0.5) # calculate absolute z-score
        drug.score[i,2]<-z.score.2
    }
    # calculating Kolmogorov-Smirnov Statistics ks.up and ks.dn
    # use "data.gene", "drug.table"
    rankMatrix<-drug.table
    for (i in 1:dim(rankMatrix)[2]){
        rankMatrix[,i]<-order(rankMatrix[,i])
    }
    tds.up<-as.matrix(subset(data.gene, logFC >0))
    tds.up<-cbind(rownames(tds.up), tds.up)
    tds.dn<-as.matrix(subset(data.gene, logFC <0))
    tds.dn<-cbind(rownames(tds.dn), tds.dn)
    tag.list.up<-tds.up
    tag.list.dn<-tds.dn
    
    scores<-matrix(rep(0,dim(rankMatrix)[2]*4), ncol=4)
    colnames(scores)<-c("Ks.up","Ks.dn","S","score")
    rownames(scores)<-colnames(rankMatrix)
    
    for (i in 1:dim(rankMatrix)[2]){ # loop each incidence
        n<-dim(rankMatrix)[1]
        t.up<-dim(tag.list.up)[1]
        t.dn<-dim(tag.list.dn)[1]
        # calculating ks.up
        a.list<-rep(0, dim(tag.list.up)[1])
        b.list<-rep(0, dim(tag.list.up)[1])
        instance <- matrix(rankMatrix[,i]) # take out instance i
        rownames(instance)<-rownames(rankMatrix)
        instance.i<-as.matrix(instance[order(instance[,1]),]) # sort instance i
        
        V<-as.matrix(sort(instance.i[rownames(tag.list.up),1]))
        for (j in 1:dim(tag.list.up)[1]){
            
            a.list[j]<-(j/t.up)-(V[j]/n)
            b.list[j]<-(V[j]/n)-((j-1)/t.up)
        }
        a<-max(a.list)
        b<-max(b.list)
        ks.up<-0
        if(a>b){ks.up=a}
        if(a<b){ks.up=-1*b}
        scores[i,1]<-ks.up
        
        #calculating ks.dn
        a.list<-rep(0, dim(tag.list.dn)[1])
        b.list<-rep(0, dim(tag.list.dn)[1])
        instance <- matrix(rankMatrix[,i]) # take out instance i
        rownames(instance)<-rownames(rankMatrix)
        instance.i<-as.matrix(instance[order(instance[,1]),]) # sort instance i
        
        V<-as.matrix(sort(instance.i[rownames(tag.list.dn),1]))
        for (j in 1:dim(tag.list.dn)[1]){
            a.list[j]<-(j/t.dn)-(V[j]/n)
            b.list[j]<-(V[j]/n)-((j-1)/t.dn)
        }
        a<-max(a.list)
        b<-max(b.list)
        ks.dn<-0
        if(a>b){ks.dn=a}
        if(a<b){ks.dn=-1*b}
        scores[i,2]<-ks.dn
    }
    
    for (x in 1:dim(scores)[1]){
        if (scores[x,1]*scores[x,2]>0){scores[x,3]<-0}
        if (scores[x,1]*scores[x,2]<0){scores[x,3]<-scores[x,1]-scores[x,2]}
    }
    p<-max(scores[,3])
    q<-min(scores[,3])
    
    for (x in 1:dim(scores)[1]){
        if (scores[x,3]>0){scores[x,4]<-scores[x,3]/p}
        if (scores[x,3]<0){scores[x,4]<-(-1)*(scores[x,3]/q)}
    }
    #scores<-cbind(scores, instances)
    #scores<-scores[order(scores[,4]),]
    
    drug.score[,3]<-scores[,1]
    drug.score[,4]<-scores[,2]
    drug.score[,5]<-scores[,3]
    drug.score[,6]<-scores[,4]
    
    
    write.csv(drug.score, file)
    })


})