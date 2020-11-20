#### LIBRARIES AND SOURCE CODE #####

# LIBRARIES
#library(tidyverse)
library(xml2)
library(igraph)
library(RCy3)

# SOURCE CODE
source("AOP-Net-Functions.R")

#### FILES #####
aopXML<-"aop-wiki-xml-2020-10-01"
AOPs_of_interest<-unlist(read.table("thyroid-AOPN-AOP-IDs.txt", header = TRUE, sep="\t"), use.names = FALSE)

#### LOAD XML AND REF ID DATA ####

# IMPORT XML
xData<-read_xml(aopXML)
xData<-xml_ns_strip(xData)

# REF ID TO AOPWIKI ID
keID<-data.frame(
  ref=xml_attr(xml_find_all(xData, "/data/vendor-specific/key-event-reference"),"id"),
  ID=as.numeric(xml_attr(xml_find_all(xData, "/data/vendor-specific/key-event-reference"),"aop-wiki-id")),
  stringsAsFactors=FALSE
)

kerID<-data.frame(
  ref=xml_attr(xml_find_all(xData, "/data/vendor-specific/key-event-relationship-reference"),"id"),
  ID=as.numeric(xml_attr(xml_find_all(xData, "/data/vendor-specific/key-event-relationship-reference"),"aop-wiki-id")),
  stringsAsFactors=FALSE
)

aopID<-data.frame(
  ref=xml_attr(xml_find_all(xData, "/data/vendor-specific/aop-reference"),"id"),
  ID=as.numeric(xml_attr(xml_find_all(xData, "/data/vendor-specific/aop-reference"),"aop-wiki-id")),
  stringsAsFactors=FALSE
)

taxID<-data.frame(
  ref=xml_attr(xml_find_all(xData, "/data/vendor-specific/taxonomy-reference"),"id"),
  ID=as.numeric(xml_attr(xml_find_all(xData, "/data/vendor-specific/taxonomy-reference"),"aop-wiki-id")),
  stringsAsFactors=FALSE
)

#### TAXONOMY DATA TABLE ####
taxData<-data.frame(
  ID=taxID$ID[match(xml_attr(xml_find_all(xData, "/data/taxonomy"), "id"),taxID$ref)],
  name=xml_text(xml_find_all(xData, "/data/taxonomy/name")),
  source=xml_text(xml_find_all(xData, "/data/taxonomy/source")),
  source_ID=xml_text(xml_find_all(xData, "/data/taxonomy/source-id")),
  stringsAsFactors=FALSE
)


#### KEY EVENT (KE) DATA TABLE ####

# KE taxonomy data
keTaxDom<-lapply(xml_find_all(xData, "/data/key-event/applicability"),FUN=function(x){
  if("taxonomy"%in%xml_name(xml_children(x))){
    tid<-taxID$ID[match(xml_attr(xml_find_all(x, "taxonomy"),"taxonomy-id"),taxID$ref)]
    return(taxData$name[taxData$ID%in%tid])
  }
})

# proportion of KEs with taxonomic information
sum(!sapply(keTaxDom, is.null)) / length(keTaxDom) # 247 / 1118 = 0.2209302 :(

# Other KE data
keData<-data.frame(
  ID=keID$ID[match(xml_attr(xml_find_all(xData, "/data/key-event"), "id"),keID$ref)],
  title=xml_text(xml_find_all(xData, "/data/key-event/title")),
  LOBO=xml_text(xml_find_all(xData, "/data/key-event/biological-organization-level")),
  taxDom=I(keTaxDom),
  stringsAsFactors=FALSE
)


#### KEY EVENT RELATIONSHIP (KER) DATA TABLE ####

# KER taxonomy data
kerTaxDom<-lapply(xml_find_all(xData, "/data/key-event-relationship/taxonomic-applicability"),FUN=function(x){
  if("taxonomy"%in%xml_name(xml_children(x))){
    tid<-taxID$ID[match(xml_attr(xml_find_all(x, "taxonomy"),"taxonomy-id"),taxID$ref)]
    return(taxData$name[taxData$ID%in%tid])
  }
})

# proportion of KERs with taxonomic information
sum(!sapply(kerTaxDom, is.null)) / length(kerTaxDom) # 230 / 1338 = 0.1718984

# all other KER Data
kerData<-data.frame(
  ID=kerID$ID[match(xml_attr(xml_find_all(xData, "/data/key-event-relationship"), "id"),kerID$ref)],
  KEup=keID$ID[match(xml_text(xml_find_all(xData, "/data/key-event-relationship/title/upstream-id")),keID$ref)],
  KEdown=keID$ID[match(xml_text(xml_find_all(xData, "/data/key-event-relationship/title/downstream-id")),keID$ref)],
  taxDom=I(kerTaxDom),
  stringsAsFactors=FALSE
)

#### AOP DATA TABLE ####

# OECD status: not all aops have an "oecd-status" xml tag, so must us "if" to return NA when missing
oecdStatus<-sapply(xml_find_all(xData, "/data/aop/status"),FUN=function(x){
  if("oecd-status"%in%xml_name(xml_children(x))){
    return(xml_text(xml_find_all(x,"oecd-status")))
  }else{
    return("not specified")
  }
})

# SAAOP status: not all aops have an "saaop-status" xml tag, so must us "if" to return NA when missing
saaopStatus<-sapply(xml_find_all(xData, "/data/aop/status"),FUN=function(x){
  if("saaop-status"%in%xml_name(xml_children(x))){
    return(xml_text(xml_find_all(x,"saaop-status")))
  }else{
    return("not specified")
  }
})

# MIEs: more than one MIE possible per aop, so must return list
mies<-lapply(xml_find_all(xData, "/data/aop"),FUN=function(x){
  if("molecular-initiating-event"%in%xml_name(xml_children(x))){
    return(keID$ID[match(xml_attr(xml_find_all(x, "molecular-initiating-event"),"key-event-id"),keID$ref)])
  }else{
    return(NULL)
  }
})

# AOs: more than one AO possible per aop, so must return list
aos<-lapply(xml_find_all(xData, "/data/aop"),FUN=function(x){
  if("adverse-outcome"%in%xml_name(xml_children(x))){
    return(keID$ID[match(xml_attr(xml_find_all(x, "adverse-outcome"),"key-event-id"),keID$ref)])
  }else{
    return(NULL)
  }
})

# KEs: more than one KE possible per aop, so must return list
kes<-lapply(xml_find_all(xData, "/data/aop/key-events"),FUN=function(x){
  if("key-event"%in%xml_name(xml_children(x))){
    return(keID$ID[match(xml_attr(xml_find_all(x, "key-event"),"id"),keID$ref)])
  }else{
    return(NULL)
  }
})

# KERs: more than one KER per aop, each with aop-specific "adjaceny", "quantitative understanding", and "WoE"
# So must return data frame of KERs
kers<-lapply(xml_find_all(xData, "/data/aop/key-event-relationships"),FUN=function(x){
  if("relationship"%in%xml_name(xml_children(x))){
    return(data.frame(
      ID=kerID$ID[match(xml_attr(xml_find_all(x, "relationship"),"id"),kerID$ref)],
      adjacency=xml_text(xml_find_all(x, "relationship/adjacency")),
      quant=xml_text(xml_find_all(x, "relationship/quantitative-understanding-value")),
      woe=xml_text(xml_find_all(x, "relationship/evidence")),
      stringsAsFactors=FALSE
    ))
  }else{
    return(NULL)
  }
})

# add kes and MIE/AO designation (which is AOP-specific) for each KER in kers data.frame
for(i in 1:length(kers)){
  if(length(kers[[i]])>0){
    KEup<-kerData$KEup[match(kers[[i]]$ID,kerData$ID)]
    KEDup<-sapply(KEup, FUN=function(x){
      if(x%in%mies[[i]]){
        return("MIE")
      }else{
        if(x%in%aos[[i]]){
          return("AO")
        }else{
          return("KE")
        }
      }
    })
    
    KEdown<-kerData$KEdown[match(kers[[i]]$ID,kerData$ID)]
    KEDdown<-sapply(KEdown, FUN=function(x){
      if(x%in%mies[[i]]){
        return("MIE")
      }else{
        if(x%in%aos[[i]]){
          return("AO")
        }else{
          return("KE")
        }
      }
    })
    
    kers[[i]]<-data.frame(
      ID=kers[[i]]$ID,
      KEup=KEup,
      KEDup=KEDup,
      KEdown=KEdown,
      KEDdown=KEDdown,
      adjacency=kers[[i]]$adjacency,
      quant=kers[[i]]$quant,
      woe=kers[[i]]$woe,
      row.names=NULL,
      stringsAsFactors = FALSE
    )
  }
}

aopData<-data.frame(
  ID=aopID$ID[match(xml_attr(xml_find_all(xData, "/data/aop"), "id"),aopID$ref)],
  oecdStatus=oecdStatus,
  saaopStatus=saaopStatus,
  mies=I(mies),
  aos=I(aos),
  kes=I(kes),
  kers=I(kers),
  stringsAsFactors=FALSE
)


#### CREATE AOP NETWORK (iGraph Object) ####

# remove all "archived" aops based on SAAOP status (and not NA)
aData<-aopData[!aopData$saaopStatus=="Archived",]

# remove all "empty" aops (i.e. aops that have ZERO mies, aos, kes, and/or ZERO kers)
mieNull<-sapply(aData$mies, is.null)
aoNull<-sapply(aData$aos, is.null)
keNull<-sapply(aData$kes, is.null)
kerNull<-sapply(aData$kers, is.null)
emptyAOP<-(mieNull & aoNull & keNull) | kerNull
aData<-aData[!emptyAOP,]

# Combine all KER tables from all AOPs
kerSet<-do.call("rbind", aData$kers) 

# Create edge list from unique "KEup" and "KEdown" pairs ("as.character" ensures the Nodes are named correctly when converting to iGraph object)
edgeList<-unique(kerSet[,c("ID","KEup","KEdown")])
edgeList$ID<-as.character(edgeList$ID)
edgeList$KEup<-as.character(edgeList$KEup)
edgeList$KEdown<-as.character(edgeList$KEdown)


# Create iGraph object from edgeList (edgelist must be a matrix)
g1<-graph_from_edgelist(as.matrix(edgeList[,c("KEup","KEdown")]), directed=TRUE)


#### MAP VERTEX (KE) ATTRIBUTES ####

# add "ID" attribute 
V(g1)$ID<-V(g1)$name

# Add KE_KED (KED=Key event designation. "MIE", "KE", or "AO") based on "once an MIE/AO always an MIE/AO" (3 steps)
# step 1:combine ke/ked data from each aop into one table
kedDat<-do.call("rbind", aData$kers) 

# step 2: combine all unique KE/ KEDs  into a single two column table
kedDat<-data.frame(KE=c(kedDat$KEup,kedDat$KEdown), KED=c(kedDat$KEDup, kedDat$KEDdown), stringsAsFactors=FALSE) 
kedDat<-unique(kedDat[order(as.numeric(kedDat$KE)),])

# step 3: assign KED to V(g1) based on table ("once an MIE/AO always an MIE/AO")
V(g1)$KE_KED<-sapply(V(g1)$ID, FUN=function(x){
  if("MIE"%in%kedDat$KED[kedDat$KE==x]){
    return("MIE")   
  }else{
    if("AO"%in%kedDat$KED[kedDat$KE==x]){
      return("AO")
    }else{
      return("KE")
    }
  }
})


# Add KE_PD (PD="Path Designation. path origin or path terminus)
g1<-add_KE_PD(g1)


# add colour attribute (MIE=green, KE=white, AOP=red)
V(g1)$col<-"white"
V(g1)$col[V(g1)$KE_KED=="MIE"]<-"green"
V(g1)$col[V(g1)$KE_KED=="AO"]<-"red"


# Add AOP IDs. Note: Each KE may belong to more than one AOP, so the AOP_IDs object is a list not just a single AOP ID
V(g1)$AOP_ID<-list(vector())
for(i in 1:nrow(aData)){
  keAll<-c(aData$mies[[i]], aData$aos[[i]], aData$kes[[i]])
  if(length(keAll)>0){
    for(j in keAll){
      if(j%in%V(g1)$ID){
        V(g1)$AOP_ID[[match(j,V(g1)$ID)]]<-c(V(g1)$AOP_ID[[match(j,V(g1)$ID)]],aData$ID[i])
      }
    }
  }
}

# map  level of biological organization (LOBO)
V(g1)$LOBO<-keData$LOBO[match(V(g1)$ID, keData$ID)]

# map titles as KE attribute
V(g1)$title<-keData$title[match(V(g1)$ID, keData$ID)]


# map taxonomic domains 
V(g1)$taxDom<-keData$taxDom[match(V(g1)$ID, keData$ID)]



#### MAP EDGE (KER) ATRTIBUTeS ####

# add KER ID as edge attribute
E(g1)$ID<-edgeList$ID

# Add AOP IDs. Note: Each KER may belong to more than one AOP, so the AOP_IDs object is a list not just a single AOP ID
E(g1)$AOP_ID<-list(vector())
for(i in 1:nrow(aData)){
  if(!is.null(aData$kers[[i]])){
    for(j in aData$kers[[i]]$ID){
      if(j%in%E(g1)$ID){
        E(g1)$AOP_ID[[match(j,E(g1)$ID)]]<-c(E(g1)$AOP_ID[[match(j,E(g1)$ID)]],aData$ID[i])
      }
    }
  }
}

# WoE and quatitative understanding

# currently WOE and Quant are AOP-SPECIFIC
# For this study, we will use the LOWEST WOE and Quant assigned to a KER, if it has multiple values

#woe
woeList<-unique(kerSet[,c("ID","KEup","KEdown", "woe")])
dupIDs<-unique(woeList$ID[duplicated(woeList$ID)])
w<-vector()
for(i in dupIDs){
  kerW<-woeList$woe[woeList$ID==i]
  if("Low"%in%kerW){
    w<-c(w, "Low")
  }else{
    if("Moderate"%in%kerW){
      w<-c(w, "Moderate")
    }else{
      if("High"%in%kerW){
        w<-c(w,"High")
      }else{
        w<-c(w,"Not Specified")
      }
    }
  }
}
# remove duplicates
woeList<-woeList[!duplicated(woeList$ID),]
# add woes for duplicates
woeList$woe[match(dupIDs, woeList$ID)]<-w
# assign as edge attribute
E(g1)$woe[match(woeList$ID, E(g1)$ID)]<-woeList$woe

#convert to numeric score
wScores<-data.frame(w=c("High","Moderate","Low","Not Specified"), score=c(3, 2, 1, 1))
E(g1)$woe_score<-wScores$score[match(E(g1)$woe, wScores$w)]


# quant
quantList<-unique(kerSet[,c("ID","KEup","KEdown", "quant")])
dupIDs<-unique(quantList$ID[duplicated(quantList$ID)])
q<-vector()
for(i in dupIDs){
  kerQ<-quantList$quant[quantList$ID==i]
  if("Low"%in%kerQ){
    q<-c(q, "Low")
  }else{
    if("Moderate"%in%kerQ){
      q<-c(q, "Moderate")
    }else{
      if("High"%in%kerQ){
        q<-c(q,"High")
      }else{
        q<-c(q,"Not Specified")
      }
    }
  }
}
# remove duplicates
quantList<-quantList[!duplicated(quantList$ID),]
# add quants for duplicates
quantList$quant[match(dupIDs, quantList$ID)]<-q
# assign as edge attribute
E(g1)$quant[match(quantList$ID, E(g1)$ID)]<-quantList$quant

#convert to numeric score
qScores<-data.frame(w=c("High","Moderate","Low","Not Specified"), score=c(3, 2, 1, 1))
E(g1)$quant_score<-wScores$score[match(E(g1)$quant, wScores$w)]

# add taxanomic domain
E(g1)$taxDom<-kerData$taxDom[match(E(g1)$ID, kerData$ID)]



#### SUBGRAPG OF AOPS OF INTEREST ####

select_KEs<-which(sapply(V(g1)$AOP_ID, function(x){
  any(x%in%unlist(AOPs_of_interest))
}))
aopNet<-induced_subgraph(g1, select_KEs)

# identify Adjacent and Non-Adjacent KERs #
aopNet<-add_KER_adjacency(aopNet)


#only endorsed AOPs
oecdAOPs<-AOPs_of_interest[aopData$oecdStatus[match(AOPs_of_interest, aopData$ID)]=="TFHA/WNT Endorsed"]

select_KEs<-which(sapply(V(g1)$AOP_ID, function(x){
  any(x%in%unlist(oecdAOPs))
}))
oecdNet<-induced_subgraph(g1, select_KEs)
oecdNet<-add_KER_adjacency(oecdNet)

g<-aopNet
#g<-oecdNet
#g<-g1


#### SPECIES FITLERS/LAYERS #####

allSpeciesKE<-unique(unlist(V(g)$taxDom))
allSpeciesKER<-unique(unlist(E(g)$taxDom))
allSpecies<-base::union(allSpeciesKE, allSpeciesKER)
allSpecies<-allSpecies[order(allSpecies)]
# problem #1: there are mores species defined in the KEs than the KERs
# problem #2: different "names" for the same species (e.g. "human" vs "humans" vs "homo sapien")

speciesKEs<-lapply(allSpecies, function(x){
  sapply(V(g)$taxDom, function(y){
    x%in%unlist(y)
  })
})
names(speciesKEs)<-allSpecies

speciesKERs<-lapply(allSpecies, function(x){
  sapply(E(g)$taxDom, function(y){
    x%in%unlist(y)
  })
})
names(speciesKERs)<-allSpecies


# species tag for filters and layers
V(g)$v_tag<-speciesKEs$human | speciesKEs$humans
E(g)$e_tag<-speciesKERs$human | speciesKERs$humans



#### PLOT GRAPH ####

# standard plot variables
set.seed(1)
l<-layout_with_fr(g, grid="nogrid")

eCol<-rep("grey50", length(E(g)))
eCol[E(g)$adjacency=="non-adjacent"]<-"orange"

eWidth<-rep(2.5,length(E(g)))
eWidth[E(g)$adjacency=="non-adjacent"]<-1.5

eType<-rep(1,length(E(g)))
eType[E(g)$adjacency=="non-adjacent"]<-3


# Plot "filter layer"
v_layer<-rep(rgb(1,1,1,alpha=0),length(V(g)))
v_layer[V(g)$v_tag]<-hsv(0.59, 0.85, 0.95)

e_layer<-rep(rgb(1,1,1,alpha=0),length(E(g)))
e_layer[E(g)$e_tag]<-hsv(0.59, 0.85, 0.95)

par(mar=c(0,0,0,0))
plot(g,
     # Vertices
     vertex.size=8, vertex.color=v_layer, vertex.frame.color=rgb(1,1,1,alpha=0), vertex.label=NA,
     #edges
     edge.width=5, edge.color=e_layer, edge.arrow.size=0.15, edge.arrow.width=3,
     #layout
     layout=l)

# plot standard layer
plot(g,
     # Vertices
     vertex.size=5.5, vertex.color=V(g)$col,
     vertex.label.cex=0.6,
     #vertex.label=V(g)$title, 
     vertex.label=NA, 
     #edges
     edge.width=eWidth, edge.color=eCol, edge.arrow.size=0.15, edge.arrow.width=3, edge.lty=eType,
     #layout
     layout=l,
     add=TRUE)


#### EXPORT IGRAPH OBJECT TO CYTOSCAPE ####
createNetworkFromIgraph(g, "thyroidNet_human_tags")

