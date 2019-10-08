wrongInput = function(reason){
	cat(reason)
	return(NULL)
	}

bmaRelationship = function(id,from,to,type){
	rel <- sprintf('{"Id":%d,"FromVariable":%d,"ToVariable":%d,"Type":"%s"}', id, from, to, type)
	return(rel)
	}

bmaVariableModel = function(id,name,granularity,formula=""){
	var <- sprintf('{"Name":"%s","Id":%d,"RangeFrom":0,"RangeTo":%d,"Formula":"%s"}', name, id, granularity, formula)
	return(var)
	}
	
bmaVariableLayout = function(id,name,x,y,description="") {
	var <- sprintf('{"Id":%d,"Name":"%s","Type":"Constant","ContainerId":0,"PositionX":%f,"PositionY":%f,"CellX":0,"CellY":0,"Angle":0,"Description":"%s"}',id,name,x,y,description)
	return(var)
	}

bmaFormula = function(inhibitor,granularity,upstream){
	f <- ifelse(inhibitor,sprintf("%d-var(%s)",granularity,upstream),"")
	return(f)
}

bmaDescription = function(e,incoming=""){
	sign <- ifelse(e$is_stimulation==1,
        	ifelse(e$is_inhibition==1,"Mixed","Activator"),
        	ifelse(e$is_inhibition==1,"Inhibitor","Unknown")) 
	refs <- paste(unlist(e$references), sep = '', collapse = ',')
	incoming <- ifelse(incoming=="","",paste("",incoming,"",sep=" "))
	return(sprintf("%s%s-PMID:%s.",incoming,sign,refs))
	}
	
bmaMotif_es = function(edgeSeq,G,granularity=2){
    if(length(edgeSeq)==0) {
        wrongInput("\nempty path\n")
    }
    #Process- 
    ## Create list of variables
    ## Create layout of variables
    ## Create list of links
    ## Print format as follows (x is a string, xN is an integer, xF is a float)
    ### {"Model": 	{"Name": "Omnipath motif",
    ###	 	 	"Variables": [{"Name":"x","Id":xN,"RangeFrom"=0,"RangeTo"=granularity},...]
    ###    	 	"Relationships": [{"Id":xN,"FromVariable":xN,"ToVariable":xN,"Type":"Activator"},...]
    ### 		}
    ###  "Layout": 	{"Variables": [{"Id":xN,"Name":"x","Type":"Constant","ContainerId":0,"PositionX":xF,"PositionY":xF,"CellX":0,"CellY":0,"Angle":0,"Description":""},...]
    ###			"Containers":[]
    ###			}
    ### }

    #Code for identifying sign
    #signs <- ifelse(edgeSeq$is_stimulation==1,
    #    ifelse(edgeSeq$is_inhibition==1,"(+/-)","( + )"),
    #    ifelse(edgeSeq$is_inhibition==1,"( - )","( ? )"))    
    #interaction <- paste0("==", signs,"==>")
    
    #relationships <- ifelse(edgeSeq$is_stimulation==1,
    #    ifelse(edgeSeq$is_inhibition==1,"Activator","Activator"),
    #    ifelse(edgeSeq$is_inhibition==1,"Inhibitor",return(wrongInput("\nUnsigned input graph\n")))    
    sources <- tail_of(G, edgeSeq)$name
    variableNames <- c(sources,head_of(G, edgeSeq)$name[length(edgeSeq)]) 
    varNum <- length(variableNames)
    
    positions <- vector('list',varNum)
    variables <- vector('list',varNum)
    relationships <- vector('list',varNum-1)
    
    
    formula = ""
    description = ""
    x=125
    y=140
    for (i in seq_along(variableNames))
    	{
	v <- bmaVariableModel(i,variableNames[i],granularity,formula)
	p <- bmaVariableLayout(i,variableNames[i],x,y,description)
	if (i < varNum){
		#Simplified sign- if inhibition, inhibitor, else (activator/mixed/unknown) activation
		r <- bmaRelationship(i+varNum,i,i+1,ifelse(edgeSeq[i]$is_inhibition==1,"Inhibitor","Activator"))
		relationships[[i]] <- r
		formula <- bmaFormula((edgeSeq[i]$is_inhibition==1),granularity,variableNames[i])
		description <- bmaDescription(edgeSeq[i])
		}
	positions[[i]] <- p
	variables[[i]] <- v

	x <- x + 86.6025404
	ymod <- ifelse(i%%2==0,50,-50)
	y <- y + ymod
	}
   result <- sprintf('{"Model": {"Name": "Omnipath motif","Variables":[%s],"Relationships":[%s]},"Layout":{"Variables":[%s],"Containers":[]}}\n',paste(variables, sep = '', collapse = ','),paste(relationships, sep = '', collapse = ','),paste(positions, sep = '', collapse = ','))
   cat(result)
    
}

bmaMotif_vs = function(nodeSeq,G){

    if(length(nodeSeq)==0){
        print("empty path")
        return(invisible(NULL))
    }
    nodeSeq_names <- unique_nodeSeq(nodeSeq)
    for(i in seq(nodeSeq_names)){
        print(paste0("pathway ", i, ": ", 
            paste(nodeSeq_names[[i]],collapse = " -> ")))
        edgeSet <- c()
        for(j in 2:length(nodeSeq_names[[i]])){
            edgeSet <- c(edgeSet, E(G)[nodeSeq_names[[i]][[j-1]]  %->% 
                nodeSeq_names[[i]][[j]]])
        }
        bmaMotif_es(E(G)[edgeSet],G)
    }
}
