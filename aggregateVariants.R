aggregateVariants <- function(variantAnnot, aggregateIf = NULL, aggregateBy = "Gene"){

	if(!is.null(aggregateIf)){

		condVar <- rep(NA, length(aggregateIf))
		condQual <- rep(NA, length(aggregateIf))
		condVal <- rep(NA, length(aggregateIf))

		# parse variable conditions
		for(i in 1:length(aggregateIf)){			
			varval <- unlist(strsplit(aggregateIf[i], split = "[[:space:]]*([!=]=)[[:space:]]*"))
			if(!(varval[1] %in% names(variantAnnot))){ stop(paste(varval[1], "is not a variable in variantAnnot")) }
			condVar[i] <- varval[1]
			condVal[i] <- varval[2]
			qual <- regexpr("[!=]=", aggregateIf[i])
			condQual[i] <- substr(aggregateIf[i], start=qual[1], stop=qual[1]+attr(qual, "match.length")-1)			
		}

		# set up text strings specifying conditons for inclusion
		conds <- NULL
		for(var in unique(condVar)){
			qual <- condQual[condVar == var]
			val <- condVal[condVar == var]

			if(all(qual == "==")){
				conds <- append(conds, paste0("(", paste0("variantAnnot$", var, qual, '"', val, '"', collapse = " | "), ")"))
			}else if(all(qual == "!=")){
				conds <- append(conds, paste0("(", paste0("variantAnnot$", var, qual, '"', val, '"', collapse = " & "), ")"))
			}
		}
		# get logical vector for inclusion
		condcol <- parse(text = paste(conds, collapse = " & "))	

		# subset
		variantAnnot <- variantAnnot[eval(condcol), ]

		message(paste("Selected", nrow(variantAnnot), "variants with", condcol))
	}

	# find variant groups
	variantGroups <- by(variantAnnot, variantAnnot[,aggregateBy], list)

	# re-order to match original ordering
	myorder <- match(unique(variantAnnot[,aggregateBy]), names(variantGroups))
	variantGroups <- variantGroups[myorder]
	
	message(paste("Aggregated variants into", length(variantGroups), "groups by", aggregateBy))

	return(variantGroups)
}
