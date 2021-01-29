"""
Taken in colptr array and a column index j, return an array representing the range of indices in rowval.
"""
function cran(colptr::Array{Int64},j::Int64)
	return (colptr[j]:(colptr[j+1]-1))::UnitRange{Int64}
end

function Mrv_Mcp_Mvl_Mm_rows_cols__slicetranspose(Arv::Array{Ti},Acp::Array{Ti},Avl::Array{Tv},Am;rows = 1:Am,cols = 1:length(Acp)-1) where {Tv, Ti}
	if rows == 1:Am && cols == 1:(length(Acp)-1)
		Crv, Ccp, Cvl = sparsetranspose(Arv::Array{Ti},Acp::Array{Ti},Avl,Am)
		return Crv, Ccp, Cvl
	end
    # Attach destination matrix
    Cm = length(cols)
    Cn = length(rows)
    Ccp = zeros(Ti,Cn+1) # we'll fill this in starting from zero
    # Compute the column counts of C and store them shifted forward by one in Ccp
	rs = rowsupportsum(Arv,Acp,Am,cols)
	for i = 1:Cn
	    Ccp[i+1] = rs[rows[i]]
	end
	Cnnz = sum(Ccp)
    Crv = Array{Ti}(undef,Cnnz)
    Cvl = Array{Tv}(undef,Cnnz)
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccp
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccp[k]
        Ccp[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crv and Cvl,
    # tracking write positions in Ccp
    rowtranslator = zeros(Int64,Am)
    rowtranslator[rows] = 1:Cn
    @inbounds for j in 1:length(cols) 		#   for each column to be transposed from A
        for Ak in cran(Acp, cols[j])        #   for each nonzero row of this column
            Ai = rowtranslator[Arv[Ak]]     #   find the corresponding column of the new matrix
            x  = Avl[Ak]                    #   and the corrsesponding scalar
            if Ai > 0
	            Ck = Ccp[Ai+1]              #   decide where to write
	            Crv[Ck] = j                 #   write the row value
                Cvl[Ck] = x                 #   write the scalar value
    	        Ccp[Ai+1] += 1              #   increment the pointer one place
    	    end
        end
    end
    # Tracking write positions in Ccp as in the last block fixes the cp shift,
    # but the first cp remains incorrect
    Ccp[1] = 1
	return Crv, Ccp, Cvl
end


function sparsetranspose(Arv::Array{Ti},Acp::Array{Ti},Avl::Array{Tv},Am::Integer) where {Tv, Ti}
    Annz = Acp[end]-1
    An = length(Acp)-1
    Cm = An
    Cn = Am
    Ccp = zeros(Ti,Am+1) # we'll fill this in starting from zero
    Crv = Array{Ti}(undef,Annz)
    Cvl = Array{Tv}(undef,Annz)
    # Compute the column counts of C and store them shifted forward by one in Ccp
    @inbounds for k in 1:Annz
        Ccp[Arv[k]+1] += 1
    end
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccp
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccp[k]
        Ccp[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crv and Cvl,
    # tracking write positions in Ccp
    @inbounds for Aj in 1:An
        for Ak in Acp[Aj]:(Acp[Aj+1]-1)
            Ai = Arv[Ak]
            Ck = Ccp[Ai+1]
            Crv[Ck] = Aj
            Cvl[Ck] = Avl[Ak]
            Ccp[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccp as in the last block fixes the cp shift,
    # but the first cp remains incorrect
    Ccp[1] = 1
	return Crv, Ccp, Cvl
end

function Mrv_Mcp_Mvl_Mm__reducedform(		Mrv,Mcp,Mvl,Mm;
											returntransform=false,
											returntempreduced=false,
											returnpivotindices=false)
	""" THINGS TO DO
		- EXPERIMENT WHETHER BETTER TO MAKE THE DYNAMIC ARRAYS (SPECIFICALLY PSTA + ACP) BIG ENOUGH THAT THEY NEVER HAVE TO GROW
	"""


	Qm 							=	Mm
	Qn 							=	length(Mcp)-1
	Qrv 						=	Mrv
	Qcp 						=	Mcp
	Qvl 						=	Mvl

	# 	REDUCED MATRIX
	Rrv 						=	Array{Int64}(undef,length(Mrv))
	Rcp 						=	Array{Int64}(undef,length(Mcp))
	Rcp[1] 						=	1
	Rvl 						=	Array{Rational{Int64}}(undef,length(Mrv))
	Rm  						=	Qm

	# 	ACCUMULATOR FOR REDUCED MATRIX
	rrv 						=	Array{Int64}(undef,length(Mcp))
	rcp 						=	Array{Int64}(undef,length(Mcp))
	rcp[1] 						=	1
	rvl 						=	Array{Rational{Int64}}(undef,length(Mcp))
	rpa 						=	Array{Int64}(undef,length(Mcp))
	rpca 						=	[0]

	if returntransform
		# 	TRANSFORM
		Trv 						=	Array{Int64}(undef,length(Mcp))
		Tcp 						=	Array{Int64}(undef,length(Mcp))
		Tcp[1] 						=	1
		Tvl 						=	Array{Rational{Int64}}(undef,length(Mcp))

		# 	ACCUMULATOR FOR TRANSFORM
		trv 						=	Array{Int64}(undef,length(Mcp))
		tcp 						=	Array{Int64}(undef,length(Mcp))
		tcp[1] 						=	1
		tvl 						=	Array{Rational{Int64}}(undef,length(Mcp))
		tpa 						=	Array{Int64}(undef,length(Mcp))
		tpca 						=	[0]

	end

	# 	PAIR TRACKER
	row__pcola::Array{Int64}	=	zeros(Int64,Qm) # we will use this vector to store the indices of paired coluMms
    for colind = 1:Qn

		# 	IF THE COLUMN IS STRUCTURALLY EMPTY:
		# 		1 	APPEND TO THE REDUCED MATRIX
		# 		3 	CONTINUE
		if Qcp[colind] 			== 	Qcp[colind+1]
			Rcp[colind+1] 		=	Rcp[colind]
			if returntransform
				Tcp[colind+1] 	=	Tcp[colind]
			end
			continue
		end

		#	STATE: THE COLUMN IS STRUCTURALLY NONZERO.
		# 	--------------------------------------------------------------------
		# 	LOCATE THE STRUCTURAL BOTTOM ELEMENT
		rv_locus_max 			=	Qcp[colind+1]-1
		# 	WALK BACKWARD THROUGH ANY ZERO ENTRIES
		#	IF ALL ENTRIES ARE ZERO:
		#		1	APPEND A ZERO COLUMN TO THE REDUCED MATRIX
		#		2	CONTINUE, WE ARE DONE WITH THIS COLUMN

		while Qvl[rv_locus_max] == 0
			rv_locus_max		-=	1
			if rv_locus_max 	<	Qcp[colind]
				break
			end
		end


		if rv_locus_max 	<	Qcp[colind]
			Rcp[colind+1] 			=	Rcp[colind]
			if returntransform
				Tcp[colind+1] 		=	Tcp[colind]
			end
			continue
		end

		#	STATE: THE BOTTOM ENTRY IS NONZERO
		# 	--------------------------------------------------------------------
		#	IF THE (NONZERO) BOTTOM ENTRY IS UNPAIRED
		# 		1 	MAKE A SPECIAL (NEGATIVE) MARK ON THE pairingtracker
		# 		2 	APPEND A ZERO COLUMN TO N
		# 		3 	CONTINUE

		if row__pcola[Qrv[rv_locus_max]] 	== 0
			row__pcola[Qrv[rv_locus_max]] 	=	-colind
			Rcp[colind+1] 					=	Rcp[colind]
			if returntransform
				Tcp[colind+1] 				=	Tcp[colind]
			end
			continue
		end
		#	STATE: THE BOTTOM ELEMENT IS NONZERO AND UNPAIRED
		#-----------------------------------------------------------------------
		# 		1	WRITE THE CURRENT COLUMN TO THE REDUCTION ACCUMULATOR MATRIX
		# 		2 	RESET THE TRANSFORM ACCUMULATOR TO ZERO
		# 		3 	INITIATE THE WHILE-LOOP

		#===========================================================================================================================#

		fromrange 						=	cran(Qcp,colind)
		rrv_rcp_rvl_rpa_rpca_rv_vl__reset!(rrv,rcp,rvl,rpa,rpca, Qrv[fromrange], Qvl[fromrange])
		if returntransform
			tpca[1] 					=	0
		end

		# 	STATE: THE COLUMN WITH THIS INDEX IN THE TRANSFORM MATRIX WILL BE NONZERO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#-----------------------------------------------------------------------
		# 	AFTER THIS WHILE-LOOP, WE WILL AGGREGATE THE CHANGES TO THE TRANSFORM
		# 	MATRIX, AND RECORD THEM TO THE MASTER COPY

		while true
			#LOCATE THE STRUCTURAL BOTTOM ELEMENT
			# lowest nonzero row???
			locus_max::Int64		=	rrv_rcp_rpa_rpca__locus_max(rrv,rcp,rpa,rpca)

			# 	IF THE COLUMN IS STRUCTURALLY ZERO
			# 		1	APPEND TO THE OUTPUT MATRIX
			# 		2 	BEAK THE WHILE-LOOP
			# 		3 	IF RECORDING THE TRANSFORM, IT'S COLUMN-ACCUMULATOR WILL BE AGGREGATED AND ADDED TO THE MASTER COPY AFTER THE WHILE-LOOP
			if locus_max 			==	0
				Rcp[colind+1] 		=	Rcp[colind]
				break
			end

			#	STATE: THERE IS A STRUCTURALLY NONZERO BOTTOM ELEMENT
			#-------------------------------------------------------------------

			# 	CALCULATE ITS VALUE
			entry_low 				=	rrv_rcp_rvl_rpa_rpc_locus__rowsum(rrv,rcp,rvl,rpa,rpca[1],locus_max;moveposts=true) # MOVES THE POSTS BACK!!!!!

			# 	IF THE BOTTOM STORED VALUE IS ZERO
			# 		1	CONTINUE TO THE NEXT ITERAITON OF THE WHILE LOOP
			if entry_low == 0
				continue
			end

			# 	STATE: THE BOTTOM STORED VALUE IS NONZERO
			#-------------------------------------------------------------------
			# 	IF THE ROW IS UNPAIRED
			# 		1	COMPLETE ALL INITIALIZED REDUCTIONS
			# 		2	ADD TO THE REDUCED MATRIX
			# 		3 	MARK THE PAIR TRACKER
			# 		4 	BREAK THE WHILE-LOOP; WE ARE DONE WITH THIS COLUMN

			if row__pcola[locus_max] == 0

				# 	ACCUMULATE + ADD TO THE REDUCED MATRIX
				rv,vl 			=	rrv_rcp_rvl_rpa_rpca__colsumIE_rv_vl!(rrv,rcp,rvl,rpa,rpca)
				Srv_Scp_Svl_colind_rv_vl__addwithextra!(Rrv,Rcp,Rvl,colind,rv,vl)

				# 	THIS OPERATION HAS FORGOTTEN ABOUT THE LOW ENTRY ITSELF
				# 	SPECIFICALLY, B/C WE MOVED THE POSTS BACK.  THEREFORE WE HAVE TO ADJUST
				# 	NOTE THAT THIS HAPPENS AT MOST ONCE PER COLUMN
				if length(Rrv) < Rcp[colind+1]
					append!(Rrv,Array{Int64}(undef,Qn))  # MAKE SURE THERE IS SPACE FOR THE NEW ENTRY
					append!(Rvl,Array{Rational{Int64}}(undef,Qn))  # MAKE SURE THERE IS SPACE FOR THE NEW ENTRY
				end
				Rrv[Rcp[colind+1]] 	= 	locus_max
				Rvl[Rcp[colind+1]] 	=	entry_low	# this is what we had to make space for
				Rcp[colind+1]  		+=	1

				# 	MARK THIS COLUMN AS THE PAIR TO ROW locus_max
				row__pcola[locus_max] 	=	colind

				break # WE WILL TEND TO THE TRANSFORM MATRIX AFTER THE WHILE-LOOP
			else

				# 	STATE: A STRICTLY LEFTWARD COLUMN PAIRS WITH THE BOTTOM ROW
				# 	--------------------------------------------------------------------
				# 		1 	USE THIS COLUMN TO CLEAR THE BOTTOM ENTRY
				# 		NOTE: WE HAVE ALREADY MOVED THE EXISTING POSTS, WHEN THE FOLLOWING PROCESS BEGINS

				# 	GET THE INDEX OF THE COLUMN USED FOR CLEARING
				addcolind					=	abs(row__pcola[locus_max]) 				# RECALL THAT NEGATIVES INDICATE UNCHANGED COLUMNS

				# 	DETERMINE THE SCALE PARAMETER + THE MATRIX THAT WILL CONTRIBUTE THE COLUMN
				if row__pcola[locus_max] > 0
					Wrv::Array{Int64},Wcp::Array{Int64},Wvl::Array{Rational{Int64}} 	=	Rrv,Rcp,Rvl
				else
					Wrv,Wcp,Wvl 														=	Qrv,Qcp,Qvl
				end
				scalar 						=	-entry_low/Wvl[Wcp[addcolind+1]-1] # = (-entrytobecleare)/(bottomentryofclearingcolumn)

				# 	THERE IS NOTHING TO DO TO THE REDUCTION COLUMN-ACCUMULATOR MATRIX IF THE ADDED COLUMN HAS ONLY ONE NONZERO
				if  Wcp[addcolind+1]-Wcp[addcolind] != 1
					rangefrom 				=	Wcp[addcolind] : (Wcp[addcolind+1]-2) # WE DON'T NEED THE LAST ENTRY SINCE WE KNOW IT WIL ANIHILATE WITH THE EXISTING STUFF
					rv 						=	Wrv[rangefrom]
					vl 						=	Wvl[rangefrom] .* scalar # scale the coluQn so as to clear the bottom element
					rrv_rcp_rvl_rpa_rpca_rv_vl__appendwithextra!(rrv,rcp,rvl,rpa,rpca,rv,vl)
					# 	NOTE: WE HAVE ALREADY PULLED BACK ALL THE REQUISITE POSTS WHEN THE VALUE OF THE BOTTOM ENTRY WAS CALCULATED; WE DON'T NEED TO DO IT HERE
				end

				# 	ADD TO THE TRANSFORM COLUMN-ACCUMULATOR
				if returntransform
					Trv_Tcp_Tvl_trv_tcp_tvl_tpa_tpca_addcolind_scalar__addcol!(Trv,Tcp,Tvl,trv,tcp,tvl,tpa,tpca,addcolind,scalar) # IT ACTUALLY MAKES SENSE TO HAVE A DEDICATED FUNCTION FOR THIS OPERATION

				end
			end
		end

		if returntransform
			rv,vl 			=	rrv_rcp_rvl_rpa_rpca__colsumIE_rv_vl!(trv,tcp,tvl,tpa,tpca)
			Srv_Scp_Svl_colind_rv_vl__addwithextra!(Trv,Tcp,Tvl,colind,rv,vl)
		end

    end

	# 	STATE: THE REDUCTION IS COMPLETE
	# 	------------------------------------------------------------------------

	# 	ASSEMBLE THE REDUCED MATRIX
	# 	---------------------------

	# 	PUT THE INFO OF WHICH COLUMN IS STORED WHERE IN A HANDY FORMAT
	colind__storedinQa  		=	falses(Qn)
	for rowind in 1:Qm
		if row__pcola[rowind] 	< 	0
			colind__storedinQa[-row__pcola[rowind]] 	=	true
		end
	end

	# 	OBTAIN THE DENSITY PATTERN
	Pcp 						=	ones(Int64, Qn+1)
	for colind in 1:Qn
		if colind__storedinQa[colind]
			Pcp[colind+1] 		=	Pcp[colind] + Qcp[colind+1] - Qcp[colind]
		else
			Pcp[colind+1] 		=	Pcp[colind] + Rcp[colind+1] - Rcp[colind]
		end
	end

	# 	FILL IN THE VALUES
	Prv 						=	Array{Int64}(undef,Pcp[end]-1)
	Pvl 						=	Array{Rational{Int64}}(undef,Pcp[end]-1)
	for colind in 1:Qn
		if colind__storedinQa[colind]
			Pran 				= 	cran(Pcp,colind)
			Qran 				=	cran(Qcp,colind)
			Prv[ Pran ] 		=	Qrv[ Qran ]
			Pvl[ Pran ] 		=	Qvl[ Qran ]
		else
			Pran 				= 	cran(Pcp,colind)
			Nran 				=	cran(Rcp,colind)
			Prv[ Pran ] 		=	Rrv[ Nran ]
			Pvl[ Pran ] 		=	Rvl[ Nran ]
		end
	end


	# 	DEFINE THE RETURN array
	#---------------------------------------------------------------------------
	returndict 					=	Dict(	"Prv"	=> 	Prv,
											"Pcp"	=> 	Pcp,
											"Pvl" 	=> 	Pvl
											)
	# printval(tcp, "tcp")
	# printval(Tcp, "Tcp")
	# printval(trv, "trv")
	# printval(Trv, "Trv")
	# printval(Tvl, "Tvl")
	# printval(tvl, "tvl")
	# 	IF REQUIRED, TRIM THE TRANSFORM
	#------------------------------------------
	if returntransform
		Trv 						= 	Trv[1:Tcp[end]-1]
		Tvl 						= 	Tvl[1:Tcp[end]-1]
		returndict["Trv"] 			=	Trv
		returndict["Tcp"] 			=	Tcp
		returndict["Tvl"] 			=	Tvl
	end

	# 	IF REQUIRED, RETURN THE STORAGE DEPOT FOR THE REDUCED MATRIX
	#------------------------------------------
	if returntempreduced
		returndict["Rrv"] 			=	Rrv
		returndict["Rcp"] 			=	Rcp
		returndict["Rvl"] 			=	Rvl
	end

	# 	IF REQUIRED, COMPUTE PIVOT ELEMENTS
	#------------------------------------------
	if returnpivotindices
		prowa 							=	[x 				for x = 1:Qm  if row__pcola[x]!= 0]
		pcola 							=	[row__pcola[x]	for x = 1:Qm  if row__pcola[x]!= 0]
		returndict["prowa"] 			=	prowa
		returndict["pcola"] 			=	pcola
	end

	return returndict
	# if returntransform
	# 	if returntempreduced
	# 		return Prv,Pcp,Pvl,Trv,Tcp,Tvl,Rrv,Rcp,Rvl
	# 	else
	# 		return Prv,Pcp,Pvl,Trv,Tcp,Tvl
	# 	end
	# else
	# 	if returntempreduced
	# 		return Prv,Pcp,Pvl,Rrv,Rcp,Rvl
	# 	else
	# 		return Prv,Pcp,Pvl
	# 	end
	# end
end

function Trv_Tcp_Tvl_trv_tcp_tvl_tpa_tpca_addcolind_scalar__addcol!(Trv,Tcp,Tvl,trv,tcp,tvl,tpa,tpca,addcolind,scalar)
	"""
	THIS ASSUMES A SQUARE MATRIX WITH SILENT 1'S ON THE DIAGONAL
	"""

	#	COUNT NEW ENTRIES
	tpcOLD 								=	tpca[1] 	# OLD NUMBER OF EXISTING POSTS
	tpcNEW 								=	tpcOLD+1 	# NEW NUMBER OF EXISTING POSTS
	newentryc 							=	Tcp[addcolind+1]-Tcp[addcolind] + 1 # NUMBER OF NEW ENTRIES TO INSERT IN rv AND vl; THE PLUS ONE ACCOUNDS FOR THE FACT THAT THE DIAGONAL ENTRIES ARE "GHOSTED"
	uboundNEW 							=	tcp[tpcNEW] + newentryc

	#	LENGTHEN THE rv, vl ARRAYS IF NEEDED (ACTUALLY WE'RE GIVING OURSELVES 1 EXTRA SPACE, BUT THIS AVOIDS AN EXTRA ARITHMETIC OPERATION)
	if length(trv) <  uboundNEW
		append!(trv,Array{Int64}(			undef,	length(Tcp)))
		append!(tvl,Array{Rational{Int64}}(	undef,	length(Tcp)))
	end

	#	LOG THE DIAGONAL ENTRY (CREATED FROM GHOSTING)
	tcp[tpcNEW+1] 						=	uboundNEW
	trv[uboundNEW-1] 					=	addcolind
	tvl[uboundNEW-1] 					=	scalar

	#	LOG THE ABOVE-DIAGONAL ENTRIES
	fillran_almost::UnitRange{Int64}	=	tcp[tpcNEW] : (uboundNEW - 2) #	THE UPPER BOUND REPRESENTED BY uboundNEW IS **STRICT**, THEREFORE WE PLACE THE LAST ENTRY IN uboundNEW-1 AND THE PENULTIMATE IN uboundNEW-2
	fromrange::UnitRange{Int64}			=	cran(Tcp,addcolind)
	trv[fillran_almost] 				=	Trv[fromrange]
	tvl[fillran_almost] 				=	Tvl[fromrange] * scalar

	# 	ADD A NEW POST TO THE LEDGER
	tpca[1] 							=	tpcNEW
	tpa[tpcNEW] 						=	uboundNEW-1
end

function rrv_rcp_rvl_rpa_rpca_rv_vl__reset!(rrv,rcp,rvl,rpa,rpca,rv,vl;addlength=length(rcp))
	"""
	COPY THE INPUT VECTOR INTO THE ARRAY, AND RESET THE NUMBER OF POSTS TO 1
	"""
	if length(rv) != length(vl)
		print("error in rrv_rcp_rvl_rpa_rpca_rv_vl__reset: rv and vl must have equal length")
		return
	end
	m 					=	length(rv)
	if m > length(rvl)
		addlength 		=	m + addlength
		append!(rvl, Array{Rational{Int64}}(undef,addlength))
		append!(rrv, Array{Int64}(undef,addlength))
	end
	rvl[1:m] 			=	vl
	rrv[1:m] 			=	rv
	rcp[2] 				=	m+1
	rpa[1] 			=	rcp[2]-1
	rpca[1] 			=	1
end

function rrv_rcp_rpa_rpca__locus_max(rrv,rcp,rpa,rpca)
	locus_max 			=	0
	for p = 1:rpca[1]
		if rcp[p] <= rpa[p]
			locus_max 	=	max(locus_max,rrv[rpa[p]]) # lowest nonzero row???
		end
	end
	return locus_max
end

function rrv_rcp_rvl_rpa_rpc_locus__rowsum(rrv,rcp,rvl,rpa,rpc,locus; moveposts=false)
	"""
	RETUNS THE SUM OF THE ENRIES IN ROW locus.  **IT ASSUMES** THAT locus IS THE
	MAXIMUM ELEMENT OF THE SUPPORT OF THE RELEVANT PART OF THE REDUCTION ARRAY
	"""
	rowsum::Rational{Int64} 	=	0
	for pstind 	= 	1:rpc
		pst 	=	rpa[pstind]
		if (pst >= rcp[pstind]) &&  (rrv[rpa[pstind]] 	==	locus) # must check that (a) we haven't already worked through the column completely, and (b) the lowest entry is in the correct row
			rowsum  			+= 	rvl[rpa[pstind]] 	# add to the entry value
			if moveposts
				rpa[pstind]	-= 	1					# move the post back one step
			end
		end
	end
	return rowsum
end

function rrv_rcp_rvl_rpa_rpca_rv_vl__appendwithextra!(rrv,rcp,rvl,rpa,rpca,rv,vl; rvgrowsize=length(rcp))
	# 	INCREASE THE RECORDED NUMBER OF POSTS BY 1
	rpc 				=	rpca[1]+1
	rpca[1]				=	rpc
	# 	TACK NEW PARTS ONTO ARRAYS
	Srv_Scp_Svl_colind_rv_vl__addwithextra!(rrv,rcp,rvl,rpc,rv,vl,rvgrowsize=rvgrowsize)
	# 	POSITION THE NEW POST AT THE BOTTOM-MOST ENTRYOF ITS COLUMN // THIS IS RIGHTLY DONE **AFTER** THE COMPLETING THE COPY, SO AS TO TAKE ADVANTAGE OF THE NEW rcp
	rpa[rpc] 			=	rcp[rpc+1]-1 # this is the location of the last recorded rowval
end

function Srv_Scp_Svl_colind_rv_vl__addwithextra!(Srv,Scp,Svl,colind,rv,vl; rvgrowsize=length(Scp))
	"""
	FILL THE GIVEN VALUES INTO COLUMN colind OF MATRIX S
	"""
	# 	IF NEED BE, GROW THE ARRAYS
	rvnewc 				=	length(rv)
	rvoldc 				=	Scp[colind] -1 # number of entries for the columns that stricely precede colind
	if rvnewc + rvoldc	> length(Srv)
		rvgrowsize 		=	rvnewc + rvoldc - length(Srv) + rvgrowsize
		append!(Srv, Array{Int64}(          undef, rvgrowsize))
		append!(Svl, Array{Rational{Int64}}(undef, rvgrowsize))
	end

	Scp[colind+1] 		=	Scp[colind] + length(rv)
	fillrange 			=	cran(Scp,colind)
	Srv[cran(Scp,colind)]  	=	rv
	Svl[cran(Scp,colind)] 	=	vl
end
function rrv_rcp_rvl_rpa_rpca__colsumIE_rv_vl!(rrv,rcp,rvl,rpa,rpca)
	"""
	SUM THE ENTRIES THAT LIE "ABOVE" THE POSTS
	THIS WILL MOVE ALL POSTS TO THE FRONT
	"""
	# rpa 				=	deepcopy(rpa) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TEMPORARY PLEASE DELETE !!!!!!!!!!!!!!!!!!!!!!!!!
	rpc 				=	rpca[1] # we do this in part to avoid overwriting rpca, since we will manipulate this value below
	locus_max::Int64 	= 	rrv_rcp_rpa_rpca__locus_max(rrv,rcp,rpa,[rpc])
	if locus_max == 0
		rv::Array{Int64}			=	zeros(Int64,0)
		vl::Array{Rational{Int64}}	=	zeros(Rational{Int64},0)
	else
		entries_maxc 	=	rcp[rpc+1]-1
		rv 				=	zeros(Int64,entries_maxc)
		vl 				=	zeros(Rational{Int64},entries_maxc)
		nzentryc 		=	0
		for rowiteration 		= 	1:entries_maxc
			rowsum 				=	rrv_rcp_rvl_rpa_rpc_locus__rowsum(rrv,rcp,rvl,rpa,rpc,locus_max; moveposts=true) # THIS SHOULD RETURN ZERO IF ALL THE POSTS ARE IN THEIR FINAL POSITIONS
			if rowsum	!=	0
				nzentryc 		+=	1
				rv[end-nzentryc+1]	=	locus_max # we do this bc we want to fill values in order back to front, high to low
				vl[end-nzentryc+1] 	=	rowsum  # we do this bc we want to fill values in order back to front, high to low
			end

			locus_max 	=	rrv_rcp_rpa_rpca__locus_max(rrv,rcp,rpa,[rpc])
		end
		rv 				=	rv[end-nzentryc+1:end] # we do this bc we want to fill values in order back to front, high to low
		vl 				=	vl[end-nzentryc+1:end] # we do this bc we want to fill values in order back to front, high to low
	end
	return rv, vl
end

function printval(var,varname)
	println(string(varname," = $(var)"))
end

function complexrank(C;dim=1)
	sd 		= 	dim+1
	if 	dim > C["input"]["maxdim"]+1 || dim < 0
		return 0
	elseif C["input"]["model"] == "complex"
		return length(C["cp"][sd])-1
	else
		return length(C["farfaces"][sd])
	end
end

function stackedsubmatrices(
	Mrowval,#::Array{Tv,1},
	Mcolptr,#::Array{Tv,1},
	rows1,#::Array{Tv,1},
	rows2,#::Array{Tv,1},
	cols,#::Array{Tv,1},
	Mm::Tv)  where Tv<:Integer

	n = length(cols)
	suppcol1 = falses(Mm)
	suppcol2 = falses(Mm)
	suppcol1[rows1].=true
	suppcol2[rows2].=true
	nz1 = 0; nz2 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,cols[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
			elseif suppcol2[i]
				nz2+=1
			end
		end
	end
	rv1 = Array{Tv}(undef,nz1)
	rv2 = Array{Tv}(undef,nz2)
	cp1 = Array{Tv}(undef,n+1); cp1[1]=1
	cp2 = Array{Tv}(undef,n+1); cp2[1]=1
	nz1 = 0; nz2 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,cols[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
				rv1[nz1]=i
			elseif suppcol2[i]
				nz2+=1
				rv2[nz2]=i
			end
		end
		cp1[jp+1] = nz1+1
		cp2[jp+1] = nz2+1
	end
	return rv1,rv2,cp1,cp2
end

function buildclosefromfar(farfaces,firstv,sd)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	# destinationmatrix = Array{Int64}(undef,sd,n)
	if sd == 1
		return Array{Int64}(undef,0,m)
	end
	lclosefaces = Array{Int64}(undef,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if sd == 2
		return lclosefaces'
	end
	for i = 3:sd
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	return lclosefaces
end
function buildallfromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr,selectedcolumnindices;verbose=false)
	if verbose
		println("PLEASE NOTE: COLUMNS MUST BE IN SORTED ORDER FOR THIS TO WORK PROPERLY")
	end
	m = length(hcolptr)-1
	numhigs = length(hrowval)
	numselected = length(selectedcolumnindices)
	rowdepth = size(lclosefaces,1)
	sd = rowdepth+1
	hclosefaces = Array{Int64}(undef,sd+1,numselected)
	if numselected == 0
		return hclosefaces
	end
	rosettacol = Array{Int64}(undef,maximum(lrowval))
	columnsupp = falses(numhigs)
	columnsupp[selectedcolumnindices].=true
	columnmarker = 0
	for i = 1:m
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			if columnsupp[j]
				columnmarker+=1
				farface = hrowval[j]
				for k = 1:rowdepth
					hclosefaces[k,columnmarker]=rosettacol[lclosefaces[k,farface]]
				end
				hclosefaces[sd,columnmarker] = rosettacol[lrowval[farface]]
				hclosefaces[sd+1,columnmarker] = farface
			end
		end
	end
	return hclosefaces
end

function buildclosefromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr;facecard = size(lclosefaces,1)+1)
	m = length(hcolptr)-1
	n = length(hrowval)
	hclosefaces = Array{Int64}(undef,facecard,n)
	if n == 0
		return hclosefaces
	else
		rowdepth = facecard-1
		rosettacol = Array{Int64}(undef,maximum(lrowval))
		for i = 1:m
			rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
			for j = cran(hcolptr,i)
				farface = hrowval[j]
				for k = 1:rowdepth
					hclosefaces[k,j]=rosettacol[lclosefaces[k,farface]]
				end
				hclosefaces[facecard,j] = rosettacol[lrowval[farface]]
			end
		end
		return hclosefaces
	end
end


function buildclosefromfar(farfaces,firstv,sd,columnsinorder)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	if sd == 1
		return Array{Int64}(undef,0,m)
	end
	lclosefaces = Array{Int64}(undef,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)].=i
	end
	if sd == 2
		return lclosefaces[columnsinorder]'
	end
	for i = 3:(sd-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	lclosefaces = buildclosefromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;facecard = sd-1)
	return lclosefaces
end

function empteval(f,a,c)
	if isempty(a)
		return c
	else
		return f(a)
	end
end

function ff2aflight_sc2(farfaces,firstv,columns)
	sd = 2
	if isempty(farfaces[sd])
		return Array{Int64}(undef,2,0)
	end
	f0faces::Array{Int64,1} = farfaces[sd]
	colptr::Array{Int64,1} = firstv[2]
	columnpost::Int64   = 1
	columnpostp1::Int64 = 2
	faces::Array{Int64,2} = Array{Int64}(undef,2,length(columns))

	for fp = 1:length(columns)
		f0 = columns[fp]
		if f0 >= colptr[columnpostp1]
			while f0 >= colptr[columnpostp1]
				columnpostp1+=1
			end
			columnpost = columnpostp1-1
		elseif f0 < colptr[columnpost]
			while f0 < colptr[columnpost]
				columnpost-=1
			end
			columnpostp1 = columnpost+1
		end
		faces[1,fp] = columnpost
		faces[2,fp] = f0faces[f0]
	end
	return faces
end

function ff2aflight_sc3(farfaces,firstv,columns)
	sd = 3

	if isempty(farfaces[sd])
		return Array{Int64}(undef,3,0)
	end

	fcfaces::Array{Int64,2} = buildclosefromfar(farfaces,firstv,sd-1,1:length(farfaces[2]))

	f0faces::Array{Int64,1} = farfaces[sd]
	f1faces::Array{Int64,1} = farfaces[sd-1]

	fvscm0::Array{Int64,1}  = firstv[sd]
	fvscm1::Array{Int64,1}  = firstv[sd-1]
	fvscm2::Array{Int64,1}  = firstv[sd-2]

	holdi=[1];holdip1=[2]
	t1::Array{Int64,1} = Array{Int64}(undef,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)

	faces::Array{Int64,2} = Array{Int64}(undef,3,length(columns))
	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		f3 = fcfaces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1} ,holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		faces[1,fp] = t1[f3]
		faces[2,fp] = t1[f2]
		faces[3,fp] = f1
	end
	return faces
end

function ff2aflight_scgt3(farfaces,firstv,sd,columns)

	if isempty(farfaces[sd])
		return Array{Int64}(undef,sd,0)
	end

	f0faces::Array{Int64,1} = farfaces[sd]
	f1faces::Array{Int64,1} = farfaces[sd-1]
	f2faces::Array{Int64,1} = farfaces[sd-2]
	fcfaces::Array{Int64,2} = buildallfromfar(farfaces,firstv,sd-2,1:(firstv[sd-2][end]-1))

	fvscm0::Array{Int64,1}  = firstv[sd]
	fvscm1::Array{Int64,1}  = firstv[sd-1]
	fvscm2::Array{Int64,1}  = firstv[sd-2]
	fvscm3::Array{Int64,1}  = firstv[sd-3]

	holdi=[1];holdip1=[2];holdj=[1];holdjp1=[2]
	t1::Array{Int64,1} = Array{Int64}(undef,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)
	t2::Array{Int64,1} = Array{Int64}(undef,fvscm3[end]-1);t2[crows(fvscm2,f2faces,1)]=cran(fvscm2,1)

	scm0::Int64 = sd; scm1::Int64 = sd-1; scm2::Int64 = sd-2
	faces::Array{Int64,2} = Array{Int64}(undef,sd,length(columns))

	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1},holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		updatetranslator!(f1::Int64,fvscm1::Array{Int64,1},holdj::Array{Int64,1},holdjp1::Array{Int64,1},t2::Array{Int64,1},fvscm2::Array{Int64,1},f2faces::Array{Int64,1})
		for i = 1:scm2
			faces[i,fp] = t1[t2[fcfaces[i,f2]]]
		end
		faces[scm1,fp] = t1[f2]
		faces[scm0,fp] = f1
	end
	return faces
end

function ff2aflight(farfaces,firstv,sd,columns)
	if sd == 1
		return Array{Int64}(undef,0,length(columns))
	elseif sd == 2
		return ff2aflight_sc2(farfaces,firstv,columns)
	elseif sd == 3
		return ff2aflight_sc3(farfaces,firstv,columns)
	else
		return ff2aflight_scgt3(farfaces,firstv,sd,columns)
	end
end

function ff2aflight(D::Dict,sd,columns)
	farfaces = D["farfaces"]; firstv = D["firstv"]
	faces = ff2aflight(farfaces,firstv,sd,columns)
	return faces
end

function crows(A::SparseMatrixCSC,j)
	return A.rowval[cran(A,j)]
end

function crows(colptr::Array,rowval::Array,j)
	return rowval[cran(colptr,j)]
end

function boundarymatrix(C;dim=1,rows="a",cols="a")
	"""
	Returns a sparse array of size (# dim-1 cells) x (# dim cells)
	"""
	crr 					= 	complexrank(C,dim=dim-1)
	crc 					= 	complexrank(C,dim=dim)
	if rows == "a"
		rows 				= 	1:crr#complexrank(C,dim=dim-1)
	end
	if cols == "a"
		cols 				= 	1:crc#complexrank(C,dim=dim)
	end
	if empteval(maximum,cols,0) > crc
		println()
		println("error: keyword argument <cols> contains an integer greater than the rank of the complex in dimension <dim>")
		println()
		return
	elseif empteval(maximum,rows,0) > crr
		println()
		print("error: keyword argument <rows> contains an integer greater than the rank of the complex in dimension <dim-1>")
		println()
		return
	end
	if isempty(rows) || isempty(cols)
		rv 					= 	zeros(Int64,0)
		cp 					= 	ones(Int64,length(cols)+1)
		return 					rv,cp
	end
	ncols 					= 	length(cols)
	nrows 					= 	length(rows)
	sd 						= 	dim+1;
	if haskey(C,"farfaces")
		rv 					= 	ff2aflight(C,dim+1,cols)
		rv  				= 	reshape(rv,length(rv))
		cp  				= 	convert(Array{Int64,1},1:sd:(ncols*sd+1))
		cols 				= 	1:ncols
	else
		rv 					= 	C["rv"][sd]
		cp 					= 	C["cp"][sd]
	end
	rv,dummrv,cp,dummycp 	= 	stackedsubmatrices(
								rv,
								cp,
								rows,
								Array{Int64}(undef,0),
								cols,
								max(empteval(maximum,rows,0),empteval(maximum,rv,0))
								)
	return 						rv,cp
end

function vertexrealization(farfaces,firstv,facecardinality,facenames)
	numfaces::Int64 = length(facenames)
	fc::Int64 = facecardinality
	m::Int64 = length(firstv[2])-1
	preallocationspace = 0
	loci::Array{Int64,1} = copy(facenames)
	vrealization = Array{Int64}(undef,facecardinality,numfaces)
	post0::Int64 = 1
	post1::Int64 = 1

	for sd = facecardinality:-1:1
		cp::Array{Int64,1} = firstv[sd]
		for i = 1:numfaces
			locus = loci[i]
			if cp[post0] > locus
				while cp[post0] > locus
					post0-=1
				end
				post1 = post0+1
			elseif cp[post1] <= locus
				while cp[post1] <= locus
					post1+=1
				end
				post0 = post1-1
			end
			loci[i] = farfaces[sd][locus]
			vrealization[sd,i]=post0
		end
	end
	return vrealization
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict,facecardinality,facenames)
	return vertexrealization(D["farfaces"],D["firstv"],facecardinality,facenames)
end
function updatetranslator!(f0,firstv0,holdi,holdip1,t,firstv1,farfaces1)
	if firstv0[holdip1[1]] <= f0
		while firstv0[holdip1[1]]<= f0
			holdip1[1]+=1
		end
		holdi[1] = holdip1[1]-1
		t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
	elseif firstv0[holdi[1]] > f0
		while firstv0[holdi[1]] > f0
			holdi[1]-=1
		end
		holdip1[1] = holdi[1]+1
		t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
	end
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict;dim = 1, class = 1)
	sd = dim+2
	facecard = dim+1

	if haskey(D,"cyclerep")
		rep = D["cyclerep"][sd][class]
	else
		cyclename = barname2cyclename(D,class;dim = dim)
		rep = getcycle(D,sd,cyclename)
	end

	vrealization = vertexrealization(D::Dict,facecard,rep)
end

function eirened_dim__signedboundarymatrix(C; dim=1, useoriginalvertexlabels=true, returnvertexrealization=true)

	# 	OBTAIN THE NONZERO PATTERN OF THE BOUNDARY MATRIX
	brv,bcp 				=	boundarymatrix(C,dim=dim)
	bvl 					=	zeros(Rational{Int64},length(brv))

	# 	COUNT THE NUMBER OF CELLS IN DIMENSION dim AND dim-1
	nhig 					=	length(C["farfaces"][dim+1])
	nlow 					=	length(C["farfaces"][dim+0])

	# 	ENUMERATE THE VERTICES IN EACH HIGH/LOW DIMENSIONAL CELL
	higverts 				=	vertexrealization(C,dim+1,1:nhig)
	lowverts 				=	vertexrealization(C,dim+0,1:nlow)

	# 	IF DESIRED, REINDEX THE VERTICES OF THE COMPLEX`
	if useoriginalvertexlabels
		lowverts 			=	C["nvl2ovl"][lowverts]
		lowverts 			=	sort(lowverts, dims=1) 	# vertices have to appear in sorted order for the algorithm to work
		higverts 			=	C["nvl2ovl"][higverts]
		higverts 			=	sort(higverts, dims=1)	# vertices have to appear in sorted order for the algorithm to work
	end

	# 	OBTAIN SIGNS FOR THE NONZERO ENTRIES
	bvlind 					=	0
	for celhig = 1:nhig  								# celhig := a cell of the higher dimension
		vtahig 				=	higverts[:,celhig] 		# vertices of the higher cell
		cellowa 			=	crows(bcp,brv,celhig) 	# nonzero rows of the column indexed by celhig
		for cellow = cellowa 							# cellow := a cell of the lower dimension (and a face of the higher cell)
			vtalow 			=	lowverts[:,cellow]		# vertices of the lower cell
			p 				=	0
			for slot = 1:dim
				if vtalow[slot] != vtahig[slot]
					p 		=	slot
					break
				end
			end
			if p == 0
				p = dim+1 # recall that the high dim cell has dim+1 vertices
			end
			if p % 2 == 0
				coefficient 	=	-1 		# if the missing vertex is in an even location, give a negative entry
			else
				coefficient 	=	1		# otherwise, give a positive
			end

			bvlind 			=	bvlind+1
			bvl[bvlind] 	=	coefficient

		end
	end

	# 	RECOVER THE SIZE OF THE MATRIX
	bm 						=	nlow
	bn 						=	nhig

	# 	RETURN
	return brv,bcp,bvl,bm,bn,lowverts,higverts
end


function buildallfromfar(farfaces,firstv,sd,columnsinorder;verbose = false)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	# destinationmatrix = Array{Int64}(undef,sd,n)
	if sd == 1
		return Array{Int64}(undef,0,m)
	end
	lclosefaces = Array{Int64}(undef,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)].=i
	end
	if sd == 2
		return vcat(lclosefaces[columnsinorder]',farfaces[sd][columnsinorder]')
	end
	for i = 3:(sd-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
		#gc()
	end
	lclosefaces = buildallfromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;verbose = verbose)
	#gc()
	return lclosefaces
end

function buildallfromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr,selectedcolumnindices;verbose=false)
	if verbose
		println("PLEASE NOTE: COLUMNS MUST BE IN SORTED ORDER FOR THIS TO WORK PROPERLY")
	end
	m = length(hcolptr)-1
	numhigs = length(hrowval)
	numselected = length(selectedcolumnindices)
	rowdepth = size(lclosefaces,1)
	sd = rowdepth+1
	hclosefaces = Array{Int64}(undef,sd+1,numselected)
	if numselected == 0
		return hclosefaces
	end
	rosettacol = Array{Int64}(undef,maximum(lrowval))
	columnsupp = falses(numhigs)
	columnsupp[selectedcolumnindices].=true
	columnmarker = 0
	for i = 1:m
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			if columnsupp[j]
				columnmarker+=1
				farface = hrowval[j]
				for k = 1:rowdepth
					hclosefaces[k,columnmarker]=rosettacol[lclosefaces[k,farface]]
				end
				hclosefaces[sd,columnmarker] = rosettacol[lrowval[farface]]
				hclosefaces[sd+1,columnmarker] = farface
			end
		end
	end
	return hclosefaces
end
