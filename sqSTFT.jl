    # Synchrosqueezing modifed from tfrrsp.m by Hau-tieng Wu 2013
    #
    #	[TFR,RTFR] = sqSTFT(X, T, N, H, DH, TRACE)
    #	computes the STFT and its SST
    # INPUTS:
    #	 X     : analysed signal.
    #	 T     : the time instants or rather indices in the form of t0:dt:t1
    #             (default : 1:length(X), and dt controls the hop size).
    #	 N     : number of frequency bins (default : length(X)).
    #	 H     : frequency smoothing window, H(0) being forced to 1.
    #  DH	   : DH (derivative of H)
    #	 TRACE : if true, the progression of the algorithm is shown (default : false).
    # OUTPUTS:
    #	 TFR   : STFT
    #	 RTFR  : reassigned version. When called without output arguments,
    #	        TFRRSP runs TFRQVIEW.
    #	Example :
    #	 sig=fmlin(128,0.1,0.4); t=1:2:128
    #	 h=tftb_window(17,"Kaiser"); tfrrsp(sig,t,64,h,1)
    #
    #	See also all the time-frequency representations listed in
    #	 the file CONTENTS (TFR*)
    
    #	F. Auger, May-July 1994, July 1995.
    #	Copyright (c) 1996 by CNRS (France).
    #
    #	------------------- CONFIDENTIAL PROGRAM --------------------
    #	This program can not be used without the authorization of its
    #	author(s). For any comment or bug report, please send e-mail to
    #	f.auger@ieee.org
    #
    # Translated into julia and modified by Naoki Saito, Dec. 10, 2016

import Pkg
Pkg.add("TickTock")
Pkg.add("JLD")
Pkg.add("SpecialFunctions")
Pkg.add("Plots")
Pkg.add("FFTW")
Pkg.add("Statistics")
Pkg.add("LinearAlgebra")
using TickTock
using FFTW
using LinearAlgebra
using Statistics
using JLD
using SpecialFunctions
using Plots
include("hermf.jl")
include("nround.jl")

function sqSTFT(x, t, N, h, Dh, trace = false)

	# If x is a matrix, it assigns the number of rows to xrow and the number of 	    columns to xcol.

	# If x is a vector, it assigns the length of x to xrow and sets xcol to 1.
	
	if ndims(x) != 1 
		
		xrow, xcol = size(x) 
		
	else
		
		xrow = length(x); xcol = 1
		
	end

	# Finding the number of rows in the time instant variable
	trow=length(size(t))

	# Checking validity of the input signal and time instant
	if xcol !=1
		
		error("X must have only one column!")
		
	elseif trow !=1
		
		error("T must only have one row!")
	
	elseif nextpow(2,N) != N && trace
		
		println("For a faster computation, N should be a power of two.")
		
	end
	

	# Calculating half-window size for FFT 
	N₂ = Int(floor(N/2))
    tcol = length(t)

	# Initiating the size of the smoothing window
	hrow, hcol = size(h)

	# Checking validity of the smoothing window
	# Ensures it only has 1 column and odd number of rows
	if hcol!=1 || rem(hrow,2)==0
		
        error("H must be a smoothing window with odd length!")
		
    end
	
	# Calculating half-length of the smoothing window
	Lh = Int((hrow-1)/2) 

	# Defining the time step based on the time instant
	
	if tcol==1
		
    	Dt=1
		
	else
		
		Δt=diff(t)
		Mini = minimum(Δt); Maxi = maximum(Δt)
	# Limiting the time step within an epsilon distance
		if Maxi-Mini > eps()
			
			error("The time instants must be regularly sampled!")
		else
			
			Dt=Mini
			
		end
	end
	
	# tfr and tf3 will be used to store intermediate results 
	# while calculating SFTF and FFT
	
	tfr = zeros(N,tcol) 
	tf3 = zeros(N,tcol)	

	
	if trace
		
		println("Spectrogram and its differentiation.")
		
	end
	
	# Measuring the elapsed time
	if trace
		
		tick()
		
	end
	

	# Computing STFT and its differentiation

	#Iterate over each column 
	for icol = 1:tcol
		
	    tᵢ = t[icol]
		
	# For each time instant, calculate the corresponding indice
	# with the in window size 'N' 		
		
    	τ = -minimum([nround(N/2)-1,Lh,tᵢ-1]):minimum([nround(N/2)-1,Lh,xrow-tᵢ])	
    	indices= rem.((N.+τ),N).+1
		norm_h=norm(h[Lh.+τ.+1])
	# Calculating the normalized STFT and normalized differentiation value
	# by multipying the input signal with the conjugate of the smoothing window and 
	# differentiation respectively
	# and dividing by the normalization factor
    	tfr[indices,icol] = x[tᵢ.+τ] .* conj(h[Lh.+τ.+1]) / norm_h
    	tf3[indices,icol] = x[tᵢ.+τ] .* conj(Dh[Lh.+τ.+1]) / norm_h

		
	end

	# Converting the time-domain representation of the signal 
	# to the frequency-domain representation by applying FFT to tfr and tf3
	
	tfr = fft(tfr,1)
	tf3 = fft(tf3,1)

	tfr_copy=deepcopy(tfr)
	tf3_copy=deepcopy(tf3)
	
	tfr_vec=vec(tfr_copy)
	tf3_vec=vec(tf3_copy)

	# Finding thelocations in the frequency spectrum
	# where the STFT values are non-zero
	
	avoid_warn = findall(!iszero,tfr_vec)

	tf3[avoid_warn] = round.(imag(N*tf3[avoid_warn]./tfr[avoid_warn]/(2.0*pi)),RoundNearestTiesUp)

	# Display the elapsed time
	if trace 
		
		QQ = tock(); println("Elapsed time = " * string(QQ)) 
		
	end
	
	if trace 
		
		println("Synchrosqueezing: ") 
		
	end

	# Creating a complex matrix to store the result of the Synchrosqueezing transform
	rtfr = zeros(Complex{Float64},N₂,tcol)

	# Calculating the mean energy of the signal
	Ex = mean(abs.(x[minimum(t):maximum(t)]).^2)

	# Defining a threshold for the SST based on the mean energy of the signal
	Threshold = (1.0e-8)*Ex

	if trace
		
		tick()
		
	end
	
	# Calculating the SST
	for icol = 1:tcol
		
		for jcol = 1:N
	# For each element in tfr that exceeds the threshold value
	# calculate the shift index
			if abs(tfr[jcol,icol]) > Threshold
				
				jcolhat = jcol - Int64(real(tf3[jcol,icol]))
	# Ensuring that jcolhat wraps around to the valid range
				jcolhat = rem(rem(jcolhat-1,N)+N,N)+1
	# If the shift index is within a valid range		
	# the corresponding element in rtfr is incremented by the value in tfr
				if jcolhat < N₂+1
					
					rtfr[jcolhat,icol] = rtfr[jcolhat,icol] + tfr[jcol,icol]
					
				end
			end
		end
	end

	if trace
		RR=tok(); println("Elapsed time = " * string(RR))
	end

	return tfr, rtfr
end


#Test Example:

#H, DH = hermf(301,1,6);
#lpfs=load("lpfs.jld","lpfs")
#tmp=lpfs[:,34];

#tmp2=zeros(1024);
#tmp2[1:1000]=tmp
#tmp2[1001:end].=tmp[end]
#tmp3=(tmp2.-mean(tmp2));
#tmp3=tmp3./norm(tmp3)

#tfr, rtfr=sqSTFT(tmp3, 1:1024, 1024, H', DH', true)

#f1=heatmap(abs.(rtfr[1:25,:]),c=:viridis)

#f2=plot(lpfs[:,34],legend=false)

#f3=heatmap(abs.(tfr[1:25,:]),c=:viridis)

#output=plot(f1,f2,f3,layout=grid(3,1,heights=[0.35,0.35,0.3]))
#savefig(output,"test.png")"

